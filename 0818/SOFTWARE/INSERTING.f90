module mod_INSERTING
! **********************************************************************
!     Author : M.YOUAN 
!     Date of last modification : 24.4.02
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
      subroutine new_i_matrix (nmb,Mb,ierr)
      implicit none
      integer ,intent(out)::nmb,ierr
      integer ,intent(in)::Mb
      type(i_matrix ),pointer::matrix_insert
      if (.NOT. iins_init ) then
         nullify(i_matrix_start )
         iins_init  = .TRUE.
      end if
      if (.not.associated(i_matrix_start )) then
         allocate(i_matrix_start ,STAT=ierr)
         i_matrix_start %number= ISP_MATRIX     
         i_matrix_start %number=- i_matrix_start %number    
         nullify(i_matrix_start %pntr)      
      else
         allocate(matrix_insert,STAT=ierr)
         matrix_insert%number= i_matrix_start %number-no_of_types
         matrix_insert%pntr=> i_matrix_start 
         i_matrix_start => matrix_insert    
      end if
      i_matrix_start %DIM=0
      i_matrix_start %property=blas_general+blas_one_base+blas_col_major
      i_matrix_start %new = 1    !new=0:blas_open_handle, new=1: blas_new_handle
      i_matrix_start %format=''
      nullify(i_matrix_start %sub_rows,i_matrix_start %sub_cols)
      nullify(i_matrix_start % i_element_start )
      allocate(i_matrix_start %trb(Mb),i_matrix_start %tre(Mb)) 
      nmb= i_matrix_start %number
      ierr=0
      end subroutine new_i_matrix       
!*
      subroutine dealloc_i_matrix  (nmb,ierr)
      implicit none
      integer ,intent(in)::nmb
      integer ,intent(out)::ierr
      type(i_matrix ),pointer ::matrix_precedent,matrix_tester
      ierr=-1
      if(.not.associated(i_matrix_start %pntr)) then
         if(i_matrix_start %number.eq.nmb) then
            deallocate(i_matrix_start %tre,i_matrix_start %trb,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            deallocate(i_matrix_start ,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            nullify(i_matrix_start )
            ierr=0
            return
         end if
      else
         matrix_tester=> i_matrix_start 
         if(matrix_tester%number.eq.nmb) then
            i_matrix_start =>matrix_tester%pntr
            deallocate(matrix_tester%tre,matrix_tester%trb,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            deallocate(matrix_tester,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            ierr=0
            return
         endif
         matrix_precedent=> i_matrix_start 
         matrix_tester=> i_matrix_start %pntr
         do while((associated(matrix_tester)))
            if(matrix_tester%number.eq.nmb) then
               matrix_precedent%pntr=>matrix_tester%pntr
               deallocate(matrix_tester%tre,matrix_tester%trb,STAT=ierr)
               if (ierr.ne.0) then
                  ierr = blas_error_memdeloc
                  return
               end if
               deallocate(matrix_tester,STAT=ierr)
               if (ierr.ne.0) then
                  ierr = blas_error_memdeloc
                  return
               end if
               ierr=0
               return
            else
               matrix_precedent=>matrix_tester
               matrix_tester=>matrix_tester%pntr
            end if
         end do
      end if
      end subroutine dealloc_i_matrix 
!*
      subroutine iaccess_matrix (pmatrix,nmb,istat)
      implicit none
      type(i_matrix ),pointer ::pmatrix
      integer,intent(out)::istat
      integer,intent(in) ::nmb
      type(i_matrix ),pointer ::matrix_tester
      istat=-1
      matrix_tester=> i_matrix_start        
      do while((matrix_tester%number.ne.nmb).and.&
               (associated(matrix_tester%pntr)))
         matrix_tester => matrix_tester%pntr
      end do
      if (matrix_tester%number.eq.nmb) then
         pmatrix => matrix_tester
         istat = 0 
         return
      else
         nullify(pmatrix)
         istat = blas_error_param
      end if
      end  subroutine iaccess_matrix 
!*
      subroutine new_i_element (pmatrix,nmb_element,istat)
      implicit none
      type(i_matrix ),pointer::pmatrix
      integer,intent(out)::nmb_element,istat
      type(i_element ),pointer::element_insert
      integer :: ierr
      istat = -1
      if (.not.associated(pmatrix% i_element_start )) then
         pmatrix%new=0 !status changed to blas_open_handle     
         allocate(pmatrix% i_element_start ,STAT=ierr)
         if (ierr.ne.0) then
            istat = blas_error_memalloc
            return
         end if
         pmatrix% i_element_start %number=1 !will certainly changed
         nullify(pmatrix% i_element_start %pntr) 
      else
         allocate(element_insert,STAT=ierr)
         if (ierr.ne.0) then
            istat = blas_error_memalloc
            return
         end if
         element_insert%pntr=>pmatrix% i_element_start 
         element_insert%number=pmatrix% i_element_start %number+1       
         pmatrix% i_element_start => element_insert 
      end if
      select case(pmatrix%format)
         case('normal')
            pmatrix% i_element_start %contents%pntin%value=0
            pmatrix% i_element_start %contents%pntin%row_ind=-1
            pmatrix% i_element_start %contents%pntin%col_ind=-1 
            nullify(pmatrix% i_element_start %contents%blin%value)
            nullify(pmatrix% i_element_start %contents%vblin%value)
         case('block')
            nullify(pmatrix% i_element_start %contents%blin%value)
            nullify(pmatrix% i_element_start %contents%vblin%value)
            pmatrix% i_element_start %contents%blin%row_block_ind=-1
            pmatrix% i_element_start %contents%blin%col_block_ind=-1
         case('vblock')
            nullify(pmatrix% i_element_start %contents%blin%value)
            nullify(pmatrix% i_element_start %contents%vblin%value)
            pmatrix% i_element_start %contents%vblin%row_vblock_ind=-1
            pmatrix% i_element_start %contents%vblin%col_vblock_ind=-1
         case default 
            istat = blas_error_param
            return
      end select
      nmb_element=pmatrix% i_element_start %number
      istat=0
      end subroutine new_i_element  
!*
      subroutine dealloc_i_element  (nmb_element,pmatrix,istat)
      implicit none
      integer ,intent(in)::nmb_element
      type(i_matrix ),pointer::pmatrix
      integer ,intent(out)::istat
      type(i_element ),pointer ::element_tester
      integer::ierr
      istat=-1
      if(.not.associated( pmatrix% i_element_start %pntr)) then
         if(pmatrix% i_element_start %number.eq.nmb_element) then
         if(associated(pmatrix% i_element_start %contents%vblin%value))&
             then
             deallocate(pmatrix% i_element_start %contents%vblin%value,&
                        STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end  if
          if(associated(pmatrix% i_element_start %contents%blin%value))& 
            then
            deallocate(pmatrix% i_element_start %contents%blin%value,&
                       STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if
            if (associated(pmatrix% i_element_start )) then
               deallocate(pmatrix% i_element_start ,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if 
            nullify(pmatrix% i_element_start )
         end if      
         istat = 0 
         return
      else 
         element_tester=>pmatrix% i_element_start 
         if(element_tester%number.eq.nmb_element) then
            pmatrix% i_element_start =>element_tester%pntr
            if(associated(element_tester%contents%vblin%value)) then
             deallocate(element_tester%contents%vblin%value,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if
            if(associated(element_tester%contents%blin%value)) then
             deallocate(element_tester%contents%blin%value,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end  if
            if (associated(element_tester)) then
               deallocate(element_tester,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if
            istat=0
            return
         endif
         element_tester=>pmatrix% i_element_start %pntr
         do while((associated(element_tester)))
            if(element_tester%number.eq.nmb_element) then
               if(associated(element_tester%contents%vblin%value)) then
                  deallocate(element_tester%contents%vblin%value,&
                             STAT=ierr)
                  if (ierr.ne.0) then
                     istat = blas_error_memdeloc
                     return
                  end if
               end if
               if(associated(element_tester%contents%blin%value)) then
                  deallocate(element_tester%contents%blin%value,&
                             STAT=ierr)
                  if (ierr.ne.0) then
                     istat = blas_error_memdeloc
                     return
                  end if
               end  if
               if (associated(element_tester)) then
                  deallocate(element_tester,STAT=ierr)
                  if (ierr.ne.0) then
                     istat = blas_error_memdeloc
                     return
                  end if
               end if
               istat=0
               return
            else
               element_tester=>element_tester%pntr
            end if
         end do
      end if
      end subroutine dealloc_i_element   
!*
      subroutine iaccess_element (pelement,nmb_element,&
                                 pmatrix,istat)
      implicit none
      type(i_inelement ),pointer::pelement
      integer,intent(in) ::nmb_element
      type(i_matrix ),pointer::pmatrix
      integer,intent(out)::istat
      type(i_element ),pointer ::element_tester
      istat=-1
      element_tester=>pmatrix% i_element_start        
      do while((element_tester%number.ne.nmb_element)&
               .and.(associated(element_tester%pntr)))
         element_tester => element_tester%pntr
      end do
      if (element_tester%number.eq.nmb_element) then
         pelement => element_tester%contents
         istat = 0 
         return
      else
         nullify(pelement)
         istat = blas_error_param
         return
      end if  
      end  subroutine  iaccess_element 
!*
      subroutine i_element_num  (nmb_element,pmatrix,i,j,istat)
      implicit none
      integer,intent(out),target::nmb_element
      type(i_matrix ),pointer::pmatrix
      integer ,intent(in)::i,j
      integer,intent(out)::istat
      type(i_element ),pointer ::element_tester
      logical:: finder
      istat = -1
      select case(pmatrix%format)
      case('normal')
         element_tester=>pmatrix% i_element_start  
         if(.not.associated( element_tester%pntr)) then 
            if((element_tester%contents%pntin%row_ind.eq.i)&
            .and.(element_tester%contents%pntin%col_ind.eq.j))then
               nmb_element=element_tester%number        
            else
               nmb_element=0 
            end if
         else  
            finder=.false.
            do while((associated(element_tester%pntr))&
                     .and.(.not.finder))
               if((element_tester%contents%pntin%row_ind.eq.i)&
              .and.(element_tester%contents%pntin%col_ind.eq.j))then 
                  finder=.true.
               else
                  element_tester => element_tester%pntr
               end if
            end do
            if((element_tester%contents%pntin%row_ind.eq.i)&
            .and.(element_tester%contents%pntin%col_ind.eq.j)) then   
               nmb_element=element_tester%number      
            else
               nmb_element=0 
            end if
         end if
      case('block')   
         element_tester=>pmatrix% i_element_start  
         if(.not.associated( element_tester%pntr)) then 
            if((element_tester%contents%blin%row_block_ind.eq.i)&
           .and.(element_tester%contents%blin%col_block_ind.eq.j))then 
               nmb_element=element_tester%number        
            else
               nmb_element=0 
            end if
         else  
            finder=.false.
            do while((associated(element_tester%pntr))&
                     .and.(.not.finder))
               if((element_tester%contents%blin%row_block_ind.eq.i)&
           .and.(element_tester%contents%blin%col_block_ind.eq.j))then 
                  finder=.true.
               else
                  element_tester => element_tester%pntr
               end if
            end do
            if((element_tester%contents%blin%row_block_ind.eq.i)&
          .and.(element_tester%contents%blin%col_block_ind.eq.j))then   
               nmb_element=element_tester%number      
            else
               nmb_element=0 
            end if
         end if
      case('vblock')
         element_tester=>pmatrix% i_element_start  
         if(.not.associated( element_tester%pntr)) then 
            if((element_tester%contents%vblin%row_vblock_ind.eq.i)&
         .and.(element_tester%contents%vblin%col_vblock_ind.eq.j))then 
               nmb_element=element_tester%number        
            else
               nmb_element=0 
            end if
         else  
            finder=.false.
            do while((associated(element_tester%pntr)).and.&
                     (.not.finder))
            if((element_tester%contents%vblin%row_vblock_ind.eq.i)&
         .and.(element_tester%contents%vblin%col_vblock_ind.eq.j))then 
                  finder=.true.
               else
                  element_tester => element_tester%pntr
               end if
            end do
            if((element_tester%contents%vblin%row_vblock_ind.eq.i)&
         .and.(element_tester%contents%vblin%col_vblock_ind.eq.j))then   
               nmb_element=element_tester%number      
            else
               nmb_element=0 
            end if
         end if
      case default 
         istat = blas_error_param
         return
      end select
      istat = 0
      end subroutine i_element_num 
!*
      subroutine i_dealloc (nmb,istat)
      implicit none
      integer,intent(in)::nmb
      integer,intent(out)::istat
      type(i_matrix ),pointer::pmatrix
      type(i_element ),pointer ::element_tester,next_element
      istat = -1
      call iaccess_matrix (pmatrix,nmb,istat)
      if (istat.ne.0) return
      element_tester=>pmatrix% i_element_start 
      if(.not.associated(element_tester%pntr)) then
         call  dealloc_i_element (element_tester%number,&
                              pmatrix,istat)
         if (istat.ne.0) return
      else
         next_element=>element_tester%pntr
         do while((associated(next_element)))
           call  dealloc_i_element (element_tester%number,&
                                 pmatrix,istat)
            if (istat.ne.0) return
            element_tester=>next_element
            next_element=>element_tester%pntr
         end do
         call  dealloc_i_element (element_tester%number,&
                              pmatrix,istat)
         if (istat.ne.0) return
      end if
      call  dealloc_i_matrix (nmb,istat)
      if (istat.ne.0) return
      istat = 0
      return
      end subroutine i_dealloc 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine new_s_matrix (nmb,Mb,ierr)
      implicit none
      integer ,intent(out)::nmb,ierr
      integer ,intent(in)::Mb
      type(s_matrix ),pointer::matrix_insert
      if (.NOT. sins_init ) then
         nullify(s_matrix_start )
         sins_init  = .TRUE.
      end if
      if (.not.associated(s_matrix_start )) then
         allocate(s_matrix_start ,STAT=ierr)
         s_matrix_start %number= SSP_MATRIX     
         s_matrix_start %number=- s_matrix_start %number    
         nullify(s_matrix_start %pntr)      
      else
         allocate(matrix_insert,STAT=ierr)
         matrix_insert%number= s_matrix_start %number-no_of_types
         matrix_insert%pntr=> s_matrix_start 
         s_matrix_start => matrix_insert    
      end if
      s_matrix_start %DIM=0
      s_matrix_start %property=blas_general+blas_one_base+blas_col_major
      s_matrix_start %new = 1    !new=0:blas_open_handle, new=1: blas_new_handle
      s_matrix_start %format=''
      nullify(s_matrix_start %sub_rows,s_matrix_start %sub_cols)
      nullify(s_matrix_start % s_element_start )
      allocate(s_matrix_start %trb(Mb),s_matrix_start %tre(Mb)) 
      nmb= s_matrix_start %number
      ierr=0
      end subroutine new_s_matrix       
!*
      subroutine dealloc_s_matrix  (nmb,ierr)
      implicit none
      integer ,intent(in)::nmb
      integer ,intent(out)::ierr
      type(s_matrix ),pointer ::matrix_precedent,matrix_tester
      ierr=-1
      if(.not.associated(s_matrix_start %pntr)) then
         if(s_matrix_start %number.eq.nmb) then
            deallocate(s_matrix_start %tre,s_matrix_start %trb,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            deallocate(s_matrix_start ,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            nullify(s_matrix_start )
            ierr=0
            return
         end if
      else
         matrix_tester=> s_matrix_start 
         if(matrix_tester%number.eq.nmb) then
            s_matrix_start =>matrix_tester%pntr
            deallocate(matrix_tester%tre,matrix_tester%trb,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            deallocate(matrix_tester,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            ierr=0
            return
         endif
         matrix_precedent=> s_matrix_start 
         matrix_tester=> s_matrix_start %pntr
         do while((associated(matrix_tester)))
            if(matrix_tester%number.eq.nmb) then
               matrix_precedent%pntr=>matrix_tester%pntr
               deallocate(matrix_tester%tre,matrix_tester%trb,STAT=ierr)
               if (ierr.ne.0) then
                  ierr = blas_error_memdeloc
                  return
               end if
               deallocate(matrix_tester,STAT=ierr)
               if (ierr.ne.0) then
                  ierr = blas_error_memdeloc
                  return
               end if
               ierr=0
               return
            else
               matrix_precedent=>matrix_tester
               matrix_tester=>matrix_tester%pntr
            end if
         end do
      end if
      end subroutine dealloc_s_matrix 
!*
      subroutine saccess_matrix (pmatrix,nmb,istat)
      implicit none
      type(s_matrix ),pointer ::pmatrix
      integer,intent(out)::istat
      integer,intent(in) ::nmb
      type(s_matrix ),pointer ::matrix_tester
      istat=-1
      matrix_tester=> s_matrix_start        
      do while((matrix_tester%number.ne.nmb).and.&
               (associated(matrix_tester%pntr)))
         matrix_tester => matrix_tester%pntr
      end do
      if (matrix_tester%number.eq.nmb) then
         pmatrix => matrix_tester
         istat = 0 
         return
      else
         nullify(pmatrix)
         istat = blas_error_param
      end if
      end  subroutine saccess_matrix 
!*
      subroutine new_s_element (pmatrix,nmb_element,istat)
      implicit none
      type(s_matrix ),pointer::pmatrix
      integer,intent(out)::nmb_element,istat
      type(s_element ),pointer::element_insert
      integer :: ierr
      istat = -1
      if (.not.associated(pmatrix% s_element_start )) then
         pmatrix%new=0 !status changed to blas_open_handle     
         allocate(pmatrix% s_element_start ,STAT=ierr)
         if (ierr.ne.0) then
            istat = blas_error_memalloc
            return
         end if
         pmatrix% s_element_start %number=1 !will certainly changed
         nullify(pmatrix% s_element_start %pntr) 
      else
         allocate(element_insert,STAT=ierr)
         if (ierr.ne.0) then
            istat = blas_error_memalloc
            return
         end if
         element_insert%pntr=>pmatrix% s_element_start 
         element_insert%number=pmatrix% s_element_start %number+1       
         pmatrix% s_element_start => element_insert 
      end if
      select case(pmatrix%format)
         case('normal')
            pmatrix% s_element_start %contents%pntin%value=0
            pmatrix% s_element_start %contents%pntin%row_ind=-1
            pmatrix% s_element_start %contents%pntin%col_ind=-1 
            nullify(pmatrix% s_element_start %contents%blin%value)
            nullify(pmatrix% s_element_start %contents%vblin%value)
         case('block')
            nullify(pmatrix% s_element_start %contents%blin%value)
            nullify(pmatrix% s_element_start %contents%vblin%value)
            pmatrix% s_element_start %contents%blin%row_block_ind=-1
            pmatrix% s_element_start %contents%blin%col_block_ind=-1
         case('vblock')
            nullify(pmatrix% s_element_start %contents%blin%value)
            nullify(pmatrix% s_element_start %contents%vblin%value)
            pmatrix% s_element_start %contents%vblin%row_vblock_ind=-1
            pmatrix% s_element_start %contents%vblin%col_vblock_ind=-1
         case default 
            istat = blas_error_param
            return
      end select
      nmb_element=pmatrix% s_element_start %number
      istat=0
      end subroutine new_s_element  
!*
      subroutine dealloc_s_element  (nmb_element,pmatrix,istat)
      implicit none
      integer ,intent(in)::nmb_element
      type(s_matrix ),pointer::pmatrix
      integer ,intent(out)::istat
      type(s_element ),pointer ::element_tester
      integer::ierr
      istat=-1
      if(.not.associated( pmatrix% s_element_start %pntr)) then
         if(pmatrix% s_element_start %number.eq.nmb_element) then
         if(associated(pmatrix% s_element_start %contents%vblin%value))&
             then
             deallocate(pmatrix% s_element_start %contents%vblin%value,&
                        STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end  if
          if(associated(pmatrix% s_element_start %contents%blin%value))& 
            then
            deallocate(pmatrix% s_element_start %contents%blin%value,&
                       STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if
            if (associated(pmatrix% s_element_start )) then
               deallocate(pmatrix% s_element_start ,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if 
            nullify(pmatrix% s_element_start )
         end if      
         istat = 0 
         return
      else 
         element_tester=>pmatrix% s_element_start 
         if(element_tester%number.eq.nmb_element) then
            pmatrix% s_element_start =>element_tester%pntr
            if(associated(element_tester%contents%vblin%value)) then
             deallocate(element_tester%contents%vblin%value,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if
            if(associated(element_tester%contents%blin%value)) then
             deallocate(element_tester%contents%blin%value,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end  if
            if (associated(element_tester)) then
               deallocate(element_tester,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if
            istat=0
            return
         endif
         element_tester=>pmatrix% s_element_start %pntr
         do while((associated(element_tester)))
            if(element_tester%number.eq.nmb_element) then
               if(associated(element_tester%contents%vblin%value)) then
                  deallocate(element_tester%contents%vblin%value,&
                             STAT=ierr)
                  if (ierr.ne.0) then
                     istat = blas_error_memdeloc
                     return
                  end if
               end if
               if(associated(element_tester%contents%blin%value)) then
                  deallocate(element_tester%contents%blin%value,&
                             STAT=ierr)
                  if (ierr.ne.0) then
                     istat = blas_error_memdeloc
                     return
                  end if
               end  if
               if (associated(element_tester)) then
                  deallocate(element_tester,STAT=ierr)
                  if (ierr.ne.0) then
                     istat = blas_error_memdeloc
                     return
                  end if
               end if
               istat=0
               return
            else
               element_tester=>element_tester%pntr
            end if
         end do
      end if
      end subroutine dealloc_s_element   
!*
      subroutine saccess_element (pelement,nmb_element,&
                                 pmatrix,istat)
      implicit none
      type(s_inelement ),pointer::pelement
      integer,intent(in) ::nmb_element
      type(s_matrix ),pointer::pmatrix
      integer,intent(out)::istat
      type(s_element ),pointer ::element_tester
      istat=-1
      element_tester=>pmatrix% s_element_start        
      do while((element_tester%number.ne.nmb_element)&
               .and.(associated(element_tester%pntr)))
         element_tester => element_tester%pntr
      end do
      if (element_tester%number.eq.nmb_element) then
         pelement => element_tester%contents
         istat = 0 
         return
      else
         nullify(pelement)
         istat = blas_error_param
         return
      end if  
      end  subroutine  saccess_element 
!*
      subroutine s_element_num  (nmb_element,pmatrix,i,j,istat)
      implicit none
      integer,intent(out),target::nmb_element
      type(s_matrix ),pointer::pmatrix
      integer ,intent(in)::i,j
      integer,intent(out)::istat
      type(s_element ),pointer ::element_tester
      logical:: finder
      istat = -1
      select case(pmatrix%format)
      case('normal')
         element_tester=>pmatrix% s_element_start  
         if(.not.associated( element_tester%pntr)) then 
            if((element_tester%contents%pntin%row_ind.eq.i)&
            .and.(element_tester%contents%pntin%col_ind.eq.j))then
               nmb_element=element_tester%number        
            else
               nmb_element=0 
            end if
         else  
            finder=.false.
            do while((associated(element_tester%pntr))&
                     .and.(.not.finder))
               if((element_tester%contents%pntin%row_ind.eq.i)&
              .and.(element_tester%contents%pntin%col_ind.eq.j))then 
                  finder=.true.
               else
                  element_tester => element_tester%pntr
               end if
            end do
            if((element_tester%contents%pntin%row_ind.eq.i)&
            .and.(element_tester%contents%pntin%col_ind.eq.j)) then   
               nmb_element=element_tester%number      
            else
               nmb_element=0 
            end if
         end if
      case('block')   
         element_tester=>pmatrix% s_element_start  
         if(.not.associated( element_tester%pntr)) then 
            if((element_tester%contents%blin%row_block_ind.eq.i)&
           .and.(element_tester%contents%blin%col_block_ind.eq.j))then 
               nmb_element=element_tester%number        
            else
               nmb_element=0 
            end if
         else  
            finder=.false.
            do while((associated(element_tester%pntr))&
                     .and.(.not.finder))
               if((element_tester%contents%blin%row_block_ind.eq.i)&
           .and.(element_tester%contents%blin%col_block_ind.eq.j))then 
                  finder=.true.
               else
                  element_tester => element_tester%pntr
               end if
            end do
            if((element_tester%contents%blin%row_block_ind.eq.i)&
          .and.(element_tester%contents%blin%col_block_ind.eq.j))then   
               nmb_element=element_tester%number      
            else
               nmb_element=0 
            end if
         end if
      case('vblock')
         element_tester=>pmatrix% s_element_start  
         if(.not.associated( element_tester%pntr)) then 
            if((element_tester%contents%vblin%row_vblock_ind.eq.i)&
         .and.(element_tester%contents%vblin%col_vblock_ind.eq.j))then 
               nmb_element=element_tester%number        
            else
               nmb_element=0 
            end if
         else  
            finder=.false.
            do while((associated(element_tester%pntr)).and.&
                     (.not.finder))
            if((element_tester%contents%vblin%row_vblock_ind.eq.i)&
         .and.(element_tester%contents%vblin%col_vblock_ind.eq.j))then 
                  finder=.true.
               else
                  element_tester => element_tester%pntr
               end if
            end do
            if((element_tester%contents%vblin%row_vblock_ind.eq.i)&
         .and.(element_tester%contents%vblin%col_vblock_ind.eq.j))then   
               nmb_element=element_tester%number      
            else
               nmb_element=0 
            end if
         end if
      case default 
         istat = blas_error_param
         return
      end select
      istat = 0
      end subroutine s_element_num 
!*
      subroutine s_dealloc (nmb,istat)
      implicit none
      integer,intent(in)::nmb
      integer,intent(out)::istat
      type(s_matrix ),pointer::pmatrix
      type(s_element ),pointer ::element_tester,next_element
      istat = -1
      call saccess_matrix (pmatrix,nmb,istat)
      if (istat.ne.0) return
      element_tester=>pmatrix% s_element_start 
      if(.not.associated(element_tester%pntr)) then
         call  dealloc_s_element (element_tester%number,&
                              pmatrix,istat)
         if (istat.ne.0) return
      else
         next_element=>element_tester%pntr
         do while((associated(next_element)))
           call  dealloc_s_element (element_tester%number,&
                                 pmatrix,istat)
            if (istat.ne.0) return
            element_tester=>next_element
            next_element=>element_tester%pntr
         end do
         call  dealloc_s_element (element_tester%number,&
                              pmatrix,istat)
         if (istat.ne.0) return
      end if
      call  dealloc_s_matrix (nmb,istat)
      if (istat.ne.0) return
      istat = 0
      return
      end subroutine s_dealloc 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine new_d_matrix (nmb,Mb,ierr)
      implicit none
      integer ,intent(out)::nmb,ierr
      integer ,intent(in)::Mb
      type(d_matrix ),pointer::matrix_insert
      if (.NOT. dins_init ) then
         nullify(d_matrix_start )
         dins_init  = .TRUE.
      end if
      if (.not.associated(d_matrix_start )) then
         allocate(d_matrix_start ,STAT=ierr)
         d_matrix_start %number= DSP_MATRIX     
         d_matrix_start %number=- d_matrix_start %number    
         nullify(d_matrix_start %pntr)      
      else
         allocate(matrix_insert,STAT=ierr)
         matrix_insert%number= d_matrix_start %number-no_of_types
         matrix_insert%pntr=> d_matrix_start 
         d_matrix_start => matrix_insert    
      end if
      d_matrix_start %DIM=0
      d_matrix_start %property=blas_general+blas_one_base+blas_col_major
      d_matrix_start %new = 1    !new=0:blas_open_handle, new=1: blas_new_handle
      d_matrix_start %format=''
      nullify(d_matrix_start %sub_rows,d_matrix_start %sub_cols)
      nullify(d_matrix_start % d_element_start )
      allocate(d_matrix_start %trb(Mb),d_matrix_start %tre(Mb)) 
      nmb= d_matrix_start %number
      ierr=0
      end subroutine new_d_matrix       
!*
      subroutine dealloc_d_matrix  (nmb,ierr)
      implicit none
      integer ,intent(in)::nmb
      integer ,intent(out)::ierr
      type(d_matrix ),pointer ::matrix_precedent,matrix_tester
      ierr=-1
      if(.not.associated(d_matrix_start %pntr)) then
         if(d_matrix_start %number.eq.nmb) then
            deallocate(d_matrix_start %tre,d_matrix_start %trb,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            deallocate(d_matrix_start ,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            nullify(d_matrix_start )
            ierr=0
            return
         end if
      else
         matrix_tester=> d_matrix_start 
         if(matrix_tester%number.eq.nmb) then
            d_matrix_start =>matrix_tester%pntr
            deallocate(matrix_tester%tre,matrix_tester%trb,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            deallocate(matrix_tester,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            ierr=0
            return
         endif
         matrix_precedent=> d_matrix_start 
         matrix_tester=> d_matrix_start %pntr
         do while((associated(matrix_tester)))
            if(matrix_tester%number.eq.nmb) then
               matrix_precedent%pntr=>matrix_tester%pntr
               deallocate(matrix_tester%tre,matrix_tester%trb,STAT=ierr)
               if (ierr.ne.0) then
                  ierr = blas_error_memdeloc
                  return
               end if
               deallocate(matrix_tester,STAT=ierr)
               if (ierr.ne.0) then
                  ierr = blas_error_memdeloc
                  return
               end if
               ierr=0
               return
            else
               matrix_precedent=>matrix_tester
               matrix_tester=>matrix_tester%pntr
            end if
         end do
      end if
      end subroutine dealloc_d_matrix 
!*
      subroutine daccess_matrix (pmatrix,nmb,istat)
      implicit none
      type(d_matrix ),pointer ::pmatrix
      integer,intent(out)::istat
      integer,intent(in) ::nmb
      type(d_matrix ),pointer ::matrix_tester
      istat=-1
      matrix_tester=> d_matrix_start        
      do while((matrix_tester%number.ne.nmb).and.&
               (associated(matrix_tester%pntr)))
         matrix_tester => matrix_tester%pntr
      end do
      if (matrix_tester%number.eq.nmb) then
         pmatrix => matrix_tester
         istat = 0 
         return
      else
         nullify(pmatrix)
         istat = blas_error_param
      end if
      end  subroutine daccess_matrix 
!*
      subroutine new_d_element (pmatrix,nmb_element,istat)
      implicit none
      type(d_matrix ),pointer::pmatrix
      integer,intent(out)::nmb_element,istat
      type(d_element ),pointer::element_insert
      integer :: ierr
      istat = -1
      if (.not.associated(pmatrix% d_element_start )) then
         pmatrix%new=0 !status changed to blas_open_handle     
         allocate(pmatrix% d_element_start ,STAT=ierr)
         if (ierr.ne.0) then
            istat = blas_error_memalloc
            return
         end if
         pmatrix% d_element_start %number=1 !will certainly changed
         nullify(pmatrix% d_element_start %pntr) 
      else
         allocate(element_insert,STAT=ierr)
         if (ierr.ne.0) then
            istat = blas_error_memalloc
            return
         end if
         element_insert%pntr=>pmatrix% d_element_start 
         element_insert%number=pmatrix% d_element_start %number+1       
         pmatrix% d_element_start => element_insert 
      end if
      select case(pmatrix%format)
         case('normal')
            pmatrix% d_element_start %contents%pntin%value=0
            pmatrix% d_element_start %contents%pntin%row_ind=-1
            pmatrix% d_element_start %contents%pntin%col_ind=-1 
            nullify(pmatrix% d_element_start %contents%blin%value)
            nullify(pmatrix% d_element_start %contents%vblin%value)
         case('block')
            nullify(pmatrix% d_element_start %contents%blin%value)
            nullify(pmatrix% d_element_start %contents%vblin%value)
            pmatrix% d_element_start %contents%blin%row_block_ind=-1
            pmatrix% d_element_start %contents%blin%col_block_ind=-1
         case('vblock')
            nullify(pmatrix% d_element_start %contents%blin%value)
            nullify(pmatrix% d_element_start %contents%vblin%value)
            pmatrix% d_element_start %contents%vblin%row_vblock_ind=-1
            pmatrix% d_element_start %contents%vblin%col_vblock_ind=-1
         case default 
            istat = blas_error_param
            return
      end select
      nmb_element=pmatrix% d_element_start %number
      istat=0
      end subroutine new_d_element  
!*
      subroutine dealloc_d_element  (nmb_element,pmatrix,istat)
      implicit none
      integer ,intent(in)::nmb_element
      type(d_matrix ),pointer::pmatrix
      integer ,intent(out)::istat
      type(d_element ),pointer ::element_tester
      integer::ierr
      istat=-1
      if(.not.associated( pmatrix% d_element_start %pntr)) then
         if(pmatrix% d_element_start %number.eq.nmb_element) then
         if(associated(pmatrix% d_element_start %contents%vblin%value))&
             then
             deallocate(pmatrix% d_element_start %contents%vblin%value,&
                        STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end  if
          if(associated(pmatrix% d_element_start %contents%blin%value))& 
            then
            deallocate(pmatrix% d_element_start %contents%blin%value,&
                       STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if
            if (associated(pmatrix% d_element_start )) then
               deallocate(pmatrix% d_element_start ,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if 
            nullify(pmatrix% d_element_start )
         end if      
         istat = 0 
         return
      else 
         element_tester=>pmatrix% d_element_start 
         if(element_tester%number.eq.nmb_element) then
            pmatrix% d_element_start =>element_tester%pntr
            if(associated(element_tester%contents%vblin%value)) then
             deallocate(element_tester%contents%vblin%value,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if
            if(associated(element_tester%contents%blin%value)) then
             deallocate(element_tester%contents%blin%value,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end  if
            if (associated(element_tester)) then
               deallocate(element_tester,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if
            istat=0
            return
         endif
         element_tester=>pmatrix% d_element_start %pntr
         do while((associated(element_tester)))
            if(element_tester%number.eq.nmb_element) then
               if(associated(element_tester%contents%vblin%value)) then
                  deallocate(element_tester%contents%vblin%value,&
                             STAT=ierr)
                  if (ierr.ne.0) then
                     istat = blas_error_memdeloc
                     return
                  end if
               end if
               if(associated(element_tester%contents%blin%value)) then
                  deallocate(element_tester%contents%blin%value,&
                             STAT=ierr)
                  if (ierr.ne.0) then
                     istat = blas_error_memdeloc
                     return
                  end if
               end  if
               if (associated(element_tester)) then
                  deallocate(element_tester,STAT=ierr)
                  if (ierr.ne.0) then
                     istat = blas_error_memdeloc
                     return
                  end if
               end if
               istat=0
               return
            else
               element_tester=>element_tester%pntr
            end if
         end do
      end if
      end subroutine dealloc_d_element   
!*
      subroutine daccess_element (pelement,nmb_element,&
                                 pmatrix,istat)
      implicit none
      type(d_inelement ),pointer::pelement
      integer,intent(in) ::nmb_element
      type(d_matrix ),pointer::pmatrix
      integer,intent(out)::istat
      type(d_element ),pointer ::element_tester
      istat=-1
      element_tester=>pmatrix% d_element_start        
      do while((element_tester%number.ne.nmb_element)&
               .and.(associated(element_tester%pntr)))
         element_tester => element_tester%pntr
      end do
      if (element_tester%number.eq.nmb_element) then
         pelement => element_tester%contents
         istat = 0 
         return
      else
         nullify(pelement)
         istat = blas_error_param
         return
      end if  
      end  subroutine  daccess_element 
!*
      subroutine d_element_num  (nmb_element,pmatrix,i,j,istat)
      implicit none
      integer,intent(out),target::nmb_element
      type(d_matrix ),pointer::pmatrix
      integer ,intent(in)::i,j
      integer,intent(out)::istat
      type(d_element ),pointer ::element_tester
      logical:: finder
      istat = -1
      select case(pmatrix%format)
      case('normal')
         element_tester=>pmatrix% d_element_start  
         if(.not.associated( element_tester%pntr)) then 
            if((element_tester%contents%pntin%row_ind.eq.i)&
            .and.(element_tester%contents%pntin%col_ind.eq.j))then
               nmb_element=element_tester%number        
            else
               nmb_element=0 
            end if
         else  
            finder=.false.
            do while((associated(element_tester%pntr))&
                     .and.(.not.finder))
               if((element_tester%contents%pntin%row_ind.eq.i)&
              .and.(element_tester%contents%pntin%col_ind.eq.j))then 
                  finder=.true.
               else
                  element_tester => element_tester%pntr
               end if
            end do
            if((element_tester%contents%pntin%row_ind.eq.i)&
            .and.(element_tester%contents%pntin%col_ind.eq.j)) then   
               nmb_element=element_tester%number      
            else
               nmb_element=0 
            end if
         end if
      case('block')   
         element_tester=>pmatrix% d_element_start  
         if(.not.associated( element_tester%pntr)) then 
            if((element_tester%contents%blin%row_block_ind.eq.i)&
           .and.(element_tester%contents%blin%col_block_ind.eq.j))then 
               nmb_element=element_tester%number        
            else
               nmb_element=0 
            end if
         else  
            finder=.false.
            do while((associated(element_tester%pntr))&
                     .and.(.not.finder))
               if((element_tester%contents%blin%row_block_ind.eq.i)&
           .and.(element_tester%contents%blin%col_block_ind.eq.j))then 
                  finder=.true.
               else
                  element_tester => element_tester%pntr
               end if
            end do
            if((element_tester%contents%blin%row_block_ind.eq.i)&
          .and.(element_tester%contents%blin%col_block_ind.eq.j))then   
               nmb_element=element_tester%number      
            else
               nmb_element=0 
            end if
         end if
      case('vblock')
         element_tester=>pmatrix% d_element_start  
         if(.not.associated( element_tester%pntr)) then 
            if((element_tester%contents%vblin%row_vblock_ind.eq.i)&
         .and.(element_tester%contents%vblin%col_vblock_ind.eq.j))then 
               nmb_element=element_tester%number        
            else
               nmb_element=0 
            end if
         else  
            finder=.false.
            do while((associated(element_tester%pntr)).and.&
                     (.not.finder))
            if((element_tester%contents%vblin%row_vblock_ind.eq.i)&
         .and.(element_tester%contents%vblin%col_vblock_ind.eq.j))then 
                  finder=.true.
               else
                  element_tester => element_tester%pntr
               end if
            end do
            if((element_tester%contents%vblin%row_vblock_ind.eq.i)&
         .and.(element_tester%contents%vblin%col_vblock_ind.eq.j))then   
               nmb_element=element_tester%number      
            else
               nmb_element=0 
            end if
         end if
      case default 
         istat = blas_error_param
         return
      end select
      istat = 0
      end subroutine d_element_num 
!*
      subroutine d_dealloc (nmb,istat)
      implicit none
      integer,intent(in)::nmb
      integer,intent(out)::istat
      type(d_matrix ),pointer::pmatrix
      type(d_element ),pointer ::element_tester,next_element
      istat = -1
      call daccess_matrix (pmatrix,nmb,istat)
      if (istat.ne.0) return
      element_tester=>pmatrix% d_element_start 
      if(.not.associated(element_tester%pntr)) then
         call  dealloc_d_element (element_tester%number,&
                              pmatrix,istat)
         if (istat.ne.0) return
      else
         next_element=>element_tester%pntr
         do while((associated(next_element)))
           call  dealloc_d_element (element_tester%number,&
                                 pmatrix,istat)
            if (istat.ne.0) return
            element_tester=>next_element
            next_element=>element_tester%pntr
         end do
         call  dealloc_d_element (element_tester%number,&
                              pmatrix,istat)
         if (istat.ne.0) return
      end if
      call  dealloc_d_matrix (nmb,istat)
      if (istat.ne.0) return
      istat = 0
      return
      end subroutine d_dealloc 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine new_c_matrix (nmb,Mb,ierr)
      implicit none
      integer ,intent(out)::nmb,ierr
      integer ,intent(in)::Mb
      type(c_matrix ),pointer::matrix_insert
      if (.NOT. cins_init ) then
         nullify(c_matrix_start )
         cins_init  = .TRUE.
      end if
      if (.not.associated(c_matrix_start )) then
         allocate(c_matrix_start ,STAT=ierr)
         c_matrix_start %number= CSP_MATRIX     
         c_matrix_start %number=- c_matrix_start %number    
         nullify(c_matrix_start %pntr)      
      else
         allocate(matrix_insert,STAT=ierr)
         matrix_insert%number= c_matrix_start %number-no_of_types
         matrix_insert%pntr=> c_matrix_start 
         c_matrix_start => matrix_insert    
      end if
      c_matrix_start %DIM=0
      c_matrix_start %property=blas_general+blas_one_base+blas_col_major
      c_matrix_start %new = 1    !new=0:blas_open_handle, new=1: blas_new_handle
      c_matrix_start %format=''
      nullify(c_matrix_start %sub_rows,c_matrix_start %sub_cols)
      nullify(c_matrix_start % c_element_start )
      allocate(c_matrix_start %trb(Mb),c_matrix_start %tre(Mb)) 
      nmb= c_matrix_start %number
      ierr=0
      end subroutine new_c_matrix       
!*
      subroutine dealloc_c_matrix  (nmb,ierr)
      implicit none
      integer ,intent(in)::nmb
      integer ,intent(out)::ierr
      type(c_matrix ),pointer ::matrix_precedent,matrix_tester
      ierr=-1
      if(.not.associated(c_matrix_start %pntr)) then
         if(c_matrix_start %number.eq.nmb) then
            deallocate(c_matrix_start %tre,c_matrix_start %trb,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            deallocate(c_matrix_start ,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            nullify(c_matrix_start )
            ierr=0
            return
         end if
      else
         matrix_tester=> c_matrix_start 
         if(matrix_tester%number.eq.nmb) then
            c_matrix_start =>matrix_tester%pntr
            deallocate(matrix_tester%tre,matrix_tester%trb,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            deallocate(matrix_tester,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            ierr=0
            return
         endif
         matrix_precedent=> c_matrix_start 
         matrix_tester=> c_matrix_start %pntr
         do while((associated(matrix_tester)))
            if(matrix_tester%number.eq.nmb) then
               matrix_precedent%pntr=>matrix_tester%pntr
               deallocate(matrix_tester%tre,matrix_tester%trb,STAT=ierr)
               if (ierr.ne.0) then
                  ierr = blas_error_memdeloc
                  return
               end if
               deallocate(matrix_tester,STAT=ierr)
               if (ierr.ne.0) then
                  ierr = blas_error_memdeloc
                  return
               end if
               ierr=0
               return
            else
               matrix_precedent=>matrix_tester
               matrix_tester=>matrix_tester%pntr
            end if
         end do
      end if
      end subroutine dealloc_c_matrix 
!*
      subroutine caccess_matrix (pmatrix,nmb,istat)
      implicit none
      type(c_matrix ),pointer ::pmatrix
      integer,intent(out)::istat
      integer,intent(in) ::nmb
      type(c_matrix ),pointer ::matrix_tester
      istat=-1
      matrix_tester=> c_matrix_start        
      do while((matrix_tester%number.ne.nmb).and.&
               (associated(matrix_tester%pntr)))
         matrix_tester => matrix_tester%pntr
      end do
      if (matrix_tester%number.eq.nmb) then
         pmatrix => matrix_tester
         istat = 0 
         return
      else
         nullify(pmatrix)
         istat = blas_error_param
      end if
      end  subroutine caccess_matrix 
!*
      subroutine new_c_element (pmatrix,nmb_element,istat)
      implicit none
      type(c_matrix ),pointer::pmatrix
      integer,intent(out)::nmb_element,istat
      type(c_element ),pointer::element_insert
      integer :: ierr
      istat = -1
      if (.not.associated(pmatrix% c_element_start )) then
         pmatrix%new=0 !status changed to blas_open_handle     
         allocate(pmatrix% c_element_start ,STAT=ierr)
         if (ierr.ne.0) then
            istat = blas_error_memalloc
            return
         end if
         pmatrix% c_element_start %number=1 !will certainly changed
         nullify(pmatrix% c_element_start %pntr) 
      else
         allocate(element_insert,STAT=ierr)
         if (ierr.ne.0) then
            istat = blas_error_memalloc
            return
         end if
         element_insert%pntr=>pmatrix% c_element_start 
         element_insert%number=pmatrix% c_element_start %number+1       
         pmatrix% c_element_start => element_insert 
      end if
      select case(pmatrix%format)
         case('normal')
            pmatrix% c_element_start %contents%pntin%value=0
            pmatrix% c_element_start %contents%pntin%row_ind=-1
            pmatrix% c_element_start %contents%pntin%col_ind=-1 
            nullify(pmatrix% c_element_start %contents%blin%value)
            nullify(pmatrix% c_element_start %contents%vblin%value)
         case('block')
            nullify(pmatrix% c_element_start %contents%blin%value)
            nullify(pmatrix% c_element_start %contents%vblin%value)
            pmatrix% c_element_start %contents%blin%row_block_ind=-1
            pmatrix% c_element_start %contents%blin%col_block_ind=-1
         case('vblock')
            nullify(pmatrix% c_element_start %contents%blin%value)
            nullify(pmatrix% c_element_start %contents%vblin%value)
            pmatrix% c_element_start %contents%vblin%row_vblock_ind=-1
            pmatrix% c_element_start %contents%vblin%col_vblock_ind=-1
         case default 
            istat = blas_error_param
            return
      end select
      nmb_element=pmatrix% c_element_start %number
      istat=0
      end subroutine new_c_element  
!*
      subroutine dealloc_c_element  (nmb_element,pmatrix,istat)
      implicit none
      integer ,intent(in)::nmb_element
      type(c_matrix ),pointer::pmatrix
      integer ,intent(out)::istat
      type(c_element ),pointer ::element_tester
      integer::ierr
      istat=-1
      if(.not.associated( pmatrix% c_element_start %pntr)) then
         if(pmatrix% c_element_start %number.eq.nmb_element) then
         if(associated(pmatrix% c_element_start %contents%vblin%value))&
             then
             deallocate(pmatrix% c_element_start %contents%vblin%value,&
                        STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end  if
          if(associated(pmatrix% c_element_start %contents%blin%value))& 
            then
            deallocate(pmatrix% c_element_start %contents%blin%value,&
                       STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if
            if (associated(pmatrix% c_element_start )) then
               deallocate(pmatrix% c_element_start ,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if 
            nullify(pmatrix% c_element_start )
         end if      
         istat = 0 
         return
      else 
         element_tester=>pmatrix% c_element_start 
         if(element_tester%number.eq.nmb_element) then
            pmatrix% c_element_start =>element_tester%pntr
            if(associated(element_tester%contents%vblin%value)) then
             deallocate(element_tester%contents%vblin%value,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if
            if(associated(element_tester%contents%blin%value)) then
             deallocate(element_tester%contents%blin%value,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end  if
            if (associated(element_tester)) then
               deallocate(element_tester,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if
            istat=0
            return
         endif
         element_tester=>pmatrix% c_element_start %pntr
         do while((associated(element_tester)))
            if(element_tester%number.eq.nmb_element) then
               if(associated(element_tester%contents%vblin%value)) then
                  deallocate(element_tester%contents%vblin%value,&
                             STAT=ierr)
                  if (ierr.ne.0) then
                     istat = blas_error_memdeloc
                     return
                  end if
               end if
               if(associated(element_tester%contents%blin%value)) then
                  deallocate(element_tester%contents%blin%value,&
                             STAT=ierr)
                  if (ierr.ne.0) then
                     istat = blas_error_memdeloc
                     return
                  end if
               end  if
               if (associated(element_tester)) then
                  deallocate(element_tester,STAT=ierr)
                  if (ierr.ne.0) then
                     istat = blas_error_memdeloc
                     return
                  end if
               end if
               istat=0
               return
            else
               element_tester=>element_tester%pntr
            end if
         end do
      end if
      end subroutine dealloc_c_element   
!*
      subroutine caccess_element (pelement,nmb_element,&
                                 pmatrix,istat)
      implicit none
      type(c_inelement ),pointer::pelement
      integer,intent(in) ::nmb_element
      type(c_matrix ),pointer::pmatrix
      integer,intent(out)::istat
      type(c_element ),pointer ::element_tester
      istat=-1
      element_tester=>pmatrix% c_element_start        
      do while((element_tester%number.ne.nmb_element)&
               .and.(associated(element_tester%pntr)))
         element_tester => element_tester%pntr
      end do
      if (element_tester%number.eq.nmb_element) then
         pelement => element_tester%contents
         istat = 0 
         return
      else
         nullify(pelement)
         istat = blas_error_param
         return
      end if  
      end  subroutine  caccess_element 
!*
      subroutine c_element_num  (nmb_element,pmatrix,i,j,istat)
      implicit none
      integer,intent(out),target::nmb_element
      type(c_matrix ),pointer::pmatrix
      integer ,intent(in)::i,j
      integer,intent(out)::istat
      type(c_element ),pointer ::element_tester
      logical:: finder
      istat = -1
      select case(pmatrix%format)
      case('normal')
         element_tester=>pmatrix% c_element_start  
         if(.not.associated( element_tester%pntr)) then 
            if((element_tester%contents%pntin%row_ind.eq.i)&
            .and.(element_tester%contents%pntin%col_ind.eq.j))then
               nmb_element=element_tester%number        
            else
               nmb_element=0 
            end if
         else  
            finder=.false.
            do while((associated(element_tester%pntr))&
                     .and.(.not.finder))
               if((element_tester%contents%pntin%row_ind.eq.i)&
              .and.(element_tester%contents%pntin%col_ind.eq.j))then 
                  finder=.true.
               else
                  element_tester => element_tester%pntr
               end if
            end do
            if((element_tester%contents%pntin%row_ind.eq.i)&
            .and.(element_tester%contents%pntin%col_ind.eq.j)) then   
               nmb_element=element_tester%number      
            else
               nmb_element=0 
            end if
         end if
      case('block')   
         element_tester=>pmatrix% c_element_start  
         if(.not.associated( element_tester%pntr)) then 
            if((element_tester%contents%blin%row_block_ind.eq.i)&
           .and.(element_tester%contents%blin%col_block_ind.eq.j))then 
               nmb_element=element_tester%number        
            else
               nmb_element=0 
            end if
         else  
            finder=.false.
            do while((associated(element_tester%pntr))&
                     .and.(.not.finder))
               if((element_tester%contents%blin%row_block_ind.eq.i)&
           .and.(element_tester%contents%blin%col_block_ind.eq.j))then 
                  finder=.true.
               else
                  element_tester => element_tester%pntr
               end if
            end do
            if((element_tester%contents%blin%row_block_ind.eq.i)&
          .and.(element_tester%contents%blin%col_block_ind.eq.j))then   
               nmb_element=element_tester%number      
            else
               nmb_element=0 
            end if
         end if
      case('vblock')
         element_tester=>pmatrix% c_element_start  
         if(.not.associated( element_tester%pntr)) then 
            if((element_tester%contents%vblin%row_vblock_ind.eq.i)&
         .and.(element_tester%contents%vblin%col_vblock_ind.eq.j))then 
               nmb_element=element_tester%number        
            else
               nmb_element=0 
            end if
         else  
            finder=.false.
            do while((associated(element_tester%pntr)).and.&
                     (.not.finder))
            if((element_tester%contents%vblin%row_vblock_ind.eq.i)&
         .and.(element_tester%contents%vblin%col_vblock_ind.eq.j))then 
                  finder=.true.
               else
                  element_tester => element_tester%pntr
               end if
            end do
            if((element_tester%contents%vblin%row_vblock_ind.eq.i)&
         .and.(element_tester%contents%vblin%col_vblock_ind.eq.j))then   
               nmb_element=element_tester%number      
            else
               nmb_element=0 
            end if
         end if
      case default 
         istat = blas_error_param
         return
      end select
      istat = 0
      end subroutine c_element_num 
!*
      subroutine c_dealloc (nmb,istat)
      implicit none
      integer,intent(in)::nmb
      integer,intent(out)::istat
      type(c_matrix ),pointer::pmatrix
      type(c_element ),pointer ::element_tester,next_element
      istat = -1
      call caccess_matrix (pmatrix,nmb,istat)
      if (istat.ne.0) return
      element_tester=>pmatrix% c_element_start 
      if(.not.associated(element_tester%pntr)) then
         call  dealloc_c_element (element_tester%number,&
                              pmatrix,istat)
         if (istat.ne.0) return
      else
         next_element=>element_tester%pntr
         do while((associated(next_element)))
           call  dealloc_c_element (element_tester%number,&
                                 pmatrix,istat)
            if (istat.ne.0) return
            element_tester=>next_element
            next_element=>element_tester%pntr
         end do
         call  dealloc_c_element (element_tester%number,&
                              pmatrix,istat)
         if (istat.ne.0) return
      end if
      call  dealloc_c_matrix (nmb,istat)
      if (istat.ne.0) return
      istat = 0
      return
      end subroutine c_dealloc 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine new_z_matrix (nmb,Mb,ierr)
      implicit none
      integer ,intent(out)::nmb,ierr
      integer ,intent(in)::Mb
      type(z_matrix ),pointer::matrix_insert
      if (.NOT. zins_init ) then
         nullify(z_matrix_start )
         zins_init  = .TRUE.
      end if
      if (.not.associated(z_matrix_start )) then
         allocate(z_matrix_start ,STAT=ierr)
         z_matrix_start %number= ZSP_MATRIX     
         z_matrix_start %number=- z_matrix_start %number    
         nullify(z_matrix_start %pntr)      
      else
         allocate(matrix_insert,STAT=ierr)
         matrix_insert%number= z_matrix_start %number-no_of_types
         matrix_insert%pntr=> z_matrix_start 
         z_matrix_start => matrix_insert    
      end if
      z_matrix_start %DIM=0
      z_matrix_start %property=blas_general+blas_one_base+blas_col_major
      z_matrix_start %new = 1    !new=0:blas_open_handle, new=1: blas_new_handle
      z_matrix_start %format=''
      nullify(z_matrix_start %sub_rows,z_matrix_start %sub_cols)
      nullify(z_matrix_start % z_element_start )
      allocate(z_matrix_start %trb(Mb),z_matrix_start %tre(Mb)) 
      nmb= z_matrix_start %number
      ierr=0
      end subroutine new_z_matrix       
!*
      subroutine dealloc_z_matrix  (nmb,ierr)
      implicit none
      integer ,intent(in)::nmb
      integer ,intent(out)::ierr
      type(z_matrix ),pointer ::matrix_precedent,matrix_tester
      ierr=-1
      if(.not.associated(z_matrix_start %pntr)) then
         if(z_matrix_start %number.eq.nmb) then
            deallocate(z_matrix_start %tre,z_matrix_start %trb,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            deallocate(z_matrix_start ,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            nullify(z_matrix_start )
            ierr=0
            return
         end if
      else
         matrix_tester=> z_matrix_start 
         if(matrix_tester%number.eq.nmb) then
            z_matrix_start =>matrix_tester%pntr
            deallocate(matrix_tester%tre,matrix_tester%trb,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            deallocate(matrix_tester,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            ierr=0
            return
         endif
         matrix_precedent=> z_matrix_start 
         matrix_tester=> z_matrix_start %pntr
         do while((associated(matrix_tester)))
            if(matrix_tester%number.eq.nmb) then
               matrix_precedent%pntr=>matrix_tester%pntr
               deallocate(matrix_tester%tre,matrix_tester%trb,STAT=ierr)
               if (ierr.ne.0) then
                  ierr = blas_error_memdeloc
                  return
               end if
               deallocate(matrix_tester,STAT=ierr)
               if (ierr.ne.0) then
                  ierr = blas_error_memdeloc
                  return
               end if
               ierr=0
               return
            else
               matrix_precedent=>matrix_tester
               matrix_tester=>matrix_tester%pntr
            end if
         end do
      end if
      end subroutine dealloc_z_matrix 
!*
      subroutine zaccess_matrix (pmatrix,nmb,istat)
      implicit none
      type(z_matrix ),pointer ::pmatrix
      integer,intent(out)::istat
      integer,intent(in) ::nmb
      type(z_matrix ),pointer ::matrix_tester
      istat=-1
      matrix_tester=> z_matrix_start        
      do while((matrix_tester%number.ne.nmb).and.&
               (associated(matrix_tester%pntr)))
         matrix_tester => matrix_tester%pntr
      end do
      if (matrix_tester%number.eq.nmb) then
         pmatrix => matrix_tester
         istat = 0 
         return
      else
         nullify(pmatrix)
         istat = blas_error_param
      end if
      end  subroutine zaccess_matrix 
!*
      subroutine new_z_element (pmatrix,nmb_element,istat)
      implicit none
      type(z_matrix ),pointer::pmatrix
      integer,intent(out)::nmb_element,istat
      type(z_element ),pointer::element_insert
      integer :: ierr
      istat = -1
      if (.not.associated(pmatrix% z_element_start )) then
         pmatrix%new=0 !status changed to blas_open_handle     
         allocate(pmatrix% z_element_start ,STAT=ierr)
         if (ierr.ne.0) then
            istat = blas_error_memalloc
            return
         end if
         pmatrix% z_element_start %number=1 !will certainly changed
         nullify(pmatrix% z_element_start %pntr) 
      else
         allocate(element_insert,STAT=ierr)
         if (ierr.ne.0) then
            istat = blas_error_memalloc
            return
         end if
         element_insert%pntr=>pmatrix% z_element_start 
         element_insert%number=pmatrix% z_element_start %number+1       
         pmatrix% z_element_start => element_insert 
      end if
      select case(pmatrix%format)
         case('normal')
            pmatrix% z_element_start %contents%pntin%value=0
            pmatrix% z_element_start %contents%pntin%row_ind=-1
            pmatrix% z_element_start %contents%pntin%col_ind=-1 
            nullify(pmatrix% z_element_start %contents%blin%value)
            nullify(pmatrix% z_element_start %contents%vblin%value)
         case('block')
            nullify(pmatrix% z_element_start %contents%blin%value)
            nullify(pmatrix% z_element_start %contents%vblin%value)
            pmatrix% z_element_start %contents%blin%row_block_ind=-1
            pmatrix% z_element_start %contents%blin%col_block_ind=-1
         case('vblock')
            nullify(pmatrix% z_element_start %contents%blin%value)
            nullify(pmatrix% z_element_start %contents%vblin%value)
            pmatrix% z_element_start %contents%vblin%row_vblock_ind=-1
            pmatrix% z_element_start %contents%vblin%col_vblock_ind=-1
         case default 
            istat = blas_error_param
            return
      end select
      nmb_element=pmatrix% z_element_start %number
      istat=0
      end subroutine new_z_element  
!*
      subroutine dealloc_z_element  (nmb_element,pmatrix,istat)
      implicit none
      integer ,intent(in)::nmb_element
      type(z_matrix ),pointer::pmatrix
      integer ,intent(out)::istat
      type(z_element ),pointer ::element_tester
      integer::ierr
      istat=-1
      if(.not.associated( pmatrix% z_element_start %pntr)) then
         if(pmatrix% z_element_start %number.eq.nmb_element) then
         if(associated(pmatrix% z_element_start %contents%vblin%value))&
             then
             deallocate(pmatrix% z_element_start %contents%vblin%value,&
                        STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end  if
          if(associated(pmatrix% z_element_start %contents%blin%value))& 
            then
            deallocate(pmatrix% z_element_start %contents%blin%value,&
                       STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if
            if (associated(pmatrix% z_element_start )) then
               deallocate(pmatrix% z_element_start ,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if 
            nullify(pmatrix% z_element_start )
         end if      
         istat = 0 
         return
      else 
         element_tester=>pmatrix% z_element_start 
         if(element_tester%number.eq.nmb_element) then
            pmatrix% z_element_start =>element_tester%pntr
            if(associated(element_tester%contents%vblin%value)) then
             deallocate(element_tester%contents%vblin%value,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if
            if(associated(element_tester%contents%blin%value)) then
             deallocate(element_tester%contents%blin%value,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end  if
            if (associated(element_tester)) then
               deallocate(element_tester,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if
            istat=0
            return
         endif
         element_tester=>pmatrix% z_element_start %pntr
         do while((associated(element_tester)))
            if(element_tester%number.eq.nmb_element) then
               if(associated(element_tester%contents%vblin%value)) then
                  deallocate(element_tester%contents%vblin%value,&
                             STAT=ierr)
                  if (ierr.ne.0) then
                     istat = blas_error_memdeloc
                     return
                  end if
               end if
               if(associated(element_tester%contents%blin%value)) then
                  deallocate(element_tester%contents%blin%value,&
                             STAT=ierr)
                  if (ierr.ne.0) then
                     istat = blas_error_memdeloc
                     return
                  end if
               end  if
               if (associated(element_tester)) then
                  deallocate(element_tester,STAT=ierr)
                  if (ierr.ne.0) then
                     istat = blas_error_memdeloc
                     return
                  end if
               end if
               istat=0
               return
            else
               element_tester=>element_tester%pntr
            end if
         end do
      end if
      end subroutine dealloc_z_element   
!*
      subroutine zaccess_element (pelement,nmb_element,&
                                 pmatrix,istat)
      implicit none
      type(z_inelement ),pointer::pelement
      integer,intent(in) ::nmb_element
      type(z_matrix ),pointer::pmatrix
      integer,intent(out)::istat
      type(z_element ),pointer ::element_tester
      istat=-1
      element_tester=>pmatrix% z_element_start        
      do while((element_tester%number.ne.nmb_element)&
               .and.(associated(element_tester%pntr)))
         element_tester => element_tester%pntr
      end do
      if (element_tester%number.eq.nmb_element) then
         pelement => element_tester%contents
         istat = 0 
         return
      else
         nullify(pelement)
         istat = blas_error_param
         return
      end if  
      end  subroutine  zaccess_element 
!*
      subroutine z_element_num  (nmb_element,pmatrix,i,j,istat)
      implicit none
      integer,intent(out),target::nmb_element
      type(z_matrix ),pointer::pmatrix
      integer ,intent(in)::i,j
      integer,intent(out)::istat
      type(z_element ),pointer ::element_tester
      logical:: finder
      istat = -1
      select case(pmatrix%format)
      case('normal')
         element_tester=>pmatrix% z_element_start  
         if(.not.associated( element_tester%pntr)) then 
            if((element_tester%contents%pntin%row_ind.eq.i)&
            .and.(element_tester%contents%pntin%col_ind.eq.j))then
               nmb_element=element_tester%number        
            else
               nmb_element=0 
            end if
         else  
            finder=.false.
            do while((associated(element_tester%pntr))&
                     .and.(.not.finder))
               if((element_tester%contents%pntin%row_ind.eq.i)&
              .and.(element_tester%contents%pntin%col_ind.eq.j))then 
                  finder=.true.
               else
                  element_tester => element_tester%pntr
               end if
            end do
            if((element_tester%contents%pntin%row_ind.eq.i)&
            .and.(element_tester%contents%pntin%col_ind.eq.j)) then   
               nmb_element=element_tester%number      
            else
               nmb_element=0 
            end if
         end if
      case('block')   
         element_tester=>pmatrix% z_element_start  
         if(.not.associated( element_tester%pntr)) then 
            if((element_tester%contents%blin%row_block_ind.eq.i)&
           .and.(element_tester%contents%blin%col_block_ind.eq.j))then 
               nmb_element=element_tester%number        
            else
               nmb_element=0 
            end if
         else  
            finder=.false.
            do while((associated(element_tester%pntr))&
                     .and.(.not.finder))
               if((element_tester%contents%blin%row_block_ind.eq.i)&
           .and.(element_tester%contents%blin%col_block_ind.eq.j))then 
                  finder=.true.
               else
                  element_tester => element_tester%pntr
               end if
            end do
            if((element_tester%contents%blin%row_block_ind.eq.i)&
          .and.(element_tester%contents%blin%col_block_ind.eq.j))then   
               nmb_element=element_tester%number      
            else
               nmb_element=0 
            end if
         end if
      case('vblock')
         element_tester=>pmatrix% z_element_start  
         if(.not.associated( element_tester%pntr)) then 
            if((element_tester%contents%vblin%row_vblock_ind.eq.i)&
         .and.(element_tester%contents%vblin%col_vblock_ind.eq.j))then 
               nmb_element=element_tester%number        
            else
               nmb_element=0 
            end if
         else  
            finder=.false.
            do while((associated(element_tester%pntr)).and.&
                     (.not.finder))
            if((element_tester%contents%vblin%row_vblock_ind.eq.i)&
         .and.(element_tester%contents%vblin%col_vblock_ind.eq.j))then 
                  finder=.true.
               else
                  element_tester => element_tester%pntr
               end if
            end do
            if((element_tester%contents%vblin%row_vblock_ind.eq.i)&
         .and.(element_tester%contents%vblin%col_vblock_ind.eq.j))then   
               nmb_element=element_tester%number      
            else
               nmb_element=0 
            end if
         end if
      case default 
         istat = blas_error_param
         return
      end select
      istat = 0
      end subroutine z_element_num 
!*
      subroutine z_dealloc (nmb,istat)
      implicit none
      integer,intent(in)::nmb
      integer,intent(out)::istat
      type(z_matrix ),pointer::pmatrix
      type(z_element ),pointer ::element_tester,next_element
      istat = -1
      call zaccess_matrix (pmatrix,nmb,istat)
      if (istat.ne.0) return
      element_tester=>pmatrix% z_element_start 
      if(.not.associated(element_tester%pntr)) then
         call  dealloc_z_element (element_tester%number,&
                              pmatrix,istat)
         if (istat.ne.0) return
      else
         next_element=>element_tester%pntr
         do while((associated(next_element)))
           call  dealloc_z_element (element_tester%number,&
                                 pmatrix,istat)
            if (istat.ne.0) return
            element_tester=>next_element
            next_element=>element_tester%pntr
         end do
         call  dealloc_z_element (element_tester%number,&
                              pmatrix,istat)
         if (istat.ne.0) return
      end if
      call  dealloc_z_matrix (nmb,istat)
      if (istat.ne.0) return
      istat = 0
      return
      end subroutine z_dealloc 
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
end module mod_INSERTING

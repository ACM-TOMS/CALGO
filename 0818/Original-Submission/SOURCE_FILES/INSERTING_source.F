      subroutine NEW_X_MATRIX(nmb,Mb,ierr)
      implicit none
      integer ,intent(out)::nmb,ierr
      integer ,intent(in)::Mb
      type(X_MATRIX),pointer::matrix_insert

      if (.NOT.XINS_INIT) then
         nullify(X_MATRIX_start)
         XINS_INIT = .TRUE.
      end if

      if (.not.associated(X_MATRIX_start)) then
         allocate(X_MATRIX_start,STAT=ierr)
         X_MATRIX_start%number=XSP_MATRIX    
         X_MATRIX_start%number=-X_MATRIX_start%number    
         nullify(X_MATRIX_start%pntr)      
      else
         allocate(matrix_insert,STAT=ierr)
         matrix_insert%number=X_MATRIX_start%number-no_of_types
         matrix_insert%pntr=>X_MATRIX_start
         X_MATRIX_start=> matrix_insert    
      end if
      X_MATRIX_start%DIM=0
      X_MATRIX_start%property=blas_general+blas_one_base+blas_col_major
      X_MATRIX_start%new = 1    !new=0:blas_open_handle, new=1: blas_new_handle
      X_MATRIX_start%format=''
      nullify(X_MATRIX_start%sub_rows,X_MATRIX_start%sub_cols)
      nullify(X_MATRIX_start%X_ELEMENT_start)
      allocate(X_MATRIX_start%trb(Mb),X_MATRIX_start%tre(Mb)) 
      nmb=X_MATRIX_start%number
      ierr=0
      end subroutine NEW_X_MATRIX      
!!!
!*
!!!
      subroutine DELOC_X_MATRIX (nmb,ierr)
      implicit none
      integer ,intent(in)::nmb
      integer ,intent(out)::ierr
      type(X_MATRIX),pointer ::matrix_precedent,matrix_tester

      ierr=-1
      if(.not.associated(X_MATRIX_start%pntr)) then
         if(X_MATRIX_start%number.eq.nmb) then
            deallocate(X_MATRIX_start%tre,X_MATRIX_start%trb,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            deallocate(X_MATRIX_start,STAT=ierr)
            if (ierr.ne.0) then
               ierr = blas_error_memdeloc
               return
            end if
            nullify(X_MATRIX_start)
            ierr=0
            return
         end if
      else
         matrix_tester=>X_MATRIX_start
         if(matrix_tester%number.eq.nmb) then
            X_MATRIX_start=>matrix_tester%pntr
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
         matrix_precedent=>X_MATRIX_start
         matrix_tester=>X_MATRIX_start%pntr
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
      end subroutine DELOC_X_MATRIX
!!!
!*
!!!
      subroutine XACCESS_MATRIX(pmatrix,nmb,istat)
      implicit none
      type(X_MATRIX),pointer ::pmatrix
      integer,intent(out)::istat
      integer,intent(in) ::nmb
      type(X_MATRIX),pointer ::matrix_tester
      
      istat=-1
      matrix_tester=>X_MATRIX_start       
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
      end  subroutine XACCESS_MATRIX
!!!
!*
!!!
      subroutine NEW_X_ELEMENT(pmatrix,nmb_element,istat)
      implicit none
      type(X_MATRIX),pointer::pmatrix
      integer,intent(out)::nmb_element,istat
      type(X_ELEMENT),pointer::element_insert
      integer :: ierr

      istat = -1
      if (.not.associated(pmatrix%X_ELEMENT_start)) then
         pmatrix%new=0 !status changed to blas_open_handle     
         allocate(pmatrix%X_ELEMENT_start,STAT=ierr)
         if (ierr.ne.0) then
            istat = blas_error_memalloc
            return
         end if
         pmatrix%X_ELEMENT_start%number=1 !will certainly changed
         nullify(pmatrix%X_ELEMENT_start%pntr) 
      else
         allocate(element_insert,STAT=ierr)
         if (ierr.ne.0) then
            istat = blas_error_memalloc
            return
         end if
         element_insert%pntr=>pmatrix%X_ELEMENT_start
         element_insert%number=pmatrix%X_ELEMENT_start%number+1       
         pmatrix%X_ELEMENT_start=> element_insert 
      end if

      select case(pmatrix%format)
         case('normal')
            pmatrix%X_ELEMENT_start%contents%pntin%value=0
            pmatrix%X_ELEMENT_start%contents%pntin%row_ind=-1
            pmatrix%X_ELEMENT_start%contents%pntin%col_ind=-1 
            nullify(pmatrix%X_ELEMENT_start%contents%blin%value)
            nullify(pmatrix%X_ELEMENT_start%contents%vblin%value)
         case('block')
            nullify(pmatrix%X_ELEMENT_start%contents%blin%value)
            nullify(pmatrix%X_ELEMENT_start%contents%vblin%value)
            pmatrix%X_ELEMENT_start%contents%blin%row_block_ind=-1
            pmatrix%X_ELEMENT_start%contents%blin%col_block_ind=-1
         case('vblock')
            nullify(pmatrix%X_ELEMENT_start%contents%blin%value)
            nullify(pmatrix%X_ELEMENT_start%contents%vblin%value)
            pmatrix%X_ELEMENT_start%contents%vblin%row_vblock_ind=-1
            pmatrix%X_ELEMENT_start%contents%vblin%col_vblock_ind=-1
         case default 
            istat = blas_error_param
            return
      end select
      nmb_element=pmatrix%X_ELEMENT_start%number
      istat=0
      
      end subroutine NEW_X_ELEMENT 
!!!
!*
!!!
      subroutine DELOC_X_ELEMENT (nmb_element,pmatrix,istat)
      implicit none
      integer ,intent(in)::nmb_element
      type(X_MATRIX),pointer::pmatrix
      integer ,intent(out)::istat
      type(X_ELEMENT),pointer ::element_tester
      integer::ierr

      istat=-1
      if(.not.associated( pmatrix%X_ELEMENT_start%pntr)) then
         if(pmatrix%X_ELEMENT_start%number.eq.nmb_element) then
         if(associated(pmatrix%X_ELEMENT_start%contents%vblin%value))&
             then
             deallocate(pmatrix%X_ELEMENT_start%contents%vblin%value,&
                        STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end  if
          if(associated(pmatrix%X_ELEMENT_start%contents%blin%value))& 
            then
            deallocate(pmatrix%X_ELEMENT_start%contents%blin%value,&
                       STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if
            if (associated(pmatrix%X_ELEMENT_start)) then
               deallocate(pmatrix%X_ELEMENT_start,STAT=ierr)
               if (ierr.ne.0) then
                  istat = blas_error_memdeloc
                  return
               end if
            end if 
            nullify(pmatrix%X_ELEMENT_start)
         end if      
         istat = 0 
         return
      else 
         element_tester=>pmatrix%X_ELEMENT_start
         if(element_tester%number.eq.nmb_element) then
            pmatrix%X_ELEMENT_start=>element_tester%pntr
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
         element_tester=>pmatrix%X_ELEMENT_start%pntr
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
      end subroutine DELOC_X_ELEMENT  
!!!
!*
!!!
      subroutine XACCESS_ELEMENT(pelement,nmb_element,&
                                 pmatrix,istat)
      implicit none
      type(X_INELEMENT),pointer::pelement
      integer,intent(in) ::nmb_element
      type(X_MATRIX),pointer::pmatrix
      integer,intent(out)::istat
      type(X_ELEMENT),pointer ::element_tester

      istat=-1
      element_tester=>pmatrix%X_ELEMENT_start       
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
      end  subroutine  XACCESS_ELEMENT
!!!
!*
!!!
      subroutine X_ELEMENT_NUM (nmb_element,pmatrix,i,j,istat)
      implicit none
      integer,intent(out),target::nmb_element
      type(X_MATRIX),pointer::pmatrix
      integer ,intent(in)::i,j
      integer,intent(out)::istat
      type(X_ELEMENT),pointer ::element_tester
      logical:: finder

      istat = -1
      select case(pmatrix%format)
      case('normal')
         element_tester=>pmatrix%X_ELEMENT_start 
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
         element_tester=>pmatrix%X_ELEMENT_start 
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
         element_tester=>pmatrix%X_ELEMENT_start 
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
      end subroutine X_ELEMENT_NUM
!!!
!*
!!!
      subroutine X_DEALLOC(nmb,istat)
      implicit none
      integer,intent(in)::nmb
      integer,intent(out)::istat
      type(X_MATRIX),pointer::pmatrix
      type(X_ELEMENT),pointer ::element_tester,next_element
      
      istat = -1
      call XACCESS_MATRIX(pmatrix,nmb,istat)
      if (istat.ne.0) return
      element_tester=>pmatrix%X_ELEMENT_start
      if(.not.associated(element_tester%pntr)) then
         call  DELOC_X_ELEMENT(element_tester%number,&
                              pmatrix,istat)
         if (istat.ne.0) return
      else
         next_element=>element_tester%pntr
         do while((associated(next_element)))
           call  DELOC_X_ELEMENT(element_tester%number,&
                                 pmatrix,istat)
            if (istat.ne.0) return
            element_tester=>next_element
            next_element=>element_tester%pntr
         end do
         call  DELOC_X_ELEMENT(element_tester%number,&
                              pmatrix,istat)
         if (istat.ne.0) return
      end if
      call  DELOC_X_MATRIX(nmb,istat)
      if (istat.ne.0) return
      istat = 0
      return

      end subroutine X_DEALLOC



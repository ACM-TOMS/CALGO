      module mod_hash
      use blas_sparse_namedconstants
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 3.1.2002
!      
!     Description : A hash table for 'COO' and 'BCO' triangular solver
! **********************************************************************

      implicit none

      type capsule
      integer :: jndx,val_pos
      type(capsule), pointer :: pntr
      end type capsule
      type cappntr
      type(capsule), pointer :: pntr
      end type cappntr
      
      type(capsule), dimension(:), target, allocatable :: hash
      type(cappntr), dimension(:), allocatable :: hash_top

      contains

      subroutine setup_hash(n,ierr)
      implicit none
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: i      
      ierr = -1
      allocate(hash(n),STAT=ierr)
      if (ierr.ne.0) then
         ierr=blas_error_memalloc 
         return
      end if
      allocate(hash_top(n),STAT=ierr)
      if (ierr.ne.0) then
         ierr=blas_error_memalloc 
         return
      end if
      do i = 1,n
         nullify(hash(i)%pntr)
         hash%jndx = -1
         hash%val_pos = -1
         hash_top(i)%pntr => hash(i)
      end do
      ierr = 0
      end subroutine setup_hash

      subroutine new_capsule_main(indx,jndx,pos,ierr)
      implicit none
      integer, intent(in) :: indx,jndx,pos
      integer, intent(out) :: ierr
      type(capsule), pointer :: cap
      ierr = -1
      if ((indx.lt.lbound(hash,1)).or.(indx.gt.ubound(hash,1))) then
         return
      end if
      if(indx.eq.jndx) then
         hash(indx)%val_pos = pos
         hash(indx)%jndx = jndx
      else
         allocate(cap,STAT=ierr)
         if (ierr.ne.0) then
            ierr=blas_error_memalloc 
            return
         end if
         cap%val_pos = pos
         cap%jndx = jndx
         nullify(cap%pntr)
         hash_top(indx)%pntr%pntr => cap
         hash_top(indx)%pntr => cap
      end if
      ierr = 0
      end subroutine new_capsule_main

      subroutine print_hash()
      implicit none
      integer :: i
      type(capsule), pointer :: dummy
      do i=lbound(hash,1),ubound(hash,1)
      write(*,*)'print hash(',i,') '
      dummy => hash(i)
      do while(associated(dummy%pntr))
         write(*,*)'jndx : ', dummy%jndx
         write(*,*)'val_pos : ',dummy%val_pos
         dummy => dummy%pntr
      end do
      write(*,*)'jndx : ', dummy%jndx      
      write(*,*)'val_pos : ',dummy%val_pos
      end do
      end subroutine print_hash

      subroutine remove_hash(ierr)
      implicit none
      integer, intent(out) :: ierr
      integer :: i
      ierr = -1
      do i=lbound(hash,1),ubound(hash,1)
         do while(.not.associated(hash_top(i)%pntr,hash(i)))
            call del_capsule(i,ierr)
            if (ierr.ne.0) then
               ierr=blas_error_memdeloc 
               return
            end if         
         end do
      end do
      deallocate(hash,hash_top,STAT=ierr)
      if (ierr.ne.0) then
         ierr=blas_error_memdeloc         
         return
      end if
      end subroutine remove_hash

      subroutine del_capsule(nmb,ierr)
      implicit none
      integer, intent(in) :: nmb
      integer, intent(out) :: ierr
      type(capsule), pointer :: dummy
      dummy => hash(nmb)
      if (associated(dummy,hash_top(nmb)%pntr)) then
         ierr = -1
         return
      end if
      do while(.not.associated(dummy%pntr,hash_top(nmb)%pntr))
         dummy => dummy%pntr
      end do
      hash_top(nmb)%pntr => dummy
      deallocate(dummy%pntr,STAT=ierr)
      if(ierr.ne.0) then
        ierr=blas_error_memdeloc
        return
      end if
      end subroutine del_capsule

      end module mod_hash

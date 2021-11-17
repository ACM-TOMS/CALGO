      module mod_lsbv_bco
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'BCO'-STORAGE
!                   lsbv = Left Solve By Vector
! **********************************************************************
      use mod_hash
      use representation_of_data
      use properties
      use mod_dense_mat_algos
      implicit none
      interface lsbv_bco
        module procedure ilsbv_bco
        module procedure slsbv_bco
        module procedure dlsbv_bco
        module procedure clsbv_bco
        module procedure zlsbv_bco
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine ilsbv_bco (mat,x,ierr)
      implicit none
      type(ispmat ), pointer :: mat
      integer , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,ofs,base,nb,mb,nnz
      integer :: mm,nn,nn_sq
      character :: diag,part,store
      type(capsule), pointer :: dummy
      integer , allocatable, dimension(:) :: y !extra stor.!
      ierr = -1
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      call get_infoa(mat%INFOA,'d',mm,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',nn,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'f',mb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'g',nb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'n',nnz,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'d',diag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'f',store,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'a',part,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if ((part.ne.'U').and.(part.ne.'L')) then
         ierr = blas_error_param
         return
      end if
      if ((mat%FIDA.ne.'BCO').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
         (mm.ne.nn).or.(nb.ne.mb)) then
         ierr = blas_error_param
         return
      end if
      nn_sq = nn*nn
      allocate(y(nn),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = 0 
      call setup_hash(nb,ierr)
      if (ierr.ne.0) then
         return
      end if
      do i = 1, nnz
         call new_capsule_main(mat%IA2(i)+ofs,mat%IA1(i)+ofs,&
              (i-1)*nn_sq+1,ierr)
         if (ierr.ne.0) then
            return
         end if
      end do
      if (part.eq.'L') then
         do i = nb,1,-1
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               call block_T_mult_vec(&
      mat%A(dummy%val_pos:dummy%val_pos+nn_sq-1),&
      x((dummy%jndx-1)*nn+1:(dummy%jndx)*nn),nn,y,nn,store,ierr)
               x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_T_left_lower(&
      mat%A(hash(i)%val_pos:hash(i)%val_pos+nn_sq-1),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do i = 1,nb
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               call block_T_mult_vec(&
      mat%A(dummy%val_pos:dummy%val_pos+nn_sq-1),&
      x((dummy%jndx-1)*nn+1:(dummy%jndx)*nn),nn,y,nn,store,ierr)
               x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_T_right_upper(&
      mat%A(hash(i)%val_pos:hash(i)%val_pos+nn_sq-1),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      call remove_hash(ierr)
      end subroutine ilsbv_bco 
! **********************************************************************
! **********************************************************************
      subroutine slsbv_bco (mat,x,ierr)
      implicit none
      type(sspmat ), pointer :: mat
      real(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,ofs,base,nb,mb,nnz
      integer :: mm,nn,nn_sq
      character :: diag,part,store
      type(capsule), pointer :: dummy
      real(KIND=sp) , allocatable, dimension(:) :: y !extra stor.!
      ierr = -1
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      call get_infoa(mat%INFOA,'d',mm,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',nn,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'f',mb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'g',nb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'n',nnz,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'d',diag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'f',store,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'a',part,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if ((part.ne.'U').and.(part.ne.'L')) then
         ierr = blas_error_param
         return
      end if
      if ((mat%FIDA.ne.'BCO').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
         (mm.ne.nn).or.(nb.ne.mb)) then
         ierr = blas_error_param
         return
      end if
      nn_sq = nn*nn
      allocate(y(nn),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = 0.0e0 
      call setup_hash(nb,ierr)
      if (ierr.ne.0) then
         return
      end if
      do i = 1, nnz
         call new_capsule_main(mat%IA2(i)+ofs,mat%IA1(i)+ofs,&
              (i-1)*nn_sq+1,ierr)
         if (ierr.ne.0) then
            return
         end if
      end do
      if (part.eq.'L') then
         do i = nb,1,-1
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               call block_T_mult_vec(&
      mat%A(dummy%val_pos:dummy%val_pos+nn_sq-1),&
      x((dummy%jndx-1)*nn+1:(dummy%jndx)*nn),nn,y,nn,store,ierr)
               x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_T_left_lower(&
      mat%A(hash(i)%val_pos:hash(i)%val_pos+nn_sq-1),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do i = 1,nb
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               call block_T_mult_vec(&
      mat%A(dummy%val_pos:dummy%val_pos+nn_sq-1),&
      x((dummy%jndx-1)*nn+1:(dummy%jndx)*nn),nn,y,nn,store,ierr)
               x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_T_right_upper(&
      mat%A(hash(i)%val_pos:hash(i)%val_pos+nn_sq-1),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      call remove_hash(ierr)
      end subroutine slsbv_bco 
! **********************************************************************
! **********************************************************************
      subroutine dlsbv_bco (mat,x,ierr)
      implicit none
      type(dspmat ), pointer :: mat
      real(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,ofs,base,nb,mb,nnz
      integer :: mm,nn,nn_sq
      character :: diag,part,store
      type(capsule), pointer :: dummy
      real(KIND=dp) , allocatable, dimension(:) :: y !extra stor.!
      ierr = -1
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      call get_infoa(mat%INFOA,'d',mm,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',nn,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'f',mb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'g',nb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'n',nnz,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'d',diag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'f',store,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'a',part,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if ((part.ne.'U').and.(part.ne.'L')) then
         ierr = blas_error_param
         return
      end if
      if ((mat%FIDA.ne.'BCO').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
         (mm.ne.nn).or.(nb.ne.mb)) then
         ierr = blas_error_param
         return
      end if
      nn_sq = nn*nn
      allocate(y(nn),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = 0.0d0 
      call setup_hash(nb,ierr)
      if (ierr.ne.0) then
         return
      end if
      do i = 1, nnz
         call new_capsule_main(mat%IA2(i)+ofs,mat%IA1(i)+ofs,&
              (i-1)*nn_sq+1,ierr)
         if (ierr.ne.0) then
            return
         end if
      end do
      if (part.eq.'L') then
         do i = nb,1,-1
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               call block_T_mult_vec(&
      mat%A(dummy%val_pos:dummy%val_pos+nn_sq-1),&
      x((dummy%jndx-1)*nn+1:(dummy%jndx)*nn),nn,y,nn,store,ierr)
               x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_T_left_lower(&
      mat%A(hash(i)%val_pos:hash(i)%val_pos+nn_sq-1),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do i = 1,nb
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               call block_T_mult_vec(&
      mat%A(dummy%val_pos:dummy%val_pos+nn_sq-1),&
      x((dummy%jndx-1)*nn+1:(dummy%jndx)*nn),nn,y,nn,store,ierr)
               x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_T_right_upper(&
      mat%A(hash(i)%val_pos:hash(i)%val_pos+nn_sq-1),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      call remove_hash(ierr)
      end subroutine dlsbv_bco 
! **********************************************************************
! **********************************************************************
      subroutine clsbv_bco (mat,x,ierr)
      implicit none
      type(cspmat ), pointer :: mat
      complex(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,ofs,base,nb,mb,nnz
      integer :: mm,nn,nn_sq
      character :: diag,part,store
      type(capsule), pointer :: dummy
      complex(KIND=sp) , allocatable, dimension(:) :: y !extra stor.!
      ierr = -1
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      call get_infoa(mat%INFOA,'d',mm,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',nn,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'f',mb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'g',nb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'n',nnz,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'d',diag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'f',store,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'a',part,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if ((part.ne.'U').and.(part.ne.'L')) then
         ierr = blas_error_param
         return
      end if
      if ((mat%FIDA.ne.'BCO').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
         (mm.ne.nn).or.(nb.ne.mb)) then
         ierr = blas_error_param
         return
      end if
      nn_sq = nn*nn
      allocate(y(nn),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = (0.0e0, 0.0e0) 
      call setup_hash(nb,ierr)
      if (ierr.ne.0) then
         return
      end if
      do i = 1, nnz
         call new_capsule_main(mat%IA2(i)+ofs,mat%IA1(i)+ofs,&
              (i-1)*nn_sq+1,ierr)
         if (ierr.ne.0) then
            return
         end if
      end do
      if (part.eq.'L') then
         do i = nb,1,-1
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               call block_T_mult_vec(&
      mat%A(dummy%val_pos:dummy%val_pos+nn_sq-1),&
      x((dummy%jndx-1)*nn+1:(dummy%jndx)*nn),nn,y,nn,store,ierr)
               x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_T_left_lower(&
      mat%A(hash(i)%val_pos:hash(i)%val_pos+nn_sq-1),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do i = 1,nb
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               call block_T_mult_vec(&
      mat%A(dummy%val_pos:dummy%val_pos+nn_sq-1),&
      x((dummy%jndx-1)*nn+1:(dummy%jndx)*nn),nn,y,nn,store,ierr)
               x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_T_right_upper(&
      mat%A(hash(i)%val_pos:hash(i)%val_pos+nn_sq-1),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      call remove_hash(ierr)
      end subroutine clsbv_bco 
! **********************************************************************
! **********************************************************************
      subroutine zlsbv_bco (mat,x,ierr)
      implicit none
      type(zspmat ), pointer :: mat
      complex(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: i,n,ofs,base,nb,mb,nnz
      integer :: mm,nn,nn_sq
      character :: diag,part,store
      type(capsule), pointer :: dummy
      complex(KIND=dp) , allocatable, dimension(:) :: y !extra stor.!
      ierr = -1
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      call get_infoa(mat%INFOA,'d',mm,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'e',nn,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'f',mb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'g',nb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_infoa(mat%INFOA,'n',nnz,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'d',diag,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'f',store,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      call get_descra(mat%DESCRA,'a',part,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      if ((part.ne.'U').and.(part.ne.'L')) then
         ierr = blas_error_param
         return
      end if
      if ((mat%FIDA.ne.'BCO').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
         (mm.ne.nn).or.(nb.ne.mb)) then
         ierr = blas_error_param
         return
      end if
      nn_sq = nn*nn
      allocate(y(nn),STAT=ierr)
      if(ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      y = (0.0d0, 0.0d0) 
      call setup_hash(nb,ierr)
      if (ierr.ne.0) then
         return
      end if
      do i = 1, nnz
         call new_capsule_main(mat%IA2(i)+ofs,mat%IA1(i)+ofs,&
              (i-1)*nn_sq+1,ierr)
         if (ierr.ne.0) then
            return
         end if
      end do
      if (part.eq.'L') then
         do i = nb,1,-1
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               call block_T_mult_vec(&
      mat%A(dummy%val_pos:dummy%val_pos+nn_sq-1),&
      x((dummy%jndx-1)*nn+1:(dummy%jndx)*nn),nn,y,nn,store,ierr)
               x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_T_left_lower(&
      mat%A(hash(i)%val_pos:hash(i)%val_pos+nn_sq-1),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do i = 1,nb
            dummy => hash(i)
            do while(associated(dummy%pntr))
               dummy => dummy%pntr
               call block_T_mult_vec(&
      mat%A(dummy%val_pos:dummy%val_pos+nn_sq-1),&
      x((dummy%jndx-1)*nn+1:(dummy%jndx)*nn),nn,y,nn,store,ierr)
               x((i-1)*nn+1:i*nn) = x((i-1)*nn+1:i*nn) - y
            end do
            if (diag.ne.'U') then
               if(hash(i)%jndx.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_T_right_upper(&
      mat%A(hash(i)%val_pos:hash(i)%val_pos+nn_sq-1),&
      x((i-1)*nn+1:i*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      call remove_hash(ierr)
      end subroutine zlsbv_bco 
! **********************************************************************
! **********************************************************************
      end module mod_lsbv_bco

      module mod_rsbv_bsc
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'BSC'-STORAGE
!                   rsbv = Right Solve By Vector
! **********************************************************************
      use representation_of_data
      use properties
      use mod_dense_mat_algos
      implicit none
      interface rsbv_bsc
        module procedure irsbv_bsc
        module procedure srsbv_bsc
        module procedure drsbv_bsc
        module procedure crsbv_bsc
        module procedure zrsbv_bsc
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine irsbv_bsc (mat,x,ierr)
      implicit none
      type(ispmat ), pointer :: mat
      integer , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: j,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq
      character :: diag,part,store
      integer , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      bofs = -base
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
      ierr = -1
      if ((mat%FIDA.ne.'BSC').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
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
      if (part.eq.'L') then
         do j = 1,nb
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr) 
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(j)
               dd = -1
               do while((pntr.lt.mat%pe(j)).and. &
                       (mat%IA1(pntr+ofs)+ofs.ne.j))
                  pntr = pntr + 1
               end do
               if (mat%IA1(pntr+ofs)+ofs.eq.j) then
                  dd = pntr
               end if
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_left_lower( &
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(mat%IA1(pntr + ofs) + ofs.ne.j) then
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr)
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  end if 
                  pntr = pntr + 1
               end do 
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do j = nb,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr) 
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(j)
               dd = -1
               do while((pntr.lt.mat%pe(j)).and. &
                        (mat%IA1(pntr+ofs)+ofs.ne.j))
                  pntr = pntr + 1
               end do
               if (mat%IA1(pntr+ofs)+ofs.eq.j) then
                  dd = pntr
               end if
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_right_upper( &
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(mat%IA1(pntr + ofs) + ofs.ne.j) then
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr)
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  end if 
                  pntr = pntr + 1
               end do 
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      end subroutine irsbv_bsc 
! **********************************************************************
! **********************************************************************
      subroutine srsbv_bsc (mat,x,ierr)
      implicit none
      type(sspmat ), pointer :: mat
      real(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: j,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq
      character :: diag,part,store
      real(KIND=sp) , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      bofs = -base
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
      ierr = -1
      if ((mat%FIDA.ne.'BSC').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
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
      if (part.eq.'L') then
         do j = 1,nb
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr) 
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(j)
               dd = -1
               do while((pntr.lt.mat%pe(j)).and. &
                       (mat%IA1(pntr+ofs)+ofs.ne.j))
                  pntr = pntr + 1
               end do
               if (mat%IA1(pntr+ofs)+ofs.eq.j) then
                  dd = pntr
               end if
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_left_lower( &
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(mat%IA1(pntr + ofs) + ofs.ne.j) then
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr)
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  end if 
                  pntr = pntr + 1
               end do 
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do j = nb,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr) 
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(j)
               dd = -1
               do while((pntr.lt.mat%pe(j)).and. &
                        (mat%IA1(pntr+ofs)+ofs.ne.j))
                  pntr = pntr + 1
               end do
               if (mat%IA1(pntr+ofs)+ofs.eq.j) then
                  dd = pntr
               end if
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_right_upper( &
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(mat%IA1(pntr + ofs) + ofs.ne.j) then
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr)
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  end if 
                  pntr = pntr + 1
               end do 
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      end subroutine srsbv_bsc 
! **********************************************************************
! **********************************************************************
      subroutine drsbv_bsc (mat,x,ierr)
      implicit none
      type(dspmat ), pointer :: mat
      real(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: j,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq
      character :: diag,part,store
      real(KIND=dp) , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      bofs = -base
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
      ierr = -1
      if ((mat%FIDA.ne.'BSC').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
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
      if (part.eq.'L') then
         do j = 1,nb
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr) 
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(j)
               dd = -1
               do while((pntr.lt.mat%pe(j)).and. &
                       (mat%IA1(pntr+ofs)+ofs.ne.j))
                  pntr = pntr + 1
               end do
               if (mat%IA1(pntr+ofs)+ofs.eq.j) then
                  dd = pntr
               end if
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_left_lower( &
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(mat%IA1(pntr + ofs) + ofs.ne.j) then
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr)
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  end if 
                  pntr = pntr + 1
               end do 
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do j = nb,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr) 
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(j)
               dd = -1
               do while((pntr.lt.mat%pe(j)).and. &
                        (mat%IA1(pntr+ofs)+ofs.ne.j))
                  pntr = pntr + 1
               end do
               if (mat%IA1(pntr+ofs)+ofs.eq.j) then
                  dd = pntr
               end if
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_right_upper( &
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(mat%IA1(pntr + ofs) + ofs.ne.j) then
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr)
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  end if 
                  pntr = pntr + 1
               end do 
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      end subroutine drsbv_bsc 
! **********************************************************************
! **********************************************************************
      subroutine crsbv_bsc (mat,x,ierr)
      implicit none
      type(cspmat ), pointer :: mat
      complex(KIND=sp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: j,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq
      character :: diag,part,store
      complex(KIND=sp) , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      bofs = -base
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
      ierr = -1
      if ((mat%FIDA.ne.'BSC').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
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
      if (part.eq.'L') then
         do j = 1,nb
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr) 
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(j)
               dd = -1
               do while((pntr.lt.mat%pe(j)).and. &
                       (mat%IA1(pntr+ofs)+ofs.ne.j))
                  pntr = pntr + 1
               end do
               if (mat%IA1(pntr+ofs)+ofs.eq.j) then
                  dd = pntr
               end if
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_left_lower( &
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(mat%IA1(pntr + ofs) + ofs.ne.j) then
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr)
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  end if 
                  pntr = pntr + 1
               end do 
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do j = nb,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr) 
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(j)
               dd = -1
               do while((pntr.lt.mat%pe(j)).and. &
                        (mat%IA1(pntr+ofs)+ofs.ne.j))
                  pntr = pntr + 1
               end do
               if (mat%IA1(pntr+ofs)+ofs.eq.j) then
                  dd = pntr
               end if
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_right_upper( &
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(mat%IA1(pntr + ofs) + ofs.ne.j) then
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr)
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  end if 
                  pntr = pntr + 1
               end do 
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      end subroutine crsbv_bsc 
! **********************************************************************
! **********************************************************************
      subroutine zrsbv_bsc (mat,x,ierr)
      implicit none
      type(zspmat ), pointer :: mat
      complex(KIND=dp) , dimension(:), intent(inout) :: x
      integer, intent(out) :: ierr
      integer :: j,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq
      character :: diag,part,store
      complex(KIND=dp) , allocatable, dimension(:) :: y
      ierr = -1
      n = size(x)
      call get_infoa(mat%INFOA,'b',base,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      ofs = 1 - base
      bofs = -base
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
      ierr = -1
      if ((mat%FIDA.ne.'BSC').or.(mat%M.ne.n).or.(mat%K.ne.n).or.&
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
      if (part.eq.'L') then
         do j = 1,nb
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr) 
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(j)
               dd = -1
               do while((pntr.lt.mat%pe(j)).and. &
                       (mat%IA1(pntr+ofs)+ofs.ne.j))
                  pntr = pntr + 1
               end do
               if (mat%IA1(pntr+ofs)+ofs.eq.j) then
                  dd = pntr
               end if
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_left_lower( &
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(mat%IA1(pntr + ofs) + ofs.ne.j) then
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr)
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  end if 
                  pntr = pntr + 1
               end do 
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      else 
         do j = nb,1,-1
            if (diag.eq.'U') then
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq),&
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr) 
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  pntr = pntr + 1
               end do 
            else
               pntr = mat%pb(j)
               dd = -1
               do while((pntr.lt.mat%pe(j)).and. &
                        (mat%IA1(pntr+ofs)+ofs.ne.j))
                  pntr = pntr + 1
               end do
               if (mat%IA1(pntr+ofs)+ofs.eq.j) then
                  dd = pntr
               end if
               if(dd.eq.-1) then
                  ierr = blas_error_singtria
                  return
               else
                  call invert_right_upper( &
      mat%A((dd + bofs)*nn_sq+1:(dd + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,store,ierr)
               end if
               if(ierr.ne.0) then
                  ierr = blas_error_singtria
                  return
               end if
               pntr = mat%pb(j)
               do while(pntr.lt.mat%pe(j))
                  if(mat%IA1(pntr + ofs) + ofs.ne.j) then
                  call block_mult_vec( &
      mat%A((pntr + bofs)*nn_sq+1:(pntr + bofs +1)*nn_sq), &
      x((j-1)*nn+1:j*nn),nn,y,nn,store,ierr)
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)=&
      x((mat%IA1(pntr+ofs)+bofs)*nn+1:(mat%IA1(pntr+ofs)+bofs+1)*nn)-y
                  end if 
                  pntr = pntr + 1
               end do 
            end if
         end do
         deallocate(y,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memdeloc
            return
         end if
      end if
      end subroutine zrsbv_bsc 
! **********************************************************************
! **********************************************************************
      end module mod_rsbv_bsc

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
      subroutine iINS_entry  (pmatrix,val,i,j,istat)
      implicit none
      type(i_matrix ),pointer ::pmatrix
      integer  ,intent(in) ::val
      integer ,intent(in) ::i,j
      integer, intent(out) :: istat
      type(i_inelement ),pointer ::pelement
      integer::nmb_element,ind
      istat=-1
      if((i.gt.pmatrix%DIM(1)).or.&
           (j.gt.pmatrix%DIM(2))) then
         istat = blas_error_param 
         return 
      else
         call new_i_element (pmatrix,nmb_element,istat)
         if (istat.ne.0) return
         call i_element_num (ind,pmatrix,i,j,istat)
         if (istat.ne.0) return
         if(ind.eq.0) then
            call access_element(pelement,nmb_element,&
                                pmatrix,istat)
            if (istat.ne.0) return
            pelement%pntin%value=val
            pelement%pntin%row_ind=i
            pelement%pntin%col_ind=j
         else 
            call access_element(pelement,ind,pmatrix,istat)
            if (istat.ne.0) return
            pelement%pntin%value= pelement%pntin%value+val
            call  dealloc_i_element (nmb_element,pmatrix,istat)
            if (istat.ne.0) return
         end if
      end if
      end subroutine iINS_entry  
!*
      subroutine iINS_block (pmatrix,val,i,j,istat)
      implicit none
      type( i_matrix ),pointer ::pmatrix
      integer  ,dimension(:,:) ,target,intent(in)::val
      integer ,intent(in)::i,j
      integer,intent(out)::istat
      integer ::nmb_element,ind,ierr
      integer   ,dimension(:,:),allocatable,target::vv
      type(i_inelement ),pointer::pelement
      integer ::s_rows,s_cols
      istat = -1
      s_rows=size(val,1)
      s_cols=size(val,2)
      allocate(vv(s_rows,s_cols),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif
      if((i.gt.pmatrix%DIM(3).or.(j.gt.pmatrix%DIM(4))&
         .or.(s_rows.ne.pmatrix%DIM(5)&
         .or.(s_cols.ne.pmatrix%DIM(6))))) then 
         istat = blas_error_param
         return 
      else
         call new_i_element (pmatrix,nmb_element,istat) 
         if (istat.ne.0) return
         call i_element_num  (ind,pmatrix,i,j,istat)
         if (istat.ne.0) return
         if(ind.eq.0) then
            call access_element(pelement,nmb_element,&
                                pmatrix,istat)
            if (istat.ne.0) return
            allocate(pelement%blin%value(s_rows,s_cols),&
                     STAT=ierr)
            if(ierr.ne.0) then
               istat = blas_error_memalloc
               return
            endif
            pelement%blin%value=val
            pelement%blin%row_block_ind=i
            pelement%blin%col_block_ind=j    
         else
            call access_element(pelement,ind,pmatrix,istat)     
            if (istat.ne.0) return
            vv=vv+val 
            pelement%blin%value=pelement%blin%value+val
            call dealloc_i_element (nmb_element,pmatrix,istat)
            if (istat.ne.0) return
         end if
      end if
      deallocate(vv,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif
      end subroutine iINS_block 
!*
      subroutine iINS_bl_entr (pmatrix,val,i,j,istat)
      implicit none
      type(i_matrix ),pointer ::pmatrix
      integer   ,intent(in)::val
      integer,intent(in)::i,j
      integer,intent(out)::istat
      integer  ,dimension(:,:),allocatable,target ::vall
      integer::ii,jj,ierr
      istat = -1
      ii=floor(real((i-1)/(pmatrix%DIM(5))))
      jj=floor(real((j-1)/(pmatrix%DIM(6))))
      allocate(vall(pmatrix%DIM(5),pmatrix%DIM(6)),&
               STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      vall= 0 
      vall(i-ii*pmatrix%DIM(5),j-jj*pmatrix%DIM(6))=val
      call iINS_block (pmatrix,vall,ii+1,jj+1,istat)
      if (istat.ne.0) return     
      deallocate(vall,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif
      end subroutine iINS_bl_entr 
!*
      subroutine iINS_varblock (vpmatrix,val,i,j,istat)
      implicit none
      type(i_matrix ),pointer ::vpmatrix
      integer  ,dimension(:,:) ,target,intent(in)::val
      integer ,intent(in)::i,j
      integer,intent(out)::istat
      integer ::nmb_element,ind
      integer::ierr
      integer  ,dimension(:,:),allocatable,target::vv
      type(i_inelement  ),pointer::pelement
      integer ::s_rows,s_cols,k
      istat = -1
      s_rows=size(val,1)
      s_cols=size(val,2)
      allocate(vv(s_rows,s_cols),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      if((i.gt.vpmatrix%DIM(3).or.j.gt.vpmatrix%DIM(4))&
        .or.(s_rows.ne.vpmatrix%sub_rows(i))&
        .or.(s_cols.ne.vpmatrix%sub_cols(j))) then 
         istat = blas_error_param
         return 
      else
         call new_i_element  (vpmatrix,nmb_element,istat) 
         if (istat.ne.0) return     
         call i_element_num  (ind,vpmatrix,i,j,istat)
         if (istat.ne.0) return     
         if(ind.eq.0) then
            call access_element(pelement,nmb_element,&
                                vpmatrix,istat)
            if (istat.ne.0) return     
            allocate(pelement%vblin%value(vpmatrix%sub_rows(i),&
                                vpmatrix%sub_cols(j)),STAT=ierr)
            if(ierr.ne.0) then
               istat = blas_error_memalloc
               return
            endif   
            pelement%vblin%value=val
            pelement%vblin%row_vblock_ind=i
            pelement%vblin%col_vblock_ind=j
            do k=i+1,vpmatrix%DIM(3)
               vpmatrix%trb(k)=vpmatrix%trb(k)+1
            end do
            do k=1,vpmatrix%DIM(3)-1
               vpmatrix%tre(k)=vpmatrix%trb(k+1)
            end do
            vpmatrix%tre(vpmatrix%DIM(3))=nmb_element+1
         else
            call access_element(pelement,ind,vpmatrix,istat)
            if (istat.ne.0) return     
            pelement%vblin%value= pelement%vblin%value+val
            call dealloc_i_element (nmb_element,vpmatrix,istat)
            if (istat.ne.0) return     
         end if
      end if
      istat = 0
      return
      end subroutine iINS_varblock 
!*
      subroutine iINS_varbl_entr (vpmatrix,val,i,j,istat)
      implicit none
      type(i_matrix ),pointer ::vpmatrix
      integer ,intent(in)::val
      integer,intent(in)::i,j
      integer,intent(out)::istat
      integer ,dimension(:,:),allocatable ::vall
      integer::ii,jj,k,p,ind1,ind2,vall_ind1,vall_ind2,ierr
      ! determine the row of block entring
      ind1=0
      do k=1,vpmatrix%DIM(3)
         ind1=ind1+vpmatrix%sub_rows(k)
         if(ind1.ge.i) exit
      end do
      if(k.le.vpmatrix%DIM(3)) then
         ii=k
         vall_ind1=i-(ind1-vpmatrix%sub_rows(k))
      else
         istat = blas_error_param
         return
      end if
      ! determine the col of block entring
      ind2=0
      do p=1,vpmatrix%DIM(3)
         ind2=ind2+vpmatrix%sub_cols(p)
         if(ind2.ge.j) exit
      end do
      if(p.le.vpmatrix%DIM(4)) then
         jj=p
         vall_ind2=j-(ind2-vpmatrix%sub_cols(p))
      else
         istat = blas_error_param
         return
      end if
      allocate(vall(vpmatrix%sub_rows(ii),&
                    vpmatrix%sub_cols(jj)),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      vall= 0 
      vall(vall_ind1,vall_ind2)=val    
      call  iINS_varblock (vpmatrix,vall,ii,jj,istat)
      if (istat.ne.0) return     
      deallocate(vall,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end subroutine iINS_varbl_entr 
!*
      subroutine iuscr_varend (a,prpty,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer,intent(inout)::prpty
      integer:: m,n,ierr,ind,kb,mb
      integer, dimension(:),allocatable :: bindx,indx,rpntr,&
                                        cpntr,bpntrb,bpntre
      integer , dimension(:),allocatable :: val
      integer :: size_val,val_ind,indx_ind,bindx_ind,ii,jj,i,j
      type(i_matrix ),pointer::pmatrix
      type(i_inelement ),pointer::pelement
      istat=-1
      call access_matrix(pmatrix,a,istat)
      if (istat.ne.0) return     
      mb=pmatrix%DIM(3)
      kb=pmatrix%DIM(4)
      ! determine  size of val,m,n 
      size_val=0
      m=0
      n=0
      do i=1,pmatrix%DIM(3)
         do j=1,pmatrix%DIM(4)
            call i_element_num (ind,pmatrix,i,j,ierr)
            if(ind.ne.0) then
               call access_element(pelement,ind,pmatrix,istat)
               if (istat.ne.0) return     
               size_val=size_val+&
                      pmatrix%sub_rows(i)*pmatrix%sub_cols(j)
            end if
         end do
      end do
      do i=1,pmatrix%DIM(3)
         m=m+pmatrix%sub_rows(i)
      end do
      do j=1,pmatrix%DIM(4)     
         n=n+pmatrix%sub_cols(j)
      end do
      allocate(val(size_val),&
              indx(pmatrix% i_element_start %number+1),&
              bindx(pmatrix% i_element_start %number),&
              STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      val= 0 
      ! fill val ,indx and bindx
      val_ind=0
      indx_ind=1
      bindx_ind=0
      indx(1)=1
      do i=1,pmatrix%DIM(3)
         do j=1,pmatrix%DIM(4) 
            call i_element_num (ind,pmatrix,i,j,istat)
            if (istat.ne.0) return     
            if(ind.ne.0) then
               call access_element(pelement,ind,pmatrix,istat)
               if (istat.ne.0) return     
               do jj=1,pmatrix%sub_cols(j)
                  do ii=1,pmatrix%sub_rows(i)
                     val_ind=val_ind+1
                     val(val_ind)=pelement%vblin%value(ii,jj)  
                  end do
               end do
               bindx_ind=bindx_ind+1
               bindx(bindx_ind)=j
               indx_ind=indx_ind+1
               indx(indx_ind)=indx(indx_ind-1)&
                      +pmatrix%sub_rows(i)*pmatrix%sub_cols(j)
            end if
         end do
      end do
      ! fill rpntr, cpntr,bpntrb,bpntre
      allocate(rpntr(pmatrix%DIM(3)+1),&
              cpntr(pmatrix%DIM(4)+1),&
              bpntrb(pmatrix%DIM(3)),&
              bpntre(pmatrix%DIM(3)),&
              STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      rpntr(1)=1
      do i=2,pmatrix%DIM(3)+1
         rpntr(i)= rpntr(i-1)+pmatrix%sub_rows(i-1)
      end do
      cpntr(1)=1
      do j=2,pmatrix%DIM(4)+1
         cpntr(j)= cpntr(j-1)+pmatrix%sub_cols(j-1)
      end do
      do i=1,pmatrix%DIM(3)
         bpntrb(i)=pmatrix%trb(i)
         bpntre(i)=pmatrix%tre(i)
      end do
      ! RELEASING 
      call  i_dealloc (a,istat)
      if (istat.ne.0) return     
      ! CREATING MATRIX IN VBR FORMAT
      istat=-1 !needed to create copy of data
      call iuscr_vbr (m,n,val,indx,bindx,rpntr,cpntr,&
                    bpntrb,bpntre,mb,kb,prpty,istat,a)
      if (istat.ne.0) return     
      deallocate(val,bindx,indx,rpntr, &
                cpntr,bpntrb,bpntre,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end subroutine  iuscr_varend  
!*
      subroutine iuscr_normend  (a,prpty,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer,intent(inout)::prpty
      integer:: m,n,nnz
      integer, dimension(:),allocatable :: indx,jndx
      integer , dimension(:),allocatable :: val 
      integer :: nmb_element,i,ierr
      type(i_matrix  ),pointer::pmatrix
      type(i_inelement ),pointer::pelement
      call access_matrix(pmatrix,a,istat)
      if (istat.ne.0) return     
      m=pmatrix%DIM(1)  !nb_of_rows
      n=pmatrix%DIM(2)  !nb_of_cols
      nnz=pmatrix% i_element_start %number
      allocate(val(nnz),indx(nnz),jndx(nnz),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      nmb_element=nnz+1
      do i=1,nnz
         call access_element(pelement,nmb_element-i,&
                             pmatrix,istat)
         if (istat.ne.0) return     
         val(i)=pelement%pntin%value
         indx(i)=pelement%pntin%row_ind
         jndx(i)=pelement%pntin%col_ind
      end do
      call  i_dealloc (a,istat)
      if (istat.ne.0) return     
      ! CREATING A MATRIX IN COO FORMAT
      istat=-1 !needed to create copy of data
      call iuscr_coo (m,n,val,indx,jndx,nnz,prpty,istat,a) 
      if (istat.ne.0) return     
      deallocate(val,indx,jndx,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end  subroutine iuscr_normend  
!*
      subroutine iuscr_blockend  (a,prpty,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer,intent(inout)::prpty
      integer:: m,n,bnnz,ierr,ind,kb,lb,dummy,k,mb,i,j
      integer, dimension(:),allocatable :: bindx,bjndx
      integer  , dimension(:),allocatable :: val
      integer :: nmb_block
      type(i_matrix ),pointer::pmatrix
      type(i_inelement ),pointer::pelement
      istat=-1
      call access_matrix(pmatrix,a,istat)
      if (istat.ne.0) return     
      lb=pmatrix%DIM(5)
      bnnz=pmatrix% i_element_start %number
      mb=pmatrix%DIM(3)
      kb=pmatrix%DIM(4)
      m=mb*lb
      n=kb*lb
      ind=0
      dummy=bnnz*lb*lb
      allocate(val(dummy),bindx(bnnz),bjndx(bnnz),&
               STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      nmb_block=bnnz+1
      do i=1,bnnz
         k=1
         call access_element(pelement,nmb_block-i,pmatrix,&
                             istat)
         if (istat.ne.0) return     
         do j=1,lb
            do k=1,lb
               ind=ind+1
               val(ind)=pelement%blin%value(k,j)          
            end do
         end do
         bindx(i)=pelement%blin%row_block_ind
         bjndx(i)=pelement%blin%col_block_ind
      end do
      ! RELEASING
      call i_dealloc (a,istat)
      if (istat.ne.0) return     
      ! CREATE A MATRIX IN BCO FORMAT
      istat=-1 !needed to create copy of data
      call  iuscr_bco (m,n,val,bindx,bjndx,bnnz,&
                      mb,kb,lb,prpty,istat,a)
      if (istat.ne.0) return     
      deallocate(val,bindx,bjndx,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end subroutine  iuscr_blockend  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine sINS_entry  (pmatrix,val,i,j,istat)
      implicit none
      type(s_matrix ),pointer ::pmatrix
      real(KIND=sp)  ,intent(in) ::val
      integer ,intent(in) ::i,j
      integer, intent(out) :: istat
      type(s_inelement ),pointer ::pelement
      integer::nmb_element,ind
      istat=-1
      if((i.gt.pmatrix%DIM(1)).or.&
           (j.gt.pmatrix%DIM(2))) then
         istat = blas_error_param 
         return 
      else
         call new_s_element (pmatrix,nmb_element,istat)
         if (istat.ne.0) return
         call s_element_num (ind,pmatrix,i,j,istat)
         if (istat.ne.0) return
         if(ind.eq.0) then
            call access_element(pelement,nmb_element,&
                                pmatrix,istat)
            if (istat.ne.0) return
            pelement%pntin%value=val
            pelement%pntin%row_ind=i
            pelement%pntin%col_ind=j
         else 
            call access_element(pelement,ind,pmatrix,istat)
            if (istat.ne.0) return
            pelement%pntin%value= pelement%pntin%value+val
            call  dealloc_s_element (nmb_element,pmatrix,istat)
            if (istat.ne.0) return
         end if
      end if
      end subroutine sINS_entry  
!*
      subroutine sINS_block (pmatrix,val,i,j,istat)
      implicit none
      type( s_matrix ),pointer ::pmatrix
      real(KIND=sp)  ,dimension(:,:) ,target,intent(in)::val
      integer ,intent(in)::i,j
      integer,intent(out)::istat
      integer ::nmb_element,ind,ierr
      real(KIND=sp)   ,dimension(:,:),allocatable,target::vv
      type(s_inelement ),pointer::pelement
      integer ::s_rows,s_cols
      istat = -1
      s_rows=size(val,1)
      s_cols=size(val,2)
      allocate(vv(s_rows,s_cols),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif
      if((i.gt.pmatrix%DIM(3).or.(j.gt.pmatrix%DIM(4))&
         .or.(s_rows.ne.pmatrix%DIM(5)&
         .or.(s_cols.ne.pmatrix%DIM(6))))) then 
         istat = blas_error_param
         return 
      else
         call new_s_element (pmatrix,nmb_element,istat) 
         if (istat.ne.0) return
         call s_element_num  (ind,pmatrix,i,j,istat)
         if (istat.ne.0) return
         if(ind.eq.0) then
            call access_element(pelement,nmb_element,&
                                pmatrix,istat)
            if (istat.ne.0) return
            allocate(pelement%blin%value(s_rows,s_cols),&
                     STAT=ierr)
            if(ierr.ne.0) then
               istat = blas_error_memalloc
               return
            endif
            pelement%blin%value=val
            pelement%blin%row_block_ind=i
            pelement%blin%col_block_ind=j    
         else
            call access_element(pelement,ind,pmatrix,istat)     
            if (istat.ne.0) return
            vv=vv+val 
            pelement%blin%value=pelement%blin%value+val
            call dealloc_s_element (nmb_element,pmatrix,istat)
            if (istat.ne.0) return
         end if
      end if
      deallocate(vv,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif
      end subroutine sINS_block 
!*
      subroutine sINS_bl_entr (pmatrix,val,i,j,istat)
      implicit none
      type(s_matrix ),pointer ::pmatrix
      real(KIND=sp)   ,intent(in)::val
      integer,intent(in)::i,j
      integer,intent(out)::istat
      real(KIND=sp)  ,dimension(:,:),allocatable,target ::vall
      integer::ii,jj,ierr
      istat = -1
      ii=floor(real((i-1)/(pmatrix%DIM(5))))
      jj=floor(real((j-1)/(pmatrix%DIM(6))))
      allocate(vall(pmatrix%DIM(5),pmatrix%DIM(6)),&
               STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      vall= 0.0e0 
      vall(i-ii*pmatrix%DIM(5),j-jj*pmatrix%DIM(6))=val
      call sINS_block (pmatrix,vall,ii+1,jj+1,istat)
      if (istat.ne.0) return     
      deallocate(vall,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif
      end subroutine sINS_bl_entr 
!*
      subroutine sINS_varblock (vpmatrix,val,i,j,istat)
      implicit none
      type(s_matrix ),pointer ::vpmatrix
      real(KIND=sp)  ,dimension(:,:) ,target,intent(in)::val
      integer ,intent(in)::i,j
      integer,intent(out)::istat
      integer ::nmb_element,ind
      integer::ierr
      real(KIND=sp)  ,dimension(:,:),allocatable,target::vv
      type(s_inelement  ),pointer::pelement
      integer ::s_rows,s_cols,k
      istat = -1
      s_rows=size(val,1)
      s_cols=size(val,2)
      allocate(vv(s_rows,s_cols),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      if((i.gt.vpmatrix%DIM(3).or.j.gt.vpmatrix%DIM(4))&
        .or.(s_rows.ne.vpmatrix%sub_rows(i))&
        .or.(s_cols.ne.vpmatrix%sub_cols(j))) then 
         istat = blas_error_param
         return 
      else
         call new_s_element  (vpmatrix,nmb_element,istat) 
         if (istat.ne.0) return     
         call s_element_num  (ind,vpmatrix,i,j,istat)
         if (istat.ne.0) return     
         if(ind.eq.0) then
            call access_element(pelement,nmb_element,&
                                vpmatrix,istat)
            if (istat.ne.0) return     
            allocate(pelement%vblin%value(vpmatrix%sub_rows(i),&
                                vpmatrix%sub_cols(j)),STAT=ierr)
            if(ierr.ne.0) then
               istat = blas_error_memalloc
               return
            endif   
            pelement%vblin%value=val
            pelement%vblin%row_vblock_ind=i
            pelement%vblin%col_vblock_ind=j
            do k=i+1,vpmatrix%DIM(3)
               vpmatrix%trb(k)=vpmatrix%trb(k)+1
            end do
            do k=1,vpmatrix%DIM(3)-1
               vpmatrix%tre(k)=vpmatrix%trb(k+1)
            end do
            vpmatrix%tre(vpmatrix%DIM(3))=nmb_element+1
         else
            call access_element(pelement,ind,vpmatrix,istat)
            if (istat.ne.0) return     
            pelement%vblin%value= pelement%vblin%value+val
            call dealloc_s_element (nmb_element,vpmatrix,istat)
            if (istat.ne.0) return     
         end if
      end if
      istat = 0
      return
      end subroutine sINS_varblock 
!*
      subroutine sINS_varbl_entr (vpmatrix,val,i,j,istat)
      implicit none
      type(s_matrix ),pointer ::vpmatrix
      real(KIND=sp) ,intent(in)::val
      integer,intent(in)::i,j
      integer,intent(out)::istat
      real(KIND=sp) ,dimension(:,:),allocatable ::vall
      integer::ii,jj,k,p,ind1,ind2,vall_ind1,vall_ind2,ierr
      ! determine the row of block entring
      ind1=0
      do k=1,vpmatrix%DIM(3)
         ind1=ind1+vpmatrix%sub_rows(k)
         if(ind1.ge.i) exit
      end do
      if(k.le.vpmatrix%DIM(3)) then
         ii=k
         vall_ind1=i-(ind1-vpmatrix%sub_rows(k))
      else
         istat = blas_error_param
         return
      end if
      ! determine the col of block entring
      ind2=0
      do p=1,vpmatrix%DIM(3)
         ind2=ind2+vpmatrix%sub_cols(p)
         if(ind2.ge.j) exit
      end do
      if(p.le.vpmatrix%DIM(4)) then
         jj=p
         vall_ind2=j-(ind2-vpmatrix%sub_cols(p))
      else
         istat = blas_error_param
         return
      end if
      allocate(vall(vpmatrix%sub_rows(ii),&
                    vpmatrix%sub_cols(jj)),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      vall= 0.0e0 
      vall(vall_ind1,vall_ind2)=val    
      call  sINS_varblock (vpmatrix,vall,ii,jj,istat)
      if (istat.ne.0) return     
      deallocate(vall,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end subroutine sINS_varbl_entr 
!*
      subroutine suscr_varend (a,prpty,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer,intent(inout)::prpty
      integer:: m,n,ierr,ind,kb,mb
      integer, dimension(:),allocatable :: bindx,indx,rpntr,&
                                        cpntr,bpntrb,bpntre
      real(KIND=sp) , dimension(:),allocatable :: val
      integer :: size_val,val_ind,indx_ind,bindx_ind,ii,jj,i,j
      type(s_matrix ),pointer::pmatrix
      type(s_inelement ),pointer::pelement
      istat=-1
      call access_matrix(pmatrix,a,istat)
      if (istat.ne.0) return     
      mb=pmatrix%DIM(3)
      kb=pmatrix%DIM(4)
      ! determine  size of val,m,n 
      size_val=0
      m=0
      n=0
      do i=1,pmatrix%DIM(3)
         do j=1,pmatrix%DIM(4)
            call s_element_num (ind,pmatrix,i,j,ierr)
            if(ind.ne.0) then
               call access_element(pelement,ind,pmatrix,istat)
               if (istat.ne.0) return     
               size_val=size_val+&
                      pmatrix%sub_rows(i)*pmatrix%sub_cols(j)
            end if
         end do
      end do
      do i=1,pmatrix%DIM(3)
         m=m+pmatrix%sub_rows(i)
      end do
      do j=1,pmatrix%DIM(4)     
         n=n+pmatrix%sub_cols(j)
      end do
      allocate(val(size_val),&
              indx(pmatrix% s_element_start %number+1),&
              bindx(pmatrix% s_element_start %number),&
              STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      val= 0.0e0 
      ! fill val ,indx and bindx
      val_ind=0
      indx_ind=1
      bindx_ind=0
      indx(1)=1
      do i=1,pmatrix%DIM(3)
         do j=1,pmatrix%DIM(4) 
            call s_element_num (ind,pmatrix,i,j,istat)
            if (istat.ne.0) return     
            if(ind.ne.0) then
               call access_element(pelement,ind,pmatrix,istat)
               if (istat.ne.0) return     
               do jj=1,pmatrix%sub_cols(j)
                  do ii=1,pmatrix%sub_rows(i)
                     val_ind=val_ind+1
                     val(val_ind)=pelement%vblin%value(ii,jj)  
                  end do
               end do
               bindx_ind=bindx_ind+1
               bindx(bindx_ind)=j
               indx_ind=indx_ind+1
               indx(indx_ind)=indx(indx_ind-1)&
                      +pmatrix%sub_rows(i)*pmatrix%sub_cols(j)
            end if
         end do
      end do
      ! fill rpntr, cpntr,bpntrb,bpntre
      allocate(rpntr(pmatrix%DIM(3)+1),&
              cpntr(pmatrix%DIM(4)+1),&
              bpntrb(pmatrix%DIM(3)),&
              bpntre(pmatrix%DIM(3)),&
              STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      rpntr(1)=1
      do i=2,pmatrix%DIM(3)+1
         rpntr(i)= rpntr(i-1)+pmatrix%sub_rows(i-1)
      end do
      cpntr(1)=1
      do j=2,pmatrix%DIM(4)+1
         cpntr(j)= cpntr(j-1)+pmatrix%sub_cols(j-1)
      end do
      do i=1,pmatrix%DIM(3)
         bpntrb(i)=pmatrix%trb(i)
         bpntre(i)=pmatrix%tre(i)
      end do
      ! RELEASING 
      call  s_dealloc (a,istat)
      if (istat.ne.0) return     
      ! CREATING MATRIX IN VBR FORMAT
      istat=-1 !needed to create copy of data
      call suscr_vbr (m,n,val,indx,bindx,rpntr,cpntr,&
                    bpntrb,bpntre,mb,kb,prpty,istat,a)
      if (istat.ne.0) return     
      deallocate(val,bindx,indx,rpntr, &
                cpntr,bpntrb,bpntre,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end subroutine  suscr_varend  
!*
      subroutine suscr_normend  (a,prpty,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer,intent(inout)::prpty
      integer:: m,n,nnz
      integer, dimension(:),allocatable :: indx,jndx
      real(KIND=sp) , dimension(:),allocatable :: val 
      integer :: nmb_element,i,ierr
      type(s_matrix  ),pointer::pmatrix
      type(s_inelement ),pointer::pelement
      call access_matrix(pmatrix,a,istat)
      if (istat.ne.0) return     
      m=pmatrix%DIM(1)  !nb_of_rows
      n=pmatrix%DIM(2)  !nb_of_cols
      nnz=pmatrix% s_element_start %number
      allocate(val(nnz),indx(nnz),jndx(nnz),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      nmb_element=nnz+1
      do i=1,nnz
         call access_element(pelement,nmb_element-i,&
                             pmatrix,istat)
         if (istat.ne.0) return     
         val(i)=pelement%pntin%value
         indx(i)=pelement%pntin%row_ind
         jndx(i)=pelement%pntin%col_ind
      end do
      call  s_dealloc (a,istat)
      if (istat.ne.0) return     
      ! CREATING A MATRIX IN COO FORMAT
      istat=-1 !needed to create copy of data
      call suscr_coo (m,n,val,indx,jndx,nnz,prpty,istat,a) 
      if (istat.ne.0) return     
      deallocate(val,indx,jndx,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end  subroutine suscr_normend  
!*
      subroutine suscr_blockend  (a,prpty,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer,intent(inout)::prpty
      integer:: m,n,bnnz,ierr,ind,kb,lb,dummy,k,mb,i,j
      integer, dimension(:),allocatable :: bindx,bjndx
      real(KIND=sp)  , dimension(:),allocatable :: val
      integer :: nmb_block
      type(s_matrix ),pointer::pmatrix
      type(s_inelement ),pointer::pelement
      istat=-1
      call access_matrix(pmatrix,a,istat)
      if (istat.ne.0) return     
      lb=pmatrix%DIM(5)
      bnnz=pmatrix% s_element_start %number
      mb=pmatrix%DIM(3)
      kb=pmatrix%DIM(4)
      m=mb*lb
      n=kb*lb
      ind=0
      dummy=bnnz*lb*lb
      allocate(val(dummy),bindx(bnnz),bjndx(bnnz),&
               STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      nmb_block=bnnz+1
      do i=1,bnnz
         k=1
         call access_element(pelement,nmb_block-i,pmatrix,&
                             istat)
         if (istat.ne.0) return     
         do j=1,lb
            do k=1,lb
               ind=ind+1
               val(ind)=pelement%blin%value(k,j)          
            end do
         end do
         bindx(i)=pelement%blin%row_block_ind
         bjndx(i)=pelement%blin%col_block_ind
      end do
      ! RELEASING
      call s_dealloc (a,istat)
      if (istat.ne.0) return     
      ! CREATE A MATRIX IN BCO FORMAT
      istat=-1 !needed to create copy of data
      call  suscr_bco (m,n,val,bindx,bjndx,bnnz,&
                      mb,kb,lb,prpty,istat,a)
      if (istat.ne.0) return     
      deallocate(val,bindx,bjndx,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end subroutine  suscr_blockend  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine dINS_entry  (pmatrix,val,i,j,istat)
      implicit none
      type(d_matrix ),pointer ::pmatrix
      real(KIND=dp)  ,intent(in) ::val
      integer ,intent(in) ::i,j
      integer, intent(out) :: istat
      type(d_inelement ),pointer ::pelement
      integer::nmb_element,ind
      istat=-1
      if((i.gt.pmatrix%DIM(1)).or.&
           (j.gt.pmatrix%DIM(2))) then
         istat = blas_error_param 
         return 
      else
         call new_d_element (pmatrix,nmb_element,istat)
         if (istat.ne.0) return
         call d_element_num (ind,pmatrix,i,j,istat)
         if (istat.ne.0) return
         if(ind.eq.0) then
            call access_element(pelement,nmb_element,&
                                pmatrix,istat)
            if (istat.ne.0) return
            pelement%pntin%value=val
            pelement%pntin%row_ind=i
            pelement%pntin%col_ind=j
         else 
            call access_element(pelement,ind,pmatrix,istat)
            if (istat.ne.0) return
            pelement%pntin%value= pelement%pntin%value+val
            call  dealloc_d_element (nmb_element,pmatrix,istat)
            if (istat.ne.0) return
         end if
      end if
      end subroutine dINS_entry  
!*
      subroutine dINS_block (pmatrix,val,i,j,istat)
      implicit none
      type( d_matrix ),pointer ::pmatrix
      real(KIND=dp)  ,dimension(:,:) ,target,intent(in)::val
      integer ,intent(in)::i,j
      integer,intent(out)::istat
      integer ::nmb_element,ind,ierr
      real(KIND=dp)   ,dimension(:,:),allocatable,target::vv
      type(d_inelement ),pointer::pelement
      integer ::s_rows,s_cols
      istat = -1
      s_rows=size(val,1)
      s_cols=size(val,2)
      allocate(vv(s_rows,s_cols),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif
      if((i.gt.pmatrix%DIM(3).or.(j.gt.pmatrix%DIM(4))&
         .or.(s_rows.ne.pmatrix%DIM(5)&
         .or.(s_cols.ne.pmatrix%DIM(6))))) then 
         istat = blas_error_param
         return 
      else
         call new_d_element (pmatrix,nmb_element,istat) 
         if (istat.ne.0) return
         call d_element_num  (ind,pmatrix,i,j,istat)
         if (istat.ne.0) return
         if(ind.eq.0) then
            call access_element(pelement,nmb_element,&
                                pmatrix,istat)
            if (istat.ne.0) return
            allocate(pelement%blin%value(s_rows,s_cols),&
                     STAT=ierr)
            if(ierr.ne.0) then
               istat = blas_error_memalloc
               return
            endif
            pelement%blin%value=val
            pelement%blin%row_block_ind=i
            pelement%blin%col_block_ind=j    
         else
            call access_element(pelement,ind,pmatrix,istat)     
            if (istat.ne.0) return
            vv=vv+val 
            pelement%blin%value=pelement%blin%value+val
            call dealloc_d_element (nmb_element,pmatrix,istat)
            if (istat.ne.0) return
         end if
      end if
      deallocate(vv,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif
      end subroutine dINS_block 
!*
      subroutine dINS_bl_entr (pmatrix,val,i,j,istat)
      implicit none
      type(d_matrix ),pointer ::pmatrix
      real(KIND=dp)   ,intent(in)::val
      integer,intent(in)::i,j
      integer,intent(out)::istat
      real(KIND=dp)  ,dimension(:,:),allocatable,target ::vall
      integer::ii,jj,ierr
      istat = -1
      ii=floor(real((i-1)/(pmatrix%DIM(5))))
      jj=floor(real((j-1)/(pmatrix%DIM(6))))
      allocate(vall(pmatrix%DIM(5),pmatrix%DIM(6)),&
               STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      vall= 0.0d0 
      vall(i-ii*pmatrix%DIM(5),j-jj*pmatrix%DIM(6))=val
      call dINS_block (pmatrix,vall,ii+1,jj+1,istat)
      if (istat.ne.0) return     
      deallocate(vall,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif
      end subroutine dINS_bl_entr 
!*
      subroutine dINS_varblock (vpmatrix,val,i,j,istat)
      implicit none
      type(d_matrix ),pointer ::vpmatrix
      real(KIND=dp)  ,dimension(:,:) ,target,intent(in)::val
      integer ,intent(in)::i,j
      integer,intent(out)::istat
      integer ::nmb_element,ind
      integer::ierr
      real(KIND=dp)  ,dimension(:,:),allocatable,target::vv
      type(d_inelement  ),pointer::pelement
      integer ::s_rows,s_cols,k
      istat = -1
      s_rows=size(val,1)
      s_cols=size(val,2)
      allocate(vv(s_rows,s_cols),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      if((i.gt.vpmatrix%DIM(3).or.j.gt.vpmatrix%DIM(4))&
        .or.(s_rows.ne.vpmatrix%sub_rows(i))&
        .or.(s_cols.ne.vpmatrix%sub_cols(j))) then 
         istat = blas_error_param
         return 
      else
         call new_d_element  (vpmatrix,nmb_element,istat) 
         if (istat.ne.0) return     
         call d_element_num  (ind,vpmatrix,i,j,istat)
         if (istat.ne.0) return     
         if(ind.eq.0) then
            call access_element(pelement,nmb_element,&
                                vpmatrix,istat)
            if (istat.ne.0) return     
            allocate(pelement%vblin%value(vpmatrix%sub_rows(i),&
                                vpmatrix%sub_cols(j)),STAT=ierr)
            if(ierr.ne.0) then
               istat = blas_error_memalloc
               return
            endif   
            pelement%vblin%value=val
            pelement%vblin%row_vblock_ind=i
            pelement%vblin%col_vblock_ind=j
            do k=i+1,vpmatrix%DIM(3)
               vpmatrix%trb(k)=vpmatrix%trb(k)+1
            end do
            do k=1,vpmatrix%DIM(3)-1
               vpmatrix%tre(k)=vpmatrix%trb(k+1)
            end do
            vpmatrix%tre(vpmatrix%DIM(3))=nmb_element+1
         else
            call access_element(pelement,ind,vpmatrix,istat)
            if (istat.ne.0) return     
            pelement%vblin%value= pelement%vblin%value+val
            call dealloc_d_element (nmb_element,vpmatrix,istat)
            if (istat.ne.0) return     
         end if
      end if
      istat = 0
      return
      end subroutine dINS_varblock 
!*
      subroutine dINS_varbl_entr (vpmatrix,val,i,j,istat)
      implicit none
      type(d_matrix ),pointer ::vpmatrix
      real(KIND=dp) ,intent(in)::val
      integer,intent(in)::i,j
      integer,intent(out)::istat
      real(KIND=dp) ,dimension(:,:),allocatable ::vall
      integer::ii,jj,k,p,ind1,ind2,vall_ind1,vall_ind2,ierr
      ! determine the row of block entring
      ind1=0
      do k=1,vpmatrix%DIM(3)
         ind1=ind1+vpmatrix%sub_rows(k)
         if(ind1.ge.i) exit
      end do
      if(k.le.vpmatrix%DIM(3)) then
         ii=k
         vall_ind1=i-(ind1-vpmatrix%sub_rows(k))
      else
         istat = blas_error_param
         return
      end if
      ! determine the col of block entring
      ind2=0
      do p=1,vpmatrix%DIM(3)
         ind2=ind2+vpmatrix%sub_cols(p)
         if(ind2.ge.j) exit
      end do
      if(p.le.vpmatrix%DIM(4)) then
         jj=p
         vall_ind2=j-(ind2-vpmatrix%sub_cols(p))
      else
         istat = blas_error_param
         return
      end if
      allocate(vall(vpmatrix%sub_rows(ii),&
                    vpmatrix%sub_cols(jj)),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      vall= 0.0d0 
      vall(vall_ind1,vall_ind2)=val    
      call  dINS_varblock (vpmatrix,vall,ii,jj,istat)
      if (istat.ne.0) return     
      deallocate(vall,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end subroutine dINS_varbl_entr 
!*
      subroutine duscr_varend (a,prpty,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer,intent(inout)::prpty
      integer:: m,n,ierr,ind,kb,mb
      integer, dimension(:),allocatable :: bindx,indx,rpntr,&
                                        cpntr,bpntrb,bpntre
      real(KIND=dp) , dimension(:),allocatable :: val
      integer :: size_val,val_ind,indx_ind,bindx_ind,ii,jj,i,j
      type(d_matrix ),pointer::pmatrix
      type(d_inelement ),pointer::pelement
      istat=-1
      call access_matrix(pmatrix,a,istat)
      if (istat.ne.0) return     
      mb=pmatrix%DIM(3)
      kb=pmatrix%DIM(4)
      ! determine  size of val,m,n 
      size_val=0
      m=0
      n=0
      do i=1,pmatrix%DIM(3)
         do j=1,pmatrix%DIM(4)
            call d_element_num (ind,pmatrix,i,j,ierr)
            if(ind.ne.0) then
               call access_element(pelement,ind,pmatrix,istat)
               if (istat.ne.0) return     
               size_val=size_val+&
                      pmatrix%sub_rows(i)*pmatrix%sub_cols(j)
            end if
         end do
      end do
      do i=1,pmatrix%DIM(3)
         m=m+pmatrix%sub_rows(i)
      end do
      do j=1,pmatrix%DIM(4)     
         n=n+pmatrix%sub_cols(j)
      end do
      allocate(val(size_val),&
              indx(pmatrix% d_element_start %number+1),&
              bindx(pmatrix% d_element_start %number),&
              STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      val= 0.0d0 
      ! fill val ,indx and bindx
      val_ind=0
      indx_ind=1
      bindx_ind=0
      indx(1)=1
      do i=1,pmatrix%DIM(3)
         do j=1,pmatrix%DIM(4) 
            call d_element_num (ind,pmatrix,i,j,istat)
            if (istat.ne.0) return     
            if(ind.ne.0) then
               call access_element(pelement,ind,pmatrix,istat)
               if (istat.ne.0) return     
               do jj=1,pmatrix%sub_cols(j)
                  do ii=1,pmatrix%sub_rows(i)
                     val_ind=val_ind+1
                     val(val_ind)=pelement%vblin%value(ii,jj)  
                  end do
               end do
               bindx_ind=bindx_ind+1
               bindx(bindx_ind)=j
               indx_ind=indx_ind+1
               indx(indx_ind)=indx(indx_ind-1)&
                      +pmatrix%sub_rows(i)*pmatrix%sub_cols(j)
            end if
         end do
      end do
      ! fill rpntr, cpntr,bpntrb,bpntre
      allocate(rpntr(pmatrix%DIM(3)+1),&
              cpntr(pmatrix%DIM(4)+1),&
              bpntrb(pmatrix%DIM(3)),&
              bpntre(pmatrix%DIM(3)),&
              STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      rpntr(1)=1
      do i=2,pmatrix%DIM(3)+1
         rpntr(i)= rpntr(i-1)+pmatrix%sub_rows(i-1)
      end do
      cpntr(1)=1
      do j=2,pmatrix%DIM(4)+1
         cpntr(j)= cpntr(j-1)+pmatrix%sub_cols(j-1)
      end do
      do i=1,pmatrix%DIM(3)
         bpntrb(i)=pmatrix%trb(i)
         bpntre(i)=pmatrix%tre(i)
      end do
      ! RELEASING 
      call  d_dealloc (a,istat)
      if (istat.ne.0) return     
      ! CREATING MATRIX IN VBR FORMAT
      istat=-1 !needed to create copy of data
      call duscr_vbr (m,n,val,indx,bindx,rpntr,cpntr,&
                    bpntrb,bpntre,mb,kb,prpty,istat,a)
      if (istat.ne.0) return     
      deallocate(val,bindx,indx,rpntr, &
                cpntr,bpntrb,bpntre,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end subroutine  duscr_varend  
!*
      subroutine duscr_normend  (a,prpty,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer,intent(inout)::prpty
      integer:: m,n,nnz
      integer, dimension(:),allocatable :: indx,jndx
      real(KIND=dp) , dimension(:),allocatable :: val 
      integer :: nmb_element,i,ierr
      type(d_matrix  ),pointer::pmatrix
      type(d_inelement ),pointer::pelement
      call access_matrix(pmatrix,a,istat)
      if (istat.ne.0) return     
      m=pmatrix%DIM(1)  !nb_of_rows
      n=pmatrix%DIM(2)  !nb_of_cols
      nnz=pmatrix% d_element_start %number
      allocate(val(nnz),indx(nnz),jndx(nnz),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      nmb_element=nnz+1
      do i=1,nnz
         call access_element(pelement,nmb_element-i,&
                             pmatrix,istat)
         if (istat.ne.0) return     
         val(i)=pelement%pntin%value
         indx(i)=pelement%pntin%row_ind
         jndx(i)=pelement%pntin%col_ind
      end do
      call  d_dealloc (a,istat)
      if (istat.ne.0) return     
      ! CREATING A MATRIX IN COO FORMAT
      istat=-1 !needed to create copy of data
      call duscr_coo (m,n,val,indx,jndx,nnz,prpty,istat,a) 
      if (istat.ne.0) return     
      deallocate(val,indx,jndx,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end  subroutine duscr_normend  
!*
      subroutine duscr_blockend  (a,prpty,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer,intent(inout)::prpty
      integer:: m,n,bnnz,ierr,ind,kb,lb,dummy,k,mb,i,j
      integer, dimension(:),allocatable :: bindx,bjndx
      real(KIND=dp)  , dimension(:),allocatable :: val
      integer :: nmb_block
      type(d_matrix ),pointer::pmatrix
      type(d_inelement ),pointer::pelement
      istat=-1
      call access_matrix(pmatrix,a,istat)
      if (istat.ne.0) return     
      lb=pmatrix%DIM(5)
      bnnz=pmatrix% d_element_start %number
      mb=pmatrix%DIM(3)
      kb=pmatrix%DIM(4)
      m=mb*lb
      n=kb*lb
      ind=0
      dummy=bnnz*lb*lb
      allocate(val(dummy),bindx(bnnz),bjndx(bnnz),&
               STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      nmb_block=bnnz+1
      do i=1,bnnz
         k=1
         call access_element(pelement,nmb_block-i,pmatrix,&
                             istat)
         if (istat.ne.0) return     
         do j=1,lb
            do k=1,lb
               ind=ind+1
               val(ind)=pelement%blin%value(k,j)          
            end do
         end do
         bindx(i)=pelement%blin%row_block_ind
         bjndx(i)=pelement%blin%col_block_ind
      end do
      ! RELEASING
      call d_dealloc (a,istat)
      if (istat.ne.0) return     
      ! CREATE A MATRIX IN BCO FORMAT
      istat=-1 !needed to create copy of data
      call  duscr_bco (m,n,val,bindx,bjndx,bnnz,&
                      mb,kb,lb,prpty,istat,a)
      if (istat.ne.0) return     
      deallocate(val,bindx,bjndx,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end subroutine  duscr_blockend  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine cINS_entry  (pmatrix,val,i,j,istat)
      implicit none
      type(c_matrix ),pointer ::pmatrix
      complex(KIND=sp)  ,intent(in) ::val
      integer ,intent(in) ::i,j
      integer, intent(out) :: istat
      type(c_inelement ),pointer ::pelement
      integer::nmb_element,ind
      istat=-1
      if((i.gt.pmatrix%DIM(1)).or.&
           (j.gt.pmatrix%DIM(2))) then
         istat = blas_error_param 
         return 
      else
         call new_c_element (pmatrix,nmb_element,istat)
         if (istat.ne.0) return
         call c_element_num (ind,pmatrix,i,j,istat)
         if (istat.ne.0) return
         if(ind.eq.0) then
            call access_element(pelement,nmb_element,&
                                pmatrix,istat)
            if (istat.ne.0) return
            pelement%pntin%value=val
            pelement%pntin%row_ind=i
            pelement%pntin%col_ind=j
         else 
            call access_element(pelement,ind,pmatrix,istat)
            if (istat.ne.0) return
            pelement%pntin%value= pelement%pntin%value+val
            call  dealloc_c_element (nmb_element,pmatrix,istat)
            if (istat.ne.0) return
         end if
      end if
      end subroutine cINS_entry  
!*
      subroutine cINS_block (pmatrix,val,i,j,istat)
      implicit none
      type( c_matrix ),pointer ::pmatrix
      complex(KIND=sp)  ,dimension(:,:) ,target,intent(in)::val
      integer ,intent(in)::i,j
      integer,intent(out)::istat
      integer ::nmb_element,ind,ierr
      complex(KIND=sp)   ,dimension(:,:),allocatable,target::vv
      type(c_inelement ),pointer::pelement
      integer ::s_rows,s_cols
      istat = -1
      s_rows=size(val,1)
      s_cols=size(val,2)
      allocate(vv(s_rows,s_cols),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif
      if((i.gt.pmatrix%DIM(3).or.(j.gt.pmatrix%DIM(4))&
         .or.(s_rows.ne.pmatrix%DIM(5)&
         .or.(s_cols.ne.pmatrix%DIM(6))))) then 
         istat = blas_error_param
         return 
      else
         call new_c_element (pmatrix,nmb_element,istat) 
         if (istat.ne.0) return
         call c_element_num  (ind,pmatrix,i,j,istat)
         if (istat.ne.0) return
         if(ind.eq.0) then
            call access_element(pelement,nmb_element,&
                                pmatrix,istat)
            if (istat.ne.0) return
            allocate(pelement%blin%value(s_rows,s_cols),&
                     STAT=ierr)
            if(ierr.ne.0) then
               istat = blas_error_memalloc
               return
            endif
            pelement%blin%value=val
            pelement%blin%row_block_ind=i
            pelement%blin%col_block_ind=j    
         else
            call access_element(pelement,ind,pmatrix,istat)     
            if (istat.ne.0) return
            vv=vv+val 
            pelement%blin%value=pelement%blin%value+val
            call dealloc_c_element (nmb_element,pmatrix,istat)
            if (istat.ne.0) return
         end if
      end if
      deallocate(vv,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif
      end subroutine cINS_block 
!*
      subroutine cINS_bl_entr (pmatrix,val,i,j,istat)
      implicit none
      type(c_matrix ),pointer ::pmatrix
      complex(KIND=sp)   ,intent(in)::val
      integer,intent(in)::i,j
      integer,intent(out)::istat
      complex(KIND=sp)  ,dimension(:,:),allocatable,target ::vall
      integer::ii,jj,ierr
      istat = -1
      ii=floor(real((i-1)/(pmatrix%DIM(5))))
      jj=floor(real((j-1)/(pmatrix%DIM(6))))
      allocate(vall(pmatrix%DIM(5),pmatrix%DIM(6)),&
               STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      vall= (0.0e0, 0.0e0) 
      vall(i-ii*pmatrix%DIM(5),j-jj*pmatrix%DIM(6))=val
      call cINS_block (pmatrix,vall,ii+1,jj+1,istat)
      if (istat.ne.0) return     
      deallocate(vall,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif
      end subroutine cINS_bl_entr 
!*
      subroutine cINS_varblock (vpmatrix,val,i,j,istat)
      implicit none
      type(c_matrix ),pointer ::vpmatrix
      complex(KIND=sp)  ,dimension(:,:) ,target,intent(in)::val
      integer ,intent(in)::i,j
      integer,intent(out)::istat
      integer ::nmb_element,ind
      integer::ierr
      complex(KIND=sp)  ,dimension(:,:),allocatable,target::vv
      type(c_inelement  ),pointer::pelement
      integer ::s_rows,s_cols,k
      istat = -1
      s_rows=size(val,1)
      s_cols=size(val,2)
      allocate(vv(s_rows,s_cols),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      if((i.gt.vpmatrix%DIM(3).or.j.gt.vpmatrix%DIM(4))&
        .or.(s_rows.ne.vpmatrix%sub_rows(i))&
        .or.(s_cols.ne.vpmatrix%sub_cols(j))) then 
         istat = blas_error_param
         return 
      else
         call new_c_element  (vpmatrix,nmb_element,istat) 
         if (istat.ne.0) return     
         call c_element_num  (ind,vpmatrix,i,j,istat)
         if (istat.ne.0) return     
         if(ind.eq.0) then
            call access_element(pelement,nmb_element,&
                                vpmatrix,istat)
            if (istat.ne.0) return     
            allocate(pelement%vblin%value(vpmatrix%sub_rows(i),&
                                vpmatrix%sub_cols(j)),STAT=ierr)
            if(ierr.ne.0) then
               istat = blas_error_memalloc
               return
            endif   
            pelement%vblin%value=val
            pelement%vblin%row_vblock_ind=i
            pelement%vblin%col_vblock_ind=j
            do k=i+1,vpmatrix%DIM(3)
               vpmatrix%trb(k)=vpmatrix%trb(k)+1
            end do
            do k=1,vpmatrix%DIM(3)-1
               vpmatrix%tre(k)=vpmatrix%trb(k+1)
            end do
            vpmatrix%tre(vpmatrix%DIM(3))=nmb_element+1
         else
            call access_element(pelement,ind,vpmatrix,istat)
            if (istat.ne.0) return     
            pelement%vblin%value= pelement%vblin%value+val
            call dealloc_c_element (nmb_element,vpmatrix,istat)
            if (istat.ne.0) return     
         end if
      end if
      istat = 0
      return
      end subroutine cINS_varblock 
!*
      subroutine cINS_varbl_entr (vpmatrix,val,i,j,istat)
      implicit none
      type(c_matrix ),pointer ::vpmatrix
      complex(KIND=sp) ,intent(in)::val
      integer,intent(in)::i,j
      integer,intent(out)::istat
      complex(KIND=sp) ,dimension(:,:),allocatable ::vall
      integer::ii,jj,k,p,ind1,ind2,vall_ind1,vall_ind2,ierr
      ! determine the row of block entring
      ind1=0
      do k=1,vpmatrix%DIM(3)
         ind1=ind1+vpmatrix%sub_rows(k)
         if(ind1.ge.i) exit
      end do
      if(k.le.vpmatrix%DIM(3)) then
         ii=k
         vall_ind1=i-(ind1-vpmatrix%sub_rows(k))
      else
         istat = blas_error_param
         return
      end if
      ! determine the col of block entring
      ind2=0
      do p=1,vpmatrix%DIM(3)
         ind2=ind2+vpmatrix%sub_cols(p)
         if(ind2.ge.j) exit
      end do
      if(p.le.vpmatrix%DIM(4)) then
         jj=p
         vall_ind2=j-(ind2-vpmatrix%sub_cols(p))
      else
         istat = blas_error_param
         return
      end if
      allocate(vall(vpmatrix%sub_rows(ii),&
                    vpmatrix%sub_cols(jj)),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      vall= (0.0e0, 0.0e0) 
      vall(vall_ind1,vall_ind2)=val    
      call  cINS_varblock (vpmatrix,vall,ii,jj,istat)
      if (istat.ne.0) return     
      deallocate(vall,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end subroutine cINS_varbl_entr 
!*
      subroutine cuscr_varend (a,prpty,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer,intent(inout)::prpty
      integer:: m,n,ierr,ind,kb,mb
      integer, dimension(:),allocatable :: bindx,indx,rpntr,&
                                        cpntr,bpntrb,bpntre
      complex(KIND=sp) , dimension(:),allocatable :: val
      integer :: size_val,val_ind,indx_ind,bindx_ind,ii,jj,i,j
      type(c_matrix ),pointer::pmatrix
      type(c_inelement ),pointer::pelement
      istat=-1
      call access_matrix(pmatrix,a,istat)
      if (istat.ne.0) return     
      mb=pmatrix%DIM(3)
      kb=pmatrix%DIM(4)
      ! determine  size of val,m,n 
      size_val=0
      m=0
      n=0
      do i=1,pmatrix%DIM(3)
         do j=1,pmatrix%DIM(4)
            call c_element_num (ind,pmatrix,i,j,ierr)
            if(ind.ne.0) then
               call access_element(pelement,ind,pmatrix,istat)
               if (istat.ne.0) return     
               size_val=size_val+&
                      pmatrix%sub_rows(i)*pmatrix%sub_cols(j)
            end if
         end do
      end do
      do i=1,pmatrix%DIM(3)
         m=m+pmatrix%sub_rows(i)
      end do
      do j=1,pmatrix%DIM(4)     
         n=n+pmatrix%sub_cols(j)
      end do
      allocate(val(size_val),&
              indx(pmatrix% c_element_start %number+1),&
              bindx(pmatrix% c_element_start %number),&
              STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      val= (0.0e0, 0.0e0) 
      ! fill val ,indx and bindx
      val_ind=0
      indx_ind=1
      bindx_ind=0
      indx(1)=1
      do i=1,pmatrix%DIM(3)
         do j=1,pmatrix%DIM(4) 
            call c_element_num (ind,pmatrix,i,j,istat)
            if (istat.ne.0) return     
            if(ind.ne.0) then
               call access_element(pelement,ind,pmatrix,istat)
               if (istat.ne.0) return     
               do jj=1,pmatrix%sub_cols(j)
                  do ii=1,pmatrix%sub_rows(i)
                     val_ind=val_ind+1
                     val(val_ind)=pelement%vblin%value(ii,jj)  
                  end do
               end do
               bindx_ind=bindx_ind+1
               bindx(bindx_ind)=j
               indx_ind=indx_ind+1
               indx(indx_ind)=indx(indx_ind-1)&
                      +pmatrix%sub_rows(i)*pmatrix%sub_cols(j)
            end if
         end do
      end do
      ! fill rpntr, cpntr,bpntrb,bpntre
      allocate(rpntr(pmatrix%DIM(3)+1),&
              cpntr(pmatrix%DIM(4)+1),&
              bpntrb(pmatrix%DIM(3)),&
              bpntre(pmatrix%DIM(3)),&
              STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      rpntr(1)=1
      do i=2,pmatrix%DIM(3)+1
         rpntr(i)= rpntr(i-1)+pmatrix%sub_rows(i-1)
      end do
      cpntr(1)=1
      do j=2,pmatrix%DIM(4)+1
         cpntr(j)= cpntr(j-1)+pmatrix%sub_cols(j-1)
      end do
      do i=1,pmatrix%DIM(3)
         bpntrb(i)=pmatrix%trb(i)
         bpntre(i)=pmatrix%tre(i)
      end do
      ! RELEASING 
      call  c_dealloc (a,istat)
      if (istat.ne.0) return     
      ! CREATING MATRIX IN VBR FORMAT
      istat=-1 !needed to create copy of data
      call cuscr_vbr (m,n,val,indx,bindx,rpntr,cpntr,&
                    bpntrb,bpntre,mb,kb,prpty,istat,a)
      if (istat.ne.0) return     
      deallocate(val,bindx,indx,rpntr, &
                cpntr,bpntrb,bpntre,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end subroutine  cuscr_varend  
!*
      subroutine cuscr_normend  (a,prpty,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer,intent(inout)::prpty
      integer:: m,n,nnz
      integer, dimension(:),allocatable :: indx,jndx
      complex(KIND=sp) , dimension(:),allocatable :: val 
      integer :: nmb_element,i,ierr
      type(c_matrix  ),pointer::pmatrix
      type(c_inelement ),pointer::pelement
      call access_matrix(pmatrix,a,istat)
      if (istat.ne.0) return     
      m=pmatrix%DIM(1)  !nb_of_rows
      n=pmatrix%DIM(2)  !nb_of_cols
      nnz=pmatrix% c_element_start %number
      allocate(val(nnz),indx(nnz),jndx(nnz),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      nmb_element=nnz+1
      do i=1,nnz
         call access_element(pelement,nmb_element-i,&
                             pmatrix,istat)
         if (istat.ne.0) return     
         val(i)=pelement%pntin%value
         indx(i)=pelement%pntin%row_ind
         jndx(i)=pelement%pntin%col_ind
      end do
      call  c_dealloc (a,istat)
      if (istat.ne.0) return     
      ! CREATING A MATRIX IN COO FORMAT
      istat=-1 !needed to create copy of data
      call cuscr_coo (m,n,val,indx,jndx,nnz,prpty,istat,a) 
      if (istat.ne.0) return     
      deallocate(val,indx,jndx,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end  subroutine cuscr_normend  
!*
      subroutine cuscr_blockend  (a,prpty,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer,intent(inout)::prpty
      integer:: m,n,bnnz,ierr,ind,kb,lb,dummy,k,mb,i,j
      integer, dimension(:),allocatable :: bindx,bjndx
      complex(KIND=sp)  , dimension(:),allocatable :: val
      integer :: nmb_block
      type(c_matrix ),pointer::pmatrix
      type(c_inelement ),pointer::pelement
      istat=-1
      call access_matrix(pmatrix,a,istat)
      if (istat.ne.0) return     
      lb=pmatrix%DIM(5)
      bnnz=pmatrix% c_element_start %number
      mb=pmatrix%DIM(3)
      kb=pmatrix%DIM(4)
      m=mb*lb
      n=kb*lb
      ind=0
      dummy=bnnz*lb*lb
      allocate(val(dummy),bindx(bnnz),bjndx(bnnz),&
               STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      nmb_block=bnnz+1
      do i=1,bnnz
         k=1
         call access_element(pelement,nmb_block-i,pmatrix,&
                             istat)
         if (istat.ne.0) return     
         do j=1,lb
            do k=1,lb
               ind=ind+1
               val(ind)=pelement%blin%value(k,j)          
            end do
         end do
         bindx(i)=pelement%blin%row_block_ind
         bjndx(i)=pelement%blin%col_block_ind
      end do
      ! RELEASING
      call c_dealloc (a,istat)
      if (istat.ne.0) return     
      ! CREATE A MATRIX IN BCO FORMAT
      istat=-1 !needed to create copy of data
      call  cuscr_bco (m,n,val,bindx,bjndx,bnnz,&
                      mb,kb,lb,prpty,istat,a)
      if (istat.ne.0) return     
      deallocate(val,bindx,bjndx,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end subroutine  cuscr_blockend  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      subroutine zINS_entry  (pmatrix,val,i,j,istat)
      implicit none
      type(z_matrix ),pointer ::pmatrix
      complex(KIND=dp)  ,intent(in) ::val
      integer ,intent(in) ::i,j
      integer, intent(out) :: istat
      type(z_inelement ),pointer ::pelement
      integer::nmb_element,ind
      istat=-1
      if((i.gt.pmatrix%DIM(1)).or.&
           (j.gt.pmatrix%DIM(2))) then
         istat = blas_error_param 
         return 
      else
         call new_z_element (pmatrix,nmb_element,istat)
         if (istat.ne.0) return
         call z_element_num (ind,pmatrix,i,j,istat)
         if (istat.ne.0) return
         if(ind.eq.0) then
            call access_element(pelement,nmb_element,&
                                pmatrix,istat)
            if (istat.ne.0) return
            pelement%pntin%value=val
            pelement%pntin%row_ind=i
            pelement%pntin%col_ind=j
         else 
            call access_element(pelement,ind,pmatrix,istat)
            if (istat.ne.0) return
            pelement%pntin%value= pelement%pntin%value+val
            call  dealloc_z_element (nmb_element,pmatrix,istat)
            if (istat.ne.0) return
         end if
      end if
      end subroutine zINS_entry  
!*
      subroutine zINS_block (pmatrix,val,i,j,istat)
      implicit none
      type( z_matrix ),pointer ::pmatrix
      complex(KIND=dp)  ,dimension(:,:) ,target,intent(in)::val
      integer ,intent(in)::i,j
      integer,intent(out)::istat
      integer ::nmb_element,ind,ierr
      complex(KIND=dp)   ,dimension(:,:),allocatable,target::vv
      type(z_inelement ),pointer::pelement
      integer ::s_rows,s_cols
      istat = -1
      s_rows=size(val,1)
      s_cols=size(val,2)
      allocate(vv(s_rows,s_cols),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif
      if((i.gt.pmatrix%DIM(3).or.(j.gt.pmatrix%DIM(4))&
         .or.(s_rows.ne.pmatrix%DIM(5)&
         .or.(s_cols.ne.pmatrix%DIM(6))))) then 
         istat = blas_error_param
         return 
      else
         call new_z_element (pmatrix,nmb_element,istat) 
         if (istat.ne.0) return
         call z_element_num  (ind,pmatrix,i,j,istat)
         if (istat.ne.0) return
         if(ind.eq.0) then
            call access_element(pelement,nmb_element,&
                                pmatrix,istat)
            if (istat.ne.0) return
            allocate(pelement%blin%value(s_rows,s_cols),&
                     STAT=ierr)
            if(ierr.ne.0) then
               istat = blas_error_memalloc
               return
            endif
            pelement%blin%value=val
            pelement%blin%row_block_ind=i
            pelement%blin%col_block_ind=j    
         else
            call access_element(pelement,ind,pmatrix,istat)     
            if (istat.ne.0) return
            vv=vv+val 
            pelement%blin%value=pelement%blin%value+val
            call dealloc_z_element (nmb_element,pmatrix,istat)
            if (istat.ne.0) return
         end if
      end if
      deallocate(vv,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif
      end subroutine zINS_block 
!*
      subroutine zINS_bl_entr (pmatrix,val,i,j,istat)
      implicit none
      type(z_matrix ),pointer ::pmatrix
      complex(KIND=dp)   ,intent(in)::val
      integer,intent(in)::i,j
      integer,intent(out)::istat
      complex(KIND=dp)  ,dimension(:,:),allocatable,target ::vall
      integer::ii,jj,ierr
      istat = -1
      ii=floor(real((i-1)/(pmatrix%DIM(5))))
      jj=floor(real((j-1)/(pmatrix%DIM(6))))
      allocate(vall(pmatrix%DIM(5),pmatrix%DIM(6)),&
               STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      vall= (0.0d0, 0.0d0) 
      vall(i-ii*pmatrix%DIM(5),j-jj*pmatrix%DIM(6))=val
      call zINS_block (pmatrix,vall,ii+1,jj+1,istat)
      if (istat.ne.0) return     
      deallocate(vall,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif
      end subroutine zINS_bl_entr 
!*
      subroutine zINS_varblock (vpmatrix,val,i,j,istat)
      implicit none
      type(z_matrix ),pointer ::vpmatrix
      complex(KIND=dp)  ,dimension(:,:) ,target,intent(in)::val
      integer ,intent(in)::i,j
      integer,intent(out)::istat
      integer ::nmb_element,ind
      integer::ierr
      complex(KIND=dp)  ,dimension(:,:),allocatable,target::vv
      type(z_inelement  ),pointer::pelement
      integer ::s_rows,s_cols,k
      istat = -1
      s_rows=size(val,1)
      s_cols=size(val,2)
      allocate(vv(s_rows,s_cols),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      if((i.gt.vpmatrix%DIM(3).or.j.gt.vpmatrix%DIM(4))&
        .or.(s_rows.ne.vpmatrix%sub_rows(i))&
        .or.(s_cols.ne.vpmatrix%sub_cols(j))) then 
         istat = blas_error_param
         return 
      else
         call new_z_element  (vpmatrix,nmb_element,istat) 
         if (istat.ne.0) return     
         call z_element_num  (ind,vpmatrix,i,j,istat)
         if (istat.ne.0) return     
         if(ind.eq.0) then
            call access_element(pelement,nmb_element,&
                                vpmatrix,istat)
            if (istat.ne.0) return     
            allocate(pelement%vblin%value(vpmatrix%sub_rows(i),&
                                vpmatrix%sub_cols(j)),STAT=ierr)
            if(ierr.ne.0) then
               istat = blas_error_memalloc
               return
            endif   
            pelement%vblin%value=val
            pelement%vblin%row_vblock_ind=i
            pelement%vblin%col_vblock_ind=j
            do k=i+1,vpmatrix%DIM(3)
               vpmatrix%trb(k)=vpmatrix%trb(k)+1
            end do
            do k=1,vpmatrix%DIM(3)-1
               vpmatrix%tre(k)=vpmatrix%trb(k+1)
            end do
            vpmatrix%tre(vpmatrix%DIM(3))=nmb_element+1
         else
            call access_element(pelement,ind,vpmatrix,istat)
            if (istat.ne.0) return     
            pelement%vblin%value= pelement%vblin%value+val
            call dealloc_z_element (nmb_element,vpmatrix,istat)
            if (istat.ne.0) return     
         end if
      end if
      istat = 0
      return
      end subroutine zINS_varblock 
!*
      subroutine zINS_varbl_entr (vpmatrix,val,i,j,istat)
      implicit none
      type(z_matrix ),pointer ::vpmatrix
      complex(KIND=dp) ,intent(in)::val
      integer,intent(in)::i,j
      integer,intent(out)::istat
      complex(KIND=dp) ,dimension(:,:),allocatable ::vall
      integer::ii,jj,k,p,ind1,ind2,vall_ind1,vall_ind2,ierr
      ! determine the row of block entring
      ind1=0
      do k=1,vpmatrix%DIM(3)
         ind1=ind1+vpmatrix%sub_rows(k)
         if(ind1.ge.i) exit
      end do
      if(k.le.vpmatrix%DIM(3)) then
         ii=k
         vall_ind1=i-(ind1-vpmatrix%sub_rows(k))
      else
         istat = blas_error_param
         return
      end if
      ! determine the col of block entring
      ind2=0
      do p=1,vpmatrix%DIM(3)
         ind2=ind2+vpmatrix%sub_cols(p)
         if(ind2.ge.j) exit
      end do
      if(p.le.vpmatrix%DIM(4)) then
         jj=p
         vall_ind2=j-(ind2-vpmatrix%sub_cols(p))
      else
         istat = blas_error_param
         return
      end if
      allocate(vall(vpmatrix%sub_rows(ii),&
                    vpmatrix%sub_cols(jj)),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      vall= (0.0d0, 0.0d0) 
      vall(vall_ind1,vall_ind2)=val    
      call  zINS_varblock (vpmatrix,vall,ii,jj,istat)
      if (istat.ne.0) return     
      deallocate(vall,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end subroutine zINS_varbl_entr 
!*
      subroutine zuscr_varend (a,prpty,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer,intent(inout)::prpty
      integer:: m,n,ierr,ind,kb,mb
      integer, dimension(:),allocatable :: bindx,indx,rpntr,&
                                        cpntr,bpntrb,bpntre
      complex(KIND=dp) , dimension(:),allocatable :: val
      integer :: size_val,val_ind,indx_ind,bindx_ind,ii,jj,i,j
      type(z_matrix ),pointer::pmatrix
      type(z_inelement ),pointer::pelement
      istat=-1
      call access_matrix(pmatrix,a,istat)
      if (istat.ne.0) return     
      mb=pmatrix%DIM(3)
      kb=pmatrix%DIM(4)
      ! determine  size of val,m,n 
      size_val=0
      m=0
      n=0
      do i=1,pmatrix%DIM(3)
         do j=1,pmatrix%DIM(4)
            call z_element_num (ind,pmatrix,i,j,ierr)
            if(ind.ne.0) then
               call access_element(pelement,ind,pmatrix,istat)
               if (istat.ne.0) return     
               size_val=size_val+&
                      pmatrix%sub_rows(i)*pmatrix%sub_cols(j)
            end if
         end do
      end do
      do i=1,pmatrix%DIM(3)
         m=m+pmatrix%sub_rows(i)
      end do
      do j=1,pmatrix%DIM(4)     
         n=n+pmatrix%sub_cols(j)
      end do
      allocate(val(size_val),&
              indx(pmatrix% z_element_start %number+1),&
              bindx(pmatrix% z_element_start %number),&
              STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      val= (0.0d0, 0.0d0) 
      ! fill val ,indx and bindx
      val_ind=0
      indx_ind=1
      bindx_ind=0
      indx(1)=1
      do i=1,pmatrix%DIM(3)
         do j=1,pmatrix%DIM(4) 
            call z_element_num (ind,pmatrix,i,j,istat)
            if (istat.ne.0) return     
            if(ind.ne.0) then
               call access_element(pelement,ind,pmatrix,istat)
               if (istat.ne.0) return     
               do jj=1,pmatrix%sub_cols(j)
                  do ii=1,pmatrix%sub_rows(i)
                     val_ind=val_ind+1
                     val(val_ind)=pelement%vblin%value(ii,jj)  
                  end do
               end do
               bindx_ind=bindx_ind+1
               bindx(bindx_ind)=j
               indx_ind=indx_ind+1
               indx(indx_ind)=indx(indx_ind-1)&
                      +pmatrix%sub_rows(i)*pmatrix%sub_cols(j)
            end if
         end do
      end do
      ! fill rpntr, cpntr,bpntrb,bpntre
      allocate(rpntr(pmatrix%DIM(3)+1),&
              cpntr(pmatrix%DIM(4)+1),&
              bpntrb(pmatrix%DIM(3)),&
              bpntre(pmatrix%DIM(3)),&
              STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      rpntr(1)=1
      do i=2,pmatrix%DIM(3)+1
         rpntr(i)= rpntr(i-1)+pmatrix%sub_rows(i-1)
      end do
      cpntr(1)=1
      do j=2,pmatrix%DIM(4)+1
         cpntr(j)= cpntr(j-1)+pmatrix%sub_cols(j-1)
      end do
      do i=1,pmatrix%DIM(3)
         bpntrb(i)=pmatrix%trb(i)
         bpntre(i)=pmatrix%tre(i)
      end do
      ! RELEASING 
      call  z_dealloc (a,istat)
      if (istat.ne.0) return     
      ! CREATING MATRIX IN VBR FORMAT
      istat=-1 !needed to create copy of data
      call zuscr_vbr (m,n,val,indx,bindx,rpntr,cpntr,&
                    bpntrb,bpntre,mb,kb,prpty,istat,a)
      if (istat.ne.0) return     
      deallocate(val,bindx,indx,rpntr, &
                cpntr,bpntrb,bpntre,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end subroutine  zuscr_varend  
!*
      subroutine zuscr_normend  (a,prpty,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer,intent(inout)::prpty
      integer:: m,n,nnz
      integer, dimension(:),allocatable :: indx,jndx
      complex(KIND=dp) , dimension(:),allocatable :: val 
      integer :: nmb_element,i,ierr
      type(z_matrix  ),pointer::pmatrix
      type(z_inelement ),pointer::pelement
      call access_matrix(pmatrix,a,istat)
      if (istat.ne.0) return     
      m=pmatrix%DIM(1)  !nb_of_rows
      n=pmatrix%DIM(2)  !nb_of_cols
      nnz=pmatrix% z_element_start %number
      allocate(val(nnz),indx(nnz),jndx(nnz),STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      nmb_element=nnz+1
      do i=1,nnz
         call access_element(pelement,nmb_element-i,&
                             pmatrix,istat)
         if (istat.ne.0) return     
         val(i)=pelement%pntin%value
         indx(i)=pelement%pntin%row_ind
         jndx(i)=pelement%pntin%col_ind
      end do
      call  z_dealloc (a,istat)
      if (istat.ne.0) return     
      ! CREATING A MATRIX IN COO FORMAT
      istat=-1 !needed to create copy of data
      call zuscr_coo (m,n,val,indx,jndx,nnz,prpty,istat,a) 
      if (istat.ne.0) return     
      deallocate(val,indx,jndx,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end  subroutine zuscr_normend  
!*
      subroutine zuscr_blockend  (a,prpty,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer,intent(inout)::prpty
      integer:: m,n,bnnz,ierr,ind,kb,lb,dummy,k,mb,i,j
      integer, dimension(:),allocatable :: bindx,bjndx
      complex(KIND=dp)  , dimension(:),allocatable :: val
      integer :: nmb_block
      type(z_matrix ),pointer::pmatrix
      type(z_inelement ),pointer::pelement
      istat=-1
      call access_matrix(pmatrix,a,istat)
      if (istat.ne.0) return     
      lb=pmatrix%DIM(5)
      bnnz=pmatrix% z_element_start %number
      mb=pmatrix%DIM(3)
      kb=pmatrix%DIM(4)
      m=mb*lb
      n=kb*lb
      ind=0
      dummy=bnnz*lb*lb
      allocate(val(dummy),bindx(bnnz),bjndx(bnnz),&
               STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      nmb_block=bnnz+1
      do i=1,bnnz
         k=1
         call access_element(pelement,nmb_block-i,pmatrix,&
                             istat)
         if (istat.ne.0) return     
         do j=1,lb
            do k=1,lb
               ind=ind+1
               val(ind)=pelement%blin%value(k,j)          
            end do
         end do
         bindx(i)=pelement%blin%row_block_ind
         bjndx(i)=pelement%blin%col_block_ind
      end do
      ! RELEASING
      call z_dealloc (a,istat)
      if (istat.ne.0) return     
      ! CREATE A MATRIX IN BCO FORMAT
      istat=-1 !needed to create copy of data
      call  zuscr_bco (m,n,val,bindx,bjndx,bnnz,&
                      mb,kb,lb,prpty,istat,a)
      if (istat.ne.0) return     
      deallocate(val,bindx,bjndx,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end subroutine  zuscr_blockend  
! **********************************************************************
!!*************************************************************************** 
! **********************************************************************
      end module mod_INS_ROUTINER

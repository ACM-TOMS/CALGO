      subroutine XINS_entry (pmatrix,val,i,j,istat)
      implicit none
      
      type(X_MATRIX),pointer ::pmatrix
      DCOMPLEX ,intent(in) ::val
      integer ,intent(in) ::i,j
      integer, intent(out) :: istat
      type(X_INELEMENT),pointer ::pelement
      integer::nmb_element,ind

      istat=-1
      if((i.gt.pmatrix%DIM(1)).or.&
           (j.gt.pmatrix%DIM(2))) then
         istat = blas_error_param 
         return 
      else
         call NEW_X_ELEMENT(pmatrix,nmb_element,istat)
         if (istat.ne.0) return
         call X_ELEMENT_NUM(ind,pmatrix,i,j,istat)
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
            call  DELOC_X_ELEMENT(nmb_element,pmatrix,istat)
            if (istat.ne.0) return
         end if
      end if
      
      end subroutine XINS_entry 
!!!
!*
!!!
      subroutine XINS_block(pmatrix,val,i,j,istat)
      implicit none
      type( X_MATRIX),pointer ::pmatrix
      DCOMPLEX ,dimension(:,:) ,target,intent(in)::val
      integer ,intent(in)::i,j
      integer,intent(out)::istat
      integer ::nmb_element,ind,ierr
      DCOMPLEX  ,dimension(:,:),allocatable,target::vv
      type(X_INELEMENT),pointer::pelement
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
         call NEW_X_ELEMENT(pmatrix,nmb_element,istat) 
         if (istat.ne.0) return
         call X_ELEMENT_NUM (ind,pmatrix,i,j,istat)
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
            call DELOC_X_ELEMENT(nmb_element,pmatrix,istat)
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
      end subroutine XINS_block
!!!
!*
!!!
      subroutine XINS_bl_entr(pmatrix,val,i,j,istat)
      implicit none
      type(X_MATRIX),pointer ::pmatrix
      DCOMPLEX  ,intent(in)::val
      integer,intent(in)::i,j
      integer,intent(out)::istat
      DCOMPLEX ,dimension(:,:),allocatable,target ::vall
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
      vall=ZZERO
      vall(i-ii*pmatrix%DIM(5),j-jj*pmatrix%DIM(6))=val
      call XINS_block(pmatrix,vall,ii+1,jj+1,istat)
      if (istat.ne.0) return     
      deallocate(vall,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif
      end subroutine XINS_bl_entr
!!!
!*
!!!
      subroutine XINS_varblock(vpmatrix,val,i,j,istat)
      implicit none
      type(X_MATRIX),pointer ::vpmatrix
      DCOMPLEX ,dimension(:,:) ,target,intent(in)::val
      integer ,intent(in)::i,j
      integer,intent(out)::istat
      integer ::nmb_element,ind
      integer::ierr
      DCOMPLEX ,dimension(:,:),allocatable,target::vv
      type(X_INELEMENT ),pointer::pelement
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
         call NEW_X_ELEMENT (vpmatrix,nmb_element,istat) 
         if (istat.ne.0) return     
         call X_ELEMENT_NUM (ind,vpmatrix,i,j,istat)
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
            call DELOC_X_ELEMENT(nmb_element,vpmatrix,istat)
            if (istat.ne.0) return     
         end if
      end if
      istat = 0
      return
      end subroutine XINS_varblock
!!!
!*
!!!
      subroutine XINS_varbl_entr(vpmatrix,val,i,j,istat)
      implicit none
      type(X_MATRIX),pointer ::vpmatrix
      DCOMPLEX,intent(in)::val
      integer,intent(in)::i,j
      integer,intent(out)::istat
      DCOMPLEX,dimension(:,:),allocatable ::vall
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

      vall=ZZERO
      vall(vall_ind1,vall_ind2)=val    
      call  XINS_varblock(vpmatrix,vall,ii,jj,istat)
      if (istat.ne.0) return     
      deallocate(vall,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   

      end subroutine XINS_varbl_entr
!!!
!*
!!!
      subroutine XUSCR_varend(a,prpty,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer,intent(inout)::prpty
      integer:: m,n,ierr,ind,kb,mb
      integer, dimension(:),allocatable :: bindx,indx,rpntr,&
                                        cpntr,bpntrb,bpntre
      DCOMPLEX, dimension(:),allocatable :: val
      integer :: size_val,val_ind,indx_ind,bindx_ind,ii,jj,i,j
      type(X_MATRIX),pointer::pmatrix
      type(X_INELEMENT),pointer::pelement
 
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
            call X_ELEMENT_NUM(ind,pmatrix,i,j,ierr)
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
              indx(pmatrix%X_ELEMENT_start%number+1),&
              bindx(pmatrix%X_ELEMENT_start%number),&
              STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memalloc
         return
      endif   
      val=ZZERO

      ! fill val ,indx and bindx
      val_ind=0
      indx_ind=1
      bindx_ind=0
      indx(1)=1
      do i=1,pmatrix%DIM(3)
         do j=1,pmatrix%DIM(4) 
            call X_ELEMENT_NUM(ind,pmatrix,i,j,istat)
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
      call  X_DEALLOC(a,istat)
      if (istat.ne.0) return     
     
      ! CREATING MATRIX IN VBR FORMAT
      istat=-1 !needed to create copy of data
      call XUSCR_VBR(m,n,val,indx,bindx,rpntr,cpntr,&
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
      end subroutine  XUSCR_varend 
!!!
!*
!!!
      subroutine XUSCR_normend (a,prpty,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer,intent(inout)::prpty
      integer:: m,n,nnz
      integer, dimension(:),allocatable :: indx,jndx
      DCOMPLEX, dimension(:),allocatable :: val 
      integer :: nmb_element,i,ierr
      type(X_MATRIX ),pointer::pmatrix
      type(X_INELEMENT),pointer::pelement
     
      call access_matrix(pmatrix,a,istat)
      if (istat.ne.0) return     

      m=pmatrix%DIM(1)  !nb_of_rows
      n=pmatrix%DIM(2)  !nb_of_cols
      nnz=pmatrix%X_ELEMENT_start%number
     
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
 
      call  X_DEALLOC(a,istat)
      if (istat.ne.0) return     
           
      ! CREATING A MATRIX IN COO FORMAT
      istat=-1 !needed to create copy of data
      call XUSCR_COO(m,n,val,indx,jndx,nnz,prpty,istat,a) 
      if (istat.ne.0) return     

      deallocate(val,indx,jndx,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   
      end  subroutine XUSCR_normend 
!!!
!*
!!!
      subroutine XUSCR_blockend (a,prpty,istat)
      implicit none
      integer ,intent(inout)::a,istat
      integer,intent(inout)::prpty
      integer:: m,n,bnnz,ierr,ind,kb,lb,dummy,k,mb,i,j
      integer, dimension(:),allocatable :: bindx,bjndx
      DCOMPLEX , dimension(:),allocatable :: val
      integer :: nmb_block
      type(X_MATRIX),pointer::pmatrix
      type(X_INELEMENT),pointer::pelement
      

      istat=-1
      call access_matrix(pmatrix,a,istat)
      if (istat.ne.0) return     

      lb=pmatrix%DIM(5)
      bnnz=pmatrix%X_ELEMENT_start%number
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
      call X_DEALLOC(a,istat)
      if (istat.ne.0) return     

      ! CREATE A MATRIX IN BCO FORMAT
      istat=-1 !needed to create copy of data
      call  XUSCR_BCO(m,n,val,bindx,bjndx,bnnz,&
                      mb,kb,lb,prpty,istat,a)
      if (istat.ne.0) return     

      deallocate(val,bindx,bjndx,STAT=ierr)
      if(ierr.ne.0) then
         istat = blas_error_memdeloc
         return
      else
         istat = 0
      endif   

      end subroutine  XUSCR_blockend 

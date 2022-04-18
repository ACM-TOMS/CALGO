      subroutine XB_UP_ORDER(VAL,lbxlb,BINDX)
      implicit none
      DCOMPLEX,pointer  ,dimension(:)::VAL
      integer ,pointer, dimension(:)::BINDX
      integer ,dimension(:),allocatable :: tes
      DCOMPLEX,dimension(:,:),allocatable::P
      integer ,intent(in) ::lbxlb
      integer ::s,i,j,k
   
      allocate(tes(lbxlb))
      do i=1,lbxlb
         tes(i)=i
      end do
      s=size(VAL)
      k=floor(real(s/lbxlb))
      allocate(P(lbxlb,k))
      do j=1,s
         i=floor(real((j-1)/lbxlb))
         P(j-i*lbxlb,i+1)=VAL(j)
      end do
      P=P(tes,BINDX)
      do j=1,s
         i=floor(real((j-1)/lbxlb))
         VAL(j)=P(j-i*lbxlb,i+1)
      end do
      deallocate(tes,P)
      end subroutine XB_UP_ORDER 

      function XA_ROW_COL(VAL,INDX,JNDX,i,j)
      DCOMPLEX,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX,JNDX
      integer,intent(in)::i,j
      integer::k
      logical::finder
      DCOMPLEX::XA_ROW_COL 
      
      finder=.false.
      k=1
      do while((k.le.size(VAL)).and.(.not.finder))
         if(INDX(k).eq.i.and.JNDX(k).eq.j) then
            finder=.true.
         else
            k=k+1
         end if
      end do
   
      if(finder) then
         XA_ROW_COL =VAL(k)
      else
         XA_ROW_COL =0
      end if
      end function XA_ROW_COL

      subroutine XDETECT_DIAG(VAL,INDX,JNDX,ind,LDA,test)
      implicit none
      DCOMPLEX,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX,JNDX
      integer,intent(in)::ind,LDA
      integer,intent(inout)::test
      logical::finder
      DCOMPLEX::val_val
      integer ::k

      test=0
      if(ind.eq.0) then ! main diag
         test=0       
         finder=.false.
         k=1
         do while((k.le.LDA).and.(.not.finder))
            val_val=XA_ROW_COL(VAL,INDX,JNDX,k,k)
            if(val_val.ne.0) then
               finder=.true.
            else
               k=k+1
            end if   
         end do
         if(finder) then
            test=1
         end if
      elseif(ind.lt.0)then ! low diag
         test=0  
         finder=.false.
         k=-ind+1
         do while((k.le.LDA).and.(.not.finder))
            val_val=A_row_col(VAL,INDX,JNDX,k,k+ind)
            if(val_val.ne.0) then
               finder=.true.
            else
               k=k+1
            end if
         end do  
         if(finder) then
            test=1
         end if
      else ! high diag
         test=0
         finder=.false.
         k=ind+1
         do while((k.le.LDA).and.(.not.finder))
            val_val=A_row_col(VAL,INDX,JNDX,k-ind,k)
            if(val_val.ne.0) then
               finder=.true.
            else
               k=k+1
            end if
         end do
         if(finder) then
            test=1
         end if
      end if
      end subroutine  XDETECT_DIAG 

      function XAB_ROW_COL (VAL,BINDX,BJNDX,i,j,sub_ind,lb)
      DCOMPLEX,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX,BJNDX
      integer,intent(in)::i,j,sub_ind,lb
      integer::k,dummy
      logical::finder
      DCOMPLEX:: XAB_ROW_COL
      dummy=lb*lb
      finder=.false.
      k=1
      do while((k.le.size(BINDX)).and.(.not.finder))
         if(BINDX(k).eq.i.and.BJNDX(k).eq.j) then
            finder=.true.
         else
            k=k+1
         end if
      end do
      if(finder) then
         XAB_ROW_COL=VAL(dummy*(k-1)+sub_ind)
      else
         XAB_ROW_COL=0.
      end if
      end function XAB_ROW_COL 

      subroutine XDETECT_BDIAG (VAL,BINDX,BJNDX,ind,BLDA,test,lb)
      implicit none
      DCOMPLEX,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX,BJNDX
      integer,intent(in)::ind,BLDA,lb
      integer,intent(inout)::test
      logical::finder,sub_finder
      DCOMPLEX::val_val
      integer ::k,sub_ind,dummy
      
      dummy=lb*lb
      test=0
      if(ind.eq.0) then ! main diag
         test=0       
         finder=.false.
         k=1
         do while((k.le.BLDA).and.(.not.finder))
            sub_ind=1
            sub_finder=.false.
            do while((sub_ind.le.dummy).and.(.not.sub_finder)) 
            val_val=Ab_row_col(VAL,BINDX,BJNDX,k,k,sub_ind,lb)
               if(val_val.ne.0) then
                  sub_finder=.true.
               else
                  sub_ind= sub_ind+1
               end if
            end do
            if(sub_finder) then
               finder=.true.
            else
               k=k+1
            end if
         end do
         if(finder) then
            test=1
         end if
      elseif(ind.lt.0)then ! low diag
         test=0  
         finder=.false.
         k=-ind+1
         !do while((k.le.BLDA+ind).and.(.not.finder))
         do while((k.le.BLDA).and.(.not.finder))
            sub_ind=1
            sub_finder=.false.
            do while((sub_ind.le.dummy).and.(.not.sub_finder)) 
        val_val=Ab_row_col(VAL,BINDX,BJNDX,k,k+ind,sub_ind,lb)
               if(val_val.ne.0) then
                  sub_finder=.true.
               else
                  sub_ind= sub_ind+1
               end if
            end do
            if(sub_finder) then
               finder=.true.
            else
               k=k+1
            end if
         end do
         if(finder) then
            test=1
         end if
      else ! high diag
         test=0
         finder=.false.
         k=ind+1
         do while((k.le.BLDA).and.(.not.finder))
            sub_ind=1
            sub_finder=.false.
            do while((sub_ind.le.dummy).and.(.not.sub_finder)) 
        val_val=Ab_row_col(VAL,BINDX,BJNDX,k-ind,k,sub_ind,lb)
               if(val_val.ne.0) then
                  sub_finder=.true.
               else
                  sub_ind= sub_ind+1
               end if
            end do
            if(sub_finder) then
               finder=.true.
            else
               k=k+1
            end if
         end do
         if(finder) then
            test=1
         end if
      end if
      end subroutine XDETECT_BDIAG

      subroutine XPRE_COO2DIA (m,n,VAL,INDX,JNDX,LDA,NDIAG)
      implicit none
      DCOMPLEX ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX
      integer,pointer,dimension(:)::JNDX
      integer,dimension(:),allocatable::IDIAG
      DCOMPLEX ,dimension(:),allocatable::VAL_DIA
      integer ,intent(in)::m,n
      integer,intent(inout)::LDA,NDIAG
      integer :: i,test,ind,j,k,IDIAG_ind
      integer :: VAL_ind,VAL_DIA_size
      intrinsic min
      
      LDA=min(m,n)
      test=0
      VAL_ind=0
      IDIAG_ind=0
      NDIAG=0
      ind=0
      call  detect_diag(VAL,INDX,JNDX,ind,LDA,test)
      if(test.eq.1) then
         NDIAG = NDIAG+1
      end if
      do i=1,m-1
         call detect_diag(VAL,INDX,JNDX,-i,LDA,test)
         if(test.eq.1) then
            NDIAG=NDIAG+1
         end if
      end do
      do j=1,n-1
         call detect_diag(VAL,INDX,JNDX,j,LDA,test)
         if(test.eq.1) then
            NDIAG=NDIAG+1
         end if
      end do
      VAL_DIA_size=NDIAG*LDA
      allocate(IDIAG(NDIAG))
      allocate(VAL_DIA(VAL_DIA_size))
      !********* main diag ************
      ind=0
      call detect_diag(VAL,INDX,JNDX,ind,LDA,test)
      if(test.eq.1) then
         do i=1,LDA
            VAL_ind=VAL_ind+1
            VAL_DIA(VAL_ind)=A_row_col(VAL,INDX,JNDX,i,i)
         end do
         IDIAG(IDIAG_ind+1)=0
         IDIAG_ind=IDIAG_ind+1
      end if
      !**********low diag ***********
      do i=1,m-1
         call detect_diag(VAL,INDX,JNDX,-i,LDA,test)
         if(test.eq.1) then 
            do k=1,i
               VAL_ind=VAL_ind+1
               VAL_DIA(VAL_ind)=0
            end do
            do k=i+1,LDA
               VAL_ind=VAL_ind+1
               VAL_DIA(VAL_ind)=A_row_col(VAL,INDX,JNDX,k,k-i)
            end do
            IDIAG(IDIAG_ind+1)=-i
            IDIAG_ind=IDIAG_ind+1
         end if
      end do
      !*********** high diag *******
      do j=1,n-1
         call detect_diag(VAL,INDX,JNDX,j,LDA,test)
         if(test.eq.1) then 
            do   k=1,LDA-j
               VAL_ind=VAL_ind+1
               VAL_DIA(VAL_ind)=A_row_col(VAL,INDX,JNDX,k,k+j)
            end do
            do k=LDA-j+1,LDA
               VAL_ind=VAL_ind+1
               VAL_DIA(VAL_ind)=0
            end do
            IDIAG(IDIAG_ind+1)=j
            IDIAG_ind=IDIAG_ind+1
         end if
      end do
      deallocate(VAL,INDX)
      allocate(VAL(VAL_DIA_size),INDX(NDIAG))
      VAL=VAL_DIA
      INDX=IDIAG
      deallocate(VAL_DIA)
      deallocate(IDIAG)
      end subroutine XPRE_COO2DIA

      subroutine XPRE_DIA2COO(VAL_DIA,IDIAG,IA2,LDA,NNZ)
      implicit none
      DCOMPLEX,pointer,dimension(:) ::VAL_DIA
      integer ,pointer,dimension(:) :: IDIAG,IA2
      integer,intent(in)::LDA,NNZ
      DCOMPLEX ,dimension(:),allocatable::VAL 
      integer,dimension(:),allocatable ::INDX,JNDX
      integer:: VAL_size,i,k
      integer::VAL_ind,IND_ind,VAL_DIA_ind

      VAL_size =NNZ
      allocate(VAL( VAL_size),INDX( VAL_size),JNDX( VAL_size))
      VAL=0.
      INDX=0.
      JNDX=0.
      VAL_ind=0
      IND_ind=0
      VAL_DIA_ind=0
      do i=1,size(IDIAG)
         if(IDIAG(i).lt.0) then
            do k=1,LDA
               VAL_DIA_ind= VAL_DIA_ind+1
               if(VAL_DIA( VAL_DIA_ind).ne.0) then
                  IND_ind=IND_ind+1
                  VAL_ind=VAL_ind+1
                  VAL(VAL_ind)=VAL_DIA( VAL_DIA_ind)
                  INDX(IND_ind)=k
                  JNDX(IND_ind)=k+IDIAG(i)           
               end if
            end do
         elseif(IDIAG(i).eq.0) then
            do k=1,LDA
               VAL_DIA_ind= VAL_DIA_ind+1
               if(VAL_DIA( VAL_DIA_ind).ne.0) then
                  IND_ind=IND_ind+1
                  VAL_ind=VAL_ind+1
                  VAL(VAL_ind)=VAL_DIA( VAL_DIA_ind)
                  INDX(IND_ind)=k
                  JNDX(IND_ind)=k
               end if
            end do
         else
            do k=1,LDA
               VAL_DIA_ind= VAL_DIA_ind+1
               if(VAL_DIA( VAL_DIA_ind).ne.0) then
                  IND_ind=IND_ind+1
                  VAL_ind=VAL_ind+1
                  VAL(VAL_ind)=VAL_DIA(VAL_DIA_ind )
                  JNDX(IND_ind)=IDIAG(i)+k
                  INDX(IND_ind)=k
               end if
            end do
         end if
      end do
      deallocate(VAL_DIA,IDIAG,IA2)
      allocate(VAL_DIA( VAL_size),IDIAG( VAL_size),IA2( VAL_size))
      VAL_DIA=VAL
      IDIAG=INDX
      IA2=JNDX
      deallocate(VAL,INDX,JNDX)

      end subroutine XPRE_DIA2COO
 
      subroutine XPRE_BCO2BDI(mb,kb,lb,VAL,BINDX,BJNDX,BLDA,BNDIAG)
      implicit none
      DCOMPLEX ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX
      integer,pointer,dimension(:)::BJNDX
      integer,dimension(:),allocatable::BIDIAG
      DCOMPLEX ,dimension(:),allocatable::VAL_DIA
      integer ,intent(in)::mb,kb,lb
      integer,intent(inout)::BLDA,BNDIAG
      integer :: i,test,ind,j,k,BIDIAG_ind,VAL_ind
      integer ::VAL_DIA_size,dummy,sub_ind
      intrinsic min
      
      BLDA=min(mb,kb)
      dummy=lb*lb
      test=0
      VAL_ind=0
      BIDIAG_ind=0
      BNDIAG=0
      ind=0
      call  detect_bdiag(VAL,BINDX,BJNDX,ind,BLDA,test,lb)
      if(test.eq.1) then
         BNDIAG = BNDIAG+1
      end if
      do i=1,mb-1
         call detect_bdiag(VAL,BINDX,BJNDX,-i,BLDA,test,lb)
         if(test.eq.1) then
            BNDIAG=BNDIAG+1
         end if
      end do
      do j=1,kb-1
         call detect_bdiag(VAL,BINDX,BJNDX,j,BLDA,test,lb)
         if(test.eq.1) then
            BNDIAG=BNDIAG+1
         end if
      end do
      VAL_DIA_size=BNDIAG*BLDA*dummy
      allocate(BIDIAG(BNDIAG))
      allocate(VAL_DIA(VAL_DIA_size))
      !********* main diag ************
      ind=0
      call detect_bdiag(VAL,BINDX,BJNDX,ind,BLDA,test,lb)
      if(test.eq.1) then
         do i=1,BLDA
            do sub_ind=1,dummy
               VAL_ind=VAL_ind+1
       VAL_DIA(VAL_ind)=Ab_row_col(VAL,BINDX,BJNDX,i,i,sub_ind,lb)
            end do
         end do
         BIDIAG(BIDIAG_ind+1)=0
         BIDIAG_ind=BIDIAG_ind+1
      end if
      !**********low diag ***********
      do i=1,mb-1
         call detect_bdiag(VAL,BINDX,BJNDX,-i,BLDA,test,lb)
         if(test.eq.1) then 
            do k=1,i
               do sub_ind=1,dummy
                  VAL_ind=VAL_ind+1
                  VAL_DIA(VAL_ind)=0
               end do
            end do
            do k=i+1,BLDA
               do sub_ind=1,dummy
                  VAL_ind=VAL_ind+1
      VAL_DIA(VAL_ind)=Ab_row_col(VAL,BINDX,BJNDX,k,k-i,sub_ind,lb)
               end do
            end do
            BIDIAG(BIDIAG_ind+1)=-i
            BIDIAG_ind=BIDIAG_ind+1
         end if
      end do
      !*********** high diag *******
      do j=1,kb-1
         call detect_bdiag(VAL,BINDX,BJNDX,j,BLDA,test,lb)
         if(test.eq.1) then 
            do   k=1,BLDA-j
               do sub_ind=1,dummy
                  VAL_ind=VAL_ind+1
      VAL_DIA(VAL_ind)=Ab_row_col(VAL,BINDX,BJNDX,k,k+j,sub_ind,lb)
               end do
            end do
            do k=BLDA-j+1,BLDA
               do sub_ind=1,dummy
                  VAL_ind=VAL_ind+1
                  VAL_DIA(VAL_ind)=0
               end do
            end do
            BIDIAG(BIDIAG_ind+1)=j
            BIDIAG_ind=BIDIAG_ind+1
         end if
      end do
      deallocate(VAL,BINDX)
      allocate(VAL(VAL_DIA_size),BINDX(BNDIAG))
      VAL=VAL_DIA
      BINDX=BIDIAG
      deallocate(VAL_DIA)
      deallocate(BIDIAG)
      end subroutine XPRE_BCO2BDI
  
      subroutine XPRE_BDI2BCO(VAL_DIA,BIDIAG,IA2,BLDA,BNNZ,lb)
      implicit none
      DCOMPLEX,pointer,dimension(:) ::VAL_DIA
      integer ,pointer,dimension(:) ::BIDIAG,IA2
      integer,intent(in)::BLDA,lb
      integer,intent(out)::BNNZ
      integer ,dimension(:),allocatable::VAL 
      integer,dimension(:),allocatable ::BINDX,BJNDX
      DCOMPLEX::val_val
      integer:: VAL_size,i,k,VAL_ind,IND_ind
      integer::VAL_DIA_ind,dummy,sub_ind,NB_BLOCKS,dummy2
      logical ::sub_finder
      intrinsic floor
      
      dummy=lb*lb
      NB_BLOCKS=floor(real(size(VAL_DIA)/dummy))
      BNNZ=0
      do VAL_DIA_ind=1,NB_BLOCKS
         sub_finder=.false.
         do sub_ind =1,dummy
            val_val=VAL_DIA(dummy*(VAL_DIA_ind-1)+sub_ind)
            if(val_val.ne.0) then 
               sub_finder=.true.
            end if
            if(sub_finder) exit
         end do
         if(sub_finder) then
            BNNZ=BNNZ+1
         end if
      end do
      VAL_size =BNNZ
      allocate(VAL(dummy*VAL_size),BINDX(VAL_size))
      allocate(BJNDX(VAL_size))
      VAL=0.
      BINDX=0
      BJNDX=0
      VAL_ind=0
      IND_ind=0
      VAL_DIA_ind=0

      do i=1,size(BIDIAG)
         if(BIDIAG(i).lt.0) then
            do k=1,BLDA
            VAL_DIA_ind= VAL_DIA_ind+1
            sub_finder=.false.
            do sub_ind =1,dummy
               val_val=VAL_DIA(dummy*(VAL_DIA_ind-1)+sub_ind)
               if(val_val.ne.0) then 
                  sub_finder=.true.
               end if
               if(sub_finder) exit
            end do
            if(sub_finder) then
               IND_ind=IND_ind+1
               VAL_ind=VAL_ind+1
               do sub_ind=1,dummy
                  dummy2=dummy*(VAL_ind-1)+sub_ind
                  VAL(dummy2)=VAL_DIA(dummy2)
               end do
               BINDX(IND_ind)=k            
               BJNDX(IND_ind)=k+BIDIAG(i)           
            end if
         end do
      elseif(BIDIAG(i).eq.0) then
         do k=1,BLDA
            VAL_DIA_ind= VAL_DIA_ind+1
            sub_finder=.false.
            do sub_ind =1,dummy
               val_val=VAL_DIA(dummy*(VAL_DIA_ind-1)+sub_ind)
               if(val_val.ne.0) then 
                  sub_finder=.true.
               end if
               if(sub_finder) exit
            end do
            if(sub_finder) then
               IND_ind=IND_ind+1
               VAL_ind=VAL_ind+1
               do sub_ind=1,dummy
                  dummy2=dummy*(VAL_ind-1)+sub_ind
                  VAL(dummy2)=VAL_DIA(dummy2)
               end do
               BINDX(IND_ind)=k
               BJNDX(IND_ind)=k
            end if
         end do
      else
         do k=1,BLDA
            VAL_DIA_ind= VAL_DIA_ind+1
            sub_finder=.false.
            do sub_ind =1,dummy
               val_val=VAL_DIA(dummy*(VAL_DIA_ind-1)+sub_ind)
               if(val_val.ne.0) then 
                  sub_finder=.true.
               end if
               if(sub_finder) exit
            end do
            if(sub_finder) then
               IND_ind=IND_ind+1
               VAL_ind=VAL_ind+1
               do sub_ind=1,dummy
                  dummy2=dummy*(VAL_ind-1)+sub_ind
                  VAL(dummy2)=VAL_DIA(dummy2)
               end do
               BJNDX(IND_ind)=BIDIAG(i)+k
               BINDX(IND_ind)=k
            end if
         end do
      end if
      end do
      deallocate(VAL_DIA,BIDIAG,IA2)
      allocate(VAL_DIA(dummy*VAL_size),BIDIAG(VAL_size))
      allocate(IA2(VAL_size))
      VAL_DIA=VAL
      BIDIAG=BINDX
      IA2=BJNDX
      deallocate(VAL,BINDX,BJNDX)
      end subroutine XPRE_BDI2BCO
 
      subroutine XPRE_COO2CSR(VAL,INDX,JNDX,M,PNTRB,PNTRE)
      implicit none
      DCOMPLEX ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX
      integer,pointer,dimension(:)::JNDX
      integer,pointer,dimension(:)::PNTRB,PNTRE
      integer ,intent(in)::M
      integer,pointer,dimension(:)::DV,ORD_RES,FNL_RES
      integer :: i,s

      s=size(INDX)
      allocate(DV(M))
      allocate(ORD_RES(s))
      allocate(FNL_RES(s))
      
      do i=1,M
         DV(i)=counter(INDX,i)
      end do
      call  PNTR(PNTRB,PNTRE,M,INDX)
      call up_order(INDX,ORD_RES)
      INDX=JNDX
      VAL=VAL(ORD_RES)
      INDX=INDX(ORD_RES)
      call final_order(INDX,FNL_RES,DV)
      VAL=VAL(FNL_RES)
      INDX=INDX(FNL_RES)
      deallocate(DV)
      deallocate(ORD_RES)
      deallocate(FNL_RES)
      end subroutine  XPRE_COO2CSR

      subroutine XPRE_COO2CSC(VAL,INDX,JNDX,K,PNTRB,PNTRE)
      implicit none
      DCOMPLEX ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX
      integer,pointer,dimension(:)::JNDX
      integer,pointer,dimension(:)::PNTRB,PNTRE
      integer ,intent(in)::K
      integer,pointer,dimension(:)::DV,ORD_RES,FNL_RES
      integer :: i,s

      s=size(JNDX)
      allocate(DV(K))
      allocate(ORD_RES(s))
      allocate(FNL_RES(s))

      do i=1,K
         DV(i)=counter(JNDX,i)
      end do
      allocate(DV(K+1))
      call  PNTR(PNTRB,PNTRE,K,JNDX)
      call up_order(JNDX,ORD_RES)
      VAL=VAL(ORD_RES)
      INDX=INDX(ORD_RES)
      call final_order(INDX,FNL_RES,DV)
      VAL=VAL(FNL_RES)
      INDX=INDX(FNL_RES)
      deallocate(DV)
      deallocate(ORD_RES)
      deallocate(FNL_RES)
      end subroutine XPRE_COO2CSC

      subroutine XPRE_BCO2BSR(VAL,BINDX,BJNDX,MB,LB,BPNTRB,BPNTRE)
      implicit none
      DCOMPLEX ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX
      integer,pointer,dimension(:)::BJNDX
      integer,pointer,dimension(:)::BPNTRB,BPNTRE
      integer ,intent(in)::MB,LB
      integer,pointer,dimension(:)::DV,ORD_RES,FNL_RES
      integer :: i,s,dummy
     
      s=size(BINDX)
      dummy=LB*LB
      allocate(DV(MB))
      allocate(ORD_RES(s))
      allocate(FNL_RES(s))

      do i=1,MB
         DV(i)=counter(BINDX,i)
      end do
      call  PNTR(BPNTRB,BPNTRE,MB,BINDX)
      call up_order(BINDX,ORD_RES)
      BINDX=BJNDX
      call  XB_UP_ORDER(VAL,dummy,ORD_RES) 
      BINDX=BINDX(ORD_RES)
      call final_order(BINDX,FNL_RES,DV)
      call  XB_UP_ORDER(VAL,dummy,FNL_RES)
      BINDX=BINDX(FNL_RES)
      deallocate(DV)
      deallocate(ORD_RES)
      deallocate(FNL_RES)
      end subroutine XPRE_BCO2BSR

      subroutine XPRE_BCO2BSC(VAL,BINDX,BJNDX,KB,LB,BPNTRB,BPNTRE)
      implicit none
      DCOMPLEX ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX
      integer,pointer,dimension(:)::BJNDX
      integer,pointer,dimension(:)::BPNTRB,BPNTRE
      integer ,intent(in)::KB,LB
      integer,pointer,dimension(:)::DV,ORD_RES,FNL_RES
      integer :: i,s,dummy

      s=size(BJNDX)
      dummy=LB*LB
      allocate(DV(KB))
      allocate(ORD_RES(s))
      allocate(FNL_RES(s))
      do i=1,KB
         DV(i)=counter(BJNDX,i)
      end do
      DV(KB+1)=counter(BJNDX,-1)
      call  PNTR(BPNTRB,BPNTRE,KB,BJNDX)
      call up_order(BJNDX,ORD_RES)
      call  XB_UP_ORDER(VAL,dummy,ORD_RES) 
      BINDX=BINDX(ORD_RES)
      call final_order(BINDX,FNL_RES,DV)
      call  XB_UP_ORDER(VAL,dummy,FNL_RES)
      BINDX=BINDX(FNL_RES)
      deallocate(DV)
      deallocate(ORD_RES)
      deallocate(FNL_RES)
      end subroutine XPRE_BCO2BSC

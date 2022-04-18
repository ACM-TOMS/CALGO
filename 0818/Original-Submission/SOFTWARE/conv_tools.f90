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
      subroutine ib_up_order (VAL,lbxlb,BINDX)
      implicit none
      integer ,pointer  ,dimension(:)::VAL
      integer ,pointer, dimension(:)::BINDX
      integer ,dimension(:),allocatable :: tes
      integer ,dimension(:,:),allocatable::P
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
      end subroutine ib_up_order  
      function iA_row_col (VAL,INDX,JNDX,i,j)
      integer ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX,JNDX
      integer,intent(in)::i,j
      integer::k
      logical::finder
      integer :: iA_row_col  
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
         iA_row_col  =VAL(k)
      else
         iA_row_col  =0
      end if
      end function iA_row_col 
      subroutine idetect_diag (VAL,INDX,JNDX,ind,LDA,test)
      implicit none
      integer ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX,JNDX
      integer,intent(in)::ind,LDA
      integer,intent(inout)::test
      logical::finder
      integer ::val_val
      integer ::k
      test=0
      if(ind.eq.0) then ! main diag
         test=0       
         finder=.false.
         k=1
         do while((k.le.LDA).and.(.not.finder))
            val_val= iA_row_col (VAL,INDX,JNDX,k,k)
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
      end subroutine  idetect_diag  
      function iAb_row_col  (VAL,BINDX,BJNDX,i,j,sub_ind,lb)
      integer ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX,BJNDX
      integer,intent(in)::i,j,sub_ind,lb
      integer::k,dummy
      logical::finder
      integer :: iAb_row_col 
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
         iAb_row_col =VAL(dummy*(k-1)+sub_ind)
      else
         iAb_row_col =0.
      end if
      end function iAb_row_col  
      subroutine idetect_bdiag  (VAL,BINDX,BJNDX,ind,BLDA,test,lb)
      implicit none
      integer ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX,BJNDX
      integer,intent(in)::ind,BLDA,lb
      integer,intent(inout)::test
      logical::finder,sub_finder
      integer ::val_val
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
      end subroutine idetect_bdiag 
      subroutine ipre_usconv_coo2dia  (m,n,VAL,INDX,JNDX,LDA,NDIAG)
      implicit none
      integer  ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX
      integer,pointer,dimension(:)::JNDX
      integer,dimension(:),allocatable::IDIAG
      integer  ,dimension(:),allocatable::VAL_DIA
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
      end subroutine ipre_usconv_coo2dia 
      subroutine ipre_usconv_dia2coo (VAL_DIA,IDIAG,IA2,LDA,NNZ)
      implicit none
      integer ,pointer,dimension(:) ::VAL_DIA
      integer ,pointer,dimension(:) :: IDIAG,IA2
      integer,intent(in)::LDA,NNZ
      integer  ,dimension(:),allocatable::VAL 
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
      end subroutine ipre_usconv_dia2coo 
      subroutine ipre_usconv_bco2bdi (mb,kb,lb,VAL,BINDX,BJNDX,BLDA,BNDIAG)
      implicit none
      integer  ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX
      integer,pointer,dimension(:)::BJNDX
      integer,dimension(:),allocatable::BIDIAG
      integer  ,dimension(:),allocatable::VAL_DIA
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
      end subroutine ipre_usconv_bco2bdi 
      subroutine ipre_usconv_bdi2bco (VAL_DIA,BIDIAG,IA2,BLDA,BNNZ,lb)
      implicit none
      integer ,pointer,dimension(:) ::VAL_DIA
      integer ,pointer,dimension(:) ::BIDIAG,IA2
      integer,intent(in)::BLDA,lb
      integer,intent(out)::BNNZ
      integer ,dimension(:),allocatable::VAL 
      integer,dimension(:),allocatable ::BINDX,BJNDX
      integer ::val_val
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
      end subroutine ipre_usconv_bdi2bco 
      subroutine ipre_usconv_coo2csr (VAL,INDX,JNDX,M,PNTRB,PNTRE)
      implicit none
      integer  ,pointer,dimension(:)::VAL
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
      end subroutine  ipre_usconv_coo2csr 
      subroutine ipre_usconv_coo2csc (VAL,INDX,JNDX,K,PNTRB,PNTRE)
      implicit none
      integer  ,pointer,dimension(:)::VAL
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
      end subroutine ipre_usconv_coo2csc 
      subroutine ipre_usconv_bco2bsr (VAL,BINDX,BJNDX,MB,LB,BPNTRB,BPNTRE)
      implicit none
      integer  ,pointer,dimension(:)::VAL
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
      call  ib_up_order (VAL,dummy,ORD_RES) 
      BINDX=BINDX(ORD_RES)
      call final_order(BINDX,FNL_RES,DV)
      call  ib_up_order (VAL,dummy,FNL_RES)
      BINDX=BINDX(FNL_RES)
      deallocate(DV)
      deallocate(ORD_RES)
      deallocate(FNL_RES)
      end subroutine ipre_usconv_bco2bsr 
      subroutine ipre_usconv_bco2bsc (VAL,BINDX,BJNDX,KB,LB,BPNTRB,BPNTRE)
      implicit none
      integer  ,pointer,dimension(:)::VAL
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
      call  ib_up_order (VAL,dummy,ORD_RES) 
      BINDX=BINDX(ORD_RES)
      call final_order(BINDX,FNL_RES,DV)
      call  ib_up_order (VAL,dummy,FNL_RES)
      BINDX=BINDX(FNL_RES)
      deallocate(DV)
      deallocate(ORD_RES)
      deallocate(FNL_RES)
      end subroutine ipre_usconv_bco2bsc 
      subroutine sb_up_order (VAL,lbxlb,BINDX)
      implicit none
      real(KIND=sp) ,pointer  ,dimension(:)::VAL
      integer ,pointer, dimension(:)::BINDX
      integer ,dimension(:),allocatable :: tes
      real(KIND=sp) ,dimension(:,:),allocatable::P
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
      end subroutine sb_up_order  
      function sA_row_col (VAL,INDX,JNDX,i,j)
      real(KIND=sp) ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX,JNDX
      integer,intent(in)::i,j
      integer::k
      logical::finder
      real(KIND=sp) :: sA_row_col  
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
         sA_row_col  =VAL(k)
      else
         sA_row_col  =0
      end if
      end function sA_row_col 
      subroutine sdetect_diag (VAL,INDX,JNDX,ind,LDA,test)
      implicit none
      real(KIND=sp) ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX,JNDX
      integer,intent(in)::ind,LDA
      integer,intent(inout)::test
      logical::finder
      real(KIND=sp) ::val_val
      integer ::k
      test=0
      if(ind.eq.0) then ! main diag
         test=0       
         finder=.false.
         k=1
         do while((k.le.LDA).and.(.not.finder))
            val_val= sA_row_col (VAL,INDX,JNDX,k,k)
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
      end subroutine  sdetect_diag  
      function sAb_row_col  (VAL,BINDX,BJNDX,i,j,sub_ind,lb)
      real(KIND=sp) ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX,BJNDX
      integer,intent(in)::i,j,sub_ind,lb
      integer::k,dummy
      logical::finder
      real(KIND=sp) :: sAb_row_col 
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
         sAb_row_col =VAL(dummy*(k-1)+sub_ind)
      else
         sAb_row_col =0.
      end if
      end function sAb_row_col  
      subroutine sdetect_bdiag  (VAL,BINDX,BJNDX,ind,BLDA,test,lb)
      implicit none
      real(KIND=sp) ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX,BJNDX
      integer,intent(in)::ind,BLDA,lb
      integer,intent(inout)::test
      logical::finder,sub_finder
      real(KIND=sp) ::val_val
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
      end subroutine sdetect_bdiag 
      subroutine spre_usconv_coo2dia  (m,n,VAL,INDX,JNDX,LDA,NDIAG)
      implicit none
      real(KIND=sp)  ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX
      integer,pointer,dimension(:)::JNDX
      integer,dimension(:),allocatable::IDIAG
      real(KIND=sp)  ,dimension(:),allocatable::VAL_DIA
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
      end subroutine spre_usconv_coo2dia 
      subroutine spre_usconv_dia2coo (VAL_DIA,IDIAG,IA2,LDA,NNZ)
      implicit none
      real(KIND=sp) ,pointer,dimension(:) ::VAL_DIA
      integer ,pointer,dimension(:) :: IDIAG,IA2
      integer,intent(in)::LDA,NNZ
      real(KIND=sp)  ,dimension(:),allocatable::VAL 
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
      end subroutine spre_usconv_dia2coo 
      subroutine spre_usconv_bco2bdi (mb,kb,lb,VAL,BINDX,BJNDX,BLDA,BNDIAG)
      implicit none
      real(KIND=sp)  ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX
      integer,pointer,dimension(:)::BJNDX
      integer,dimension(:),allocatable::BIDIAG
      real(KIND=sp)  ,dimension(:),allocatable::VAL_DIA
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
      end subroutine spre_usconv_bco2bdi 
      subroutine spre_usconv_bdi2bco (VAL_DIA,BIDIAG,IA2,BLDA,BNNZ,lb)
      implicit none
      real(KIND=sp) ,pointer,dimension(:) ::VAL_DIA
      integer ,pointer,dimension(:) ::BIDIAG,IA2
      integer,intent(in)::BLDA,lb
      integer,intent(out)::BNNZ
      integer ,dimension(:),allocatable::VAL 
      integer,dimension(:),allocatable ::BINDX,BJNDX
      real(KIND=sp) ::val_val
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
      end subroutine spre_usconv_bdi2bco 
      subroutine spre_usconv_coo2csr (VAL,INDX,JNDX,M,PNTRB,PNTRE)
      implicit none
      real(KIND=sp)  ,pointer,dimension(:)::VAL
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
      end subroutine  spre_usconv_coo2csr 
      subroutine spre_usconv_coo2csc (VAL,INDX,JNDX,K,PNTRB,PNTRE)
      implicit none
      real(KIND=sp)  ,pointer,dimension(:)::VAL
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
      end subroutine spre_usconv_coo2csc 
      subroutine spre_usconv_bco2bsr (VAL,BINDX,BJNDX,MB,LB,BPNTRB,BPNTRE)
      implicit none
      real(KIND=sp)  ,pointer,dimension(:)::VAL
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
      call  sb_up_order (VAL,dummy,ORD_RES) 
      BINDX=BINDX(ORD_RES)
      call final_order(BINDX,FNL_RES,DV)
      call  sb_up_order (VAL,dummy,FNL_RES)
      BINDX=BINDX(FNL_RES)
      deallocate(DV)
      deallocate(ORD_RES)
      deallocate(FNL_RES)
      end subroutine spre_usconv_bco2bsr 
      subroutine spre_usconv_bco2bsc (VAL,BINDX,BJNDX,KB,LB,BPNTRB,BPNTRE)
      implicit none
      real(KIND=sp)  ,pointer,dimension(:)::VAL
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
      call  sb_up_order (VAL,dummy,ORD_RES) 
      BINDX=BINDX(ORD_RES)
      call final_order(BINDX,FNL_RES,DV)
      call  sb_up_order (VAL,dummy,FNL_RES)
      BINDX=BINDX(FNL_RES)
      deallocate(DV)
      deallocate(ORD_RES)
      deallocate(FNL_RES)
      end subroutine spre_usconv_bco2bsc 
      subroutine db_up_order (VAL,lbxlb,BINDX)
      implicit none
      real(KIND=dp) ,pointer  ,dimension(:)::VAL
      integer ,pointer, dimension(:)::BINDX
      integer ,dimension(:),allocatable :: tes
      real(KIND=dp) ,dimension(:,:),allocatable::P
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
      end subroutine db_up_order  
      function dA_row_col (VAL,INDX,JNDX,i,j)
      real(KIND=dp) ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX,JNDX
      integer,intent(in)::i,j
      integer::k
      logical::finder
      real(KIND=dp) :: dA_row_col  
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
         dA_row_col  =VAL(k)
      else
         dA_row_col  =0
      end if
      end function dA_row_col 
      subroutine ddetect_diag (VAL,INDX,JNDX,ind,LDA,test)
      implicit none
      real(KIND=dp) ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX,JNDX
      integer,intent(in)::ind,LDA
      integer,intent(inout)::test
      logical::finder
      real(KIND=dp) ::val_val
      integer ::k
      test=0
      if(ind.eq.0) then ! main diag
         test=0       
         finder=.false.
         k=1
         do while((k.le.LDA).and.(.not.finder))
            val_val= dA_row_col (VAL,INDX,JNDX,k,k)
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
      end subroutine  ddetect_diag  
      function dAb_row_col  (VAL,BINDX,BJNDX,i,j,sub_ind,lb)
      real(KIND=dp) ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX,BJNDX
      integer,intent(in)::i,j,sub_ind,lb
      integer::k,dummy
      logical::finder
      real(KIND=dp) :: dAb_row_col 
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
         dAb_row_col =VAL(dummy*(k-1)+sub_ind)
      else
         dAb_row_col =0.
      end if
      end function dAb_row_col  
      subroutine ddetect_bdiag  (VAL,BINDX,BJNDX,ind,BLDA,test,lb)
      implicit none
      real(KIND=dp) ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX,BJNDX
      integer,intent(in)::ind,BLDA,lb
      integer,intent(inout)::test
      logical::finder,sub_finder
      real(KIND=dp) ::val_val
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
      end subroutine ddetect_bdiag 
      subroutine dpre_usconv_coo2dia  (m,n,VAL,INDX,JNDX,LDA,NDIAG)
      implicit none
      real(KIND=dp)  ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX
      integer,pointer,dimension(:)::JNDX
      integer,dimension(:),allocatable::IDIAG
      real(KIND=dp)  ,dimension(:),allocatable::VAL_DIA
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
      end subroutine dpre_usconv_coo2dia 
      subroutine dpre_usconv_dia2coo (VAL_DIA,IDIAG,IA2,LDA,NNZ)
      implicit none
      real(KIND=dp) ,pointer,dimension(:) ::VAL_DIA
      integer ,pointer,dimension(:) :: IDIAG,IA2
      integer,intent(in)::LDA,NNZ
      real(KIND=dp)  ,dimension(:),allocatable::VAL 
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
      end subroutine dpre_usconv_dia2coo 
      subroutine dpre_usconv_bco2bdi (mb,kb,lb,VAL,BINDX,BJNDX,BLDA,BNDIAG)
      implicit none
      real(KIND=dp)  ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX
      integer,pointer,dimension(:)::BJNDX
      integer,dimension(:),allocatable::BIDIAG
      real(KIND=dp)  ,dimension(:),allocatable::VAL_DIA
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
      end subroutine dpre_usconv_bco2bdi 
      subroutine dpre_usconv_bdi2bco (VAL_DIA,BIDIAG,IA2,BLDA,BNNZ,lb)
      implicit none
      real(KIND=dp) ,pointer,dimension(:) ::VAL_DIA
      integer ,pointer,dimension(:) ::BIDIAG,IA2
      integer,intent(in)::BLDA,lb
      integer,intent(out)::BNNZ
      integer ,dimension(:),allocatable::VAL 
      integer,dimension(:),allocatable ::BINDX,BJNDX
      real(KIND=dp) ::val_val
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
      end subroutine dpre_usconv_bdi2bco 
      subroutine dpre_usconv_coo2csr (VAL,INDX,JNDX,M,PNTRB,PNTRE)
      implicit none
      real(KIND=dp)  ,pointer,dimension(:)::VAL
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
      end subroutine  dpre_usconv_coo2csr 
      subroutine dpre_usconv_coo2csc (VAL,INDX,JNDX,K,PNTRB,PNTRE)
      implicit none
      real(KIND=dp)  ,pointer,dimension(:)::VAL
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
      end subroutine dpre_usconv_coo2csc 
      subroutine dpre_usconv_bco2bsr (VAL,BINDX,BJNDX,MB,LB,BPNTRB,BPNTRE)
      implicit none
      real(KIND=dp)  ,pointer,dimension(:)::VAL
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
      call  db_up_order (VAL,dummy,ORD_RES) 
      BINDX=BINDX(ORD_RES)
      call final_order(BINDX,FNL_RES,DV)
      call  db_up_order (VAL,dummy,FNL_RES)
      BINDX=BINDX(FNL_RES)
      deallocate(DV)
      deallocate(ORD_RES)
      deallocate(FNL_RES)
      end subroutine dpre_usconv_bco2bsr 
      subroutine dpre_usconv_bco2bsc (VAL,BINDX,BJNDX,KB,LB,BPNTRB,BPNTRE)
      implicit none
      real(KIND=dp)  ,pointer,dimension(:)::VAL
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
      call  db_up_order (VAL,dummy,ORD_RES) 
      BINDX=BINDX(ORD_RES)
      call final_order(BINDX,FNL_RES,DV)
      call  db_up_order (VAL,dummy,FNL_RES)
      BINDX=BINDX(FNL_RES)
      deallocate(DV)
      deallocate(ORD_RES)
      deallocate(FNL_RES)
      end subroutine dpre_usconv_bco2bsc 
      subroutine cb_up_order (VAL,lbxlb,BINDX)
      implicit none
      complex(KIND=sp) ,pointer  ,dimension(:)::VAL
      integer ,pointer, dimension(:)::BINDX
      integer ,dimension(:),allocatable :: tes
      complex(KIND=sp) ,dimension(:,:),allocatable::P
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
      end subroutine cb_up_order  
      function cA_row_col (VAL,INDX,JNDX,i,j)
      complex(KIND=sp) ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX,JNDX
      integer,intent(in)::i,j
      integer::k
      logical::finder
      complex(KIND=sp) :: cA_row_col  
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
         cA_row_col  =VAL(k)
      else
         cA_row_col  =0
      end if
      end function cA_row_col 
      subroutine cdetect_diag (VAL,INDX,JNDX,ind,LDA,test)
      implicit none
      complex(KIND=sp) ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX,JNDX
      integer,intent(in)::ind,LDA
      integer,intent(inout)::test
      logical::finder
      complex(KIND=sp) ::val_val
      integer ::k
      test=0
      if(ind.eq.0) then ! main diag
         test=0       
         finder=.false.
         k=1
         do while((k.le.LDA).and.(.not.finder))
            val_val= cA_row_col (VAL,INDX,JNDX,k,k)
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
      end subroutine  cdetect_diag  
      function cAb_row_col  (VAL,BINDX,BJNDX,i,j,sub_ind,lb)
      complex(KIND=sp) ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX,BJNDX
      integer,intent(in)::i,j,sub_ind,lb
      integer::k,dummy
      logical::finder
      complex(KIND=sp) :: cAb_row_col 
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
         cAb_row_col =VAL(dummy*(k-1)+sub_ind)
      else
         cAb_row_col =0.
      end if
      end function cAb_row_col  
      subroutine cdetect_bdiag  (VAL,BINDX,BJNDX,ind,BLDA,test,lb)
      implicit none
      complex(KIND=sp) ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX,BJNDX
      integer,intent(in)::ind,BLDA,lb
      integer,intent(inout)::test
      logical::finder,sub_finder
      complex(KIND=sp) ::val_val
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
      end subroutine cdetect_bdiag 
      subroutine cpre_usconv_coo2dia  (m,n,VAL,INDX,JNDX,LDA,NDIAG)
      implicit none
      complex(KIND=sp)  ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX
      integer,pointer,dimension(:)::JNDX
      integer,dimension(:),allocatable::IDIAG
      complex(KIND=sp)  ,dimension(:),allocatable::VAL_DIA
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
      end subroutine cpre_usconv_coo2dia 
      subroutine cpre_usconv_dia2coo (VAL_DIA,IDIAG,IA2,LDA,NNZ)
      implicit none
      complex(KIND=sp) ,pointer,dimension(:) ::VAL_DIA
      integer ,pointer,dimension(:) :: IDIAG,IA2
      integer,intent(in)::LDA,NNZ
      complex(KIND=sp)  ,dimension(:),allocatable::VAL 
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
      end subroutine cpre_usconv_dia2coo 
      subroutine cpre_usconv_bco2bdi (mb,kb,lb,VAL,BINDX,BJNDX,BLDA,BNDIAG)
      implicit none
      complex(KIND=sp)  ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX
      integer,pointer,dimension(:)::BJNDX
      integer,dimension(:),allocatable::BIDIAG
      complex(KIND=sp)  ,dimension(:),allocatable::VAL_DIA
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
      end subroutine cpre_usconv_bco2bdi 
      subroutine cpre_usconv_bdi2bco (VAL_DIA,BIDIAG,IA2,BLDA,BNNZ,lb)
      implicit none
      complex(KIND=sp) ,pointer,dimension(:) ::VAL_DIA
      integer ,pointer,dimension(:) ::BIDIAG,IA2
      integer,intent(in)::BLDA,lb
      integer,intent(out)::BNNZ
      integer ,dimension(:),allocatable::VAL 
      integer,dimension(:),allocatable ::BINDX,BJNDX
      complex(KIND=sp) ::val_val
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
      end subroutine cpre_usconv_bdi2bco 
      subroutine cpre_usconv_coo2csr (VAL,INDX,JNDX,M,PNTRB,PNTRE)
      implicit none
      complex(KIND=sp)  ,pointer,dimension(:)::VAL
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
      end subroutine  cpre_usconv_coo2csr 
      subroutine cpre_usconv_coo2csc (VAL,INDX,JNDX,K,PNTRB,PNTRE)
      implicit none
      complex(KIND=sp)  ,pointer,dimension(:)::VAL
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
      end subroutine cpre_usconv_coo2csc 
      subroutine cpre_usconv_bco2bsr (VAL,BINDX,BJNDX,MB,LB,BPNTRB,BPNTRE)
      implicit none
      complex(KIND=sp)  ,pointer,dimension(:)::VAL
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
      call  cb_up_order (VAL,dummy,ORD_RES) 
      BINDX=BINDX(ORD_RES)
      call final_order(BINDX,FNL_RES,DV)
      call  cb_up_order (VAL,dummy,FNL_RES)
      BINDX=BINDX(FNL_RES)
      deallocate(DV)
      deallocate(ORD_RES)
      deallocate(FNL_RES)
      end subroutine cpre_usconv_bco2bsr 
      subroutine cpre_usconv_bco2bsc (VAL,BINDX,BJNDX,KB,LB,BPNTRB,BPNTRE)
      implicit none
      complex(KIND=sp)  ,pointer,dimension(:)::VAL
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
      call  cb_up_order (VAL,dummy,ORD_RES) 
      BINDX=BINDX(ORD_RES)
      call final_order(BINDX,FNL_RES,DV)
      call  cb_up_order (VAL,dummy,FNL_RES)
      BINDX=BINDX(FNL_RES)
      deallocate(DV)
      deallocate(ORD_RES)
      deallocate(FNL_RES)
      end subroutine cpre_usconv_bco2bsc 
      subroutine zb_up_order (VAL,lbxlb,BINDX)
      implicit none
      complex(KIND=dp) ,pointer  ,dimension(:)::VAL
      integer ,pointer, dimension(:)::BINDX
      integer ,dimension(:),allocatable :: tes
      complex(KIND=dp) ,dimension(:,:),allocatable::P
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
      end subroutine zb_up_order  
      function zA_row_col (VAL,INDX,JNDX,i,j)
      complex(KIND=dp) ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX,JNDX
      integer,intent(in)::i,j
      integer::k
      logical::finder
      complex(KIND=dp) :: zA_row_col  
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
         zA_row_col  =VAL(k)
      else
         zA_row_col  =0
      end if
      end function zA_row_col 
      subroutine zdetect_diag (VAL,INDX,JNDX,ind,LDA,test)
      implicit none
      complex(KIND=dp) ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX,JNDX
      integer,intent(in)::ind,LDA
      integer,intent(inout)::test
      logical::finder
      complex(KIND=dp) ::val_val
      integer ::k
      test=0
      if(ind.eq.0) then ! main diag
         test=0       
         finder=.false.
         k=1
         do while((k.le.LDA).and.(.not.finder))
            val_val= zA_row_col (VAL,INDX,JNDX,k,k)
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
      end subroutine  zdetect_diag  
      function zAb_row_col  (VAL,BINDX,BJNDX,i,j,sub_ind,lb)
      complex(KIND=dp) ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX,BJNDX
      integer,intent(in)::i,j,sub_ind,lb
      integer::k,dummy
      logical::finder
      complex(KIND=dp) :: zAb_row_col 
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
         zAb_row_col =VAL(dummy*(k-1)+sub_ind)
      else
         zAb_row_col =0.
      end if
      end function zAb_row_col  
      subroutine zdetect_bdiag  (VAL,BINDX,BJNDX,ind,BLDA,test,lb)
      implicit none
      complex(KIND=dp) ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX,BJNDX
      integer,intent(in)::ind,BLDA,lb
      integer,intent(inout)::test
      logical::finder,sub_finder
      complex(KIND=dp) ::val_val
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
      end subroutine zdetect_bdiag 
      subroutine zpre_usconv_coo2dia  (m,n,VAL,INDX,JNDX,LDA,NDIAG)
      implicit none
      complex(KIND=dp)  ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: INDX
      integer,pointer,dimension(:)::JNDX
      integer,dimension(:),allocatable::IDIAG
      complex(KIND=dp)  ,dimension(:),allocatable::VAL_DIA
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
      end subroutine zpre_usconv_coo2dia 
      subroutine zpre_usconv_dia2coo (VAL_DIA,IDIAG,IA2,LDA,NNZ)
      implicit none
      complex(KIND=dp) ,pointer,dimension(:) ::VAL_DIA
      integer ,pointer,dimension(:) :: IDIAG,IA2
      integer,intent(in)::LDA,NNZ
      complex(KIND=dp)  ,dimension(:),allocatable::VAL 
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
      end subroutine zpre_usconv_dia2coo 
      subroutine zpre_usconv_bco2bdi (mb,kb,lb,VAL,BINDX,BJNDX,BLDA,BNDIAG)
      implicit none
      complex(KIND=dp)  ,pointer,dimension(:)::VAL
      integer ,pointer,dimension(:):: BINDX
      integer,pointer,dimension(:)::BJNDX
      integer,dimension(:),allocatable::BIDIAG
      complex(KIND=dp)  ,dimension(:),allocatable::VAL_DIA
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
      end subroutine zpre_usconv_bco2bdi 
      subroutine zpre_usconv_bdi2bco (VAL_DIA,BIDIAG,IA2,BLDA,BNNZ,lb)
      implicit none
      complex(KIND=dp) ,pointer,dimension(:) ::VAL_DIA
      integer ,pointer,dimension(:) ::BIDIAG,IA2
      integer,intent(in)::BLDA,lb
      integer,intent(out)::BNNZ
      integer ,dimension(:),allocatable::VAL 
      integer,dimension(:),allocatable ::BINDX,BJNDX
      complex(KIND=dp) ::val_val
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
      end subroutine zpre_usconv_bdi2bco 
      subroutine zpre_usconv_coo2csr (VAL,INDX,JNDX,M,PNTRB,PNTRE)
      implicit none
      complex(KIND=dp)  ,pointer,dimension(:)::VAL
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
      end subroutine  zpre_usconv_coo2csr 
      subroutine zpre_usconv_coo2csc (VAL,INDX,JNDX,K,PNTRB,PNTRE)
      implicit none
      complex(KIND=dp)  ,pointer,dimension(:)::VAL
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
      end subroutine zpre_usconv_coo2csc 
      subroutine zpre_usconv_bco2bsr (VAL,BINDX,BJNDX,MB,LB,BPNTRB,BPNTRE)
      implicit none
      complex(KIND=dp)  ,pointer,dimension(:)::VAL
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
      call  zb_up_order (VAL,dummy,ORD_RES) 
      BINDX=BINDX(ORD_RES)
      call final_order(BINDX,FNL_RES,DV)
      call  zb_up_order (VAL,dummy,FNL_RES)
      BINDX=BINDX(FNL_RES)
      deallocate(DV)
      deallocate(ORD_RES)
      deallocate(FNL_RES)
      end subroutine zpre_usconv_bco2bsr 
      subroutine zpre_usconv_bco2bsc (VAL,BINDX,BJNDX,KB,LB,BPNTRB,BPNTRE)
      implicit none
      complex(KIND=dp)  ,pointer,dimension(:)::VAL
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
      call  zb_up_order (VAL,dummy,ORD_RES) 
      BINDX=BINDX(ORD_RES)
      call final_order(BINDX,FNL_RES,DV)
      call  zb_up_order (VAL,dummy,FNL_RES)
      BINDX=BINDX(FNL_RES)
      deallocate(DV)
      deallocate(ORD_RES)
      deallocate(FNL_RES)
      end subroutine zpre_usconv_bco2bsc 
      end module mod_conv_tools

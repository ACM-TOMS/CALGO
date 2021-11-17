      module mod_uscr_bsr
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : CREATION ROUTINE FOR MATRIX HANDLE FROM 'BSR'-FORMAT
! **********************************************************************
      use representation_of_data
      use properties
      use mod_usds
      implicit none
      interface uscr_bsr
        module procedure iuscr_bsr
        module procedure suscr_bsr
        module procedure duscr_bsr
        module procedure cuscr_bsr
        module procedure zuscr_bsr
      end interface
      contains
! **********************************************************************
! **********************************************************************
      subroutine iuscr_bsr (m,n,val,bindx,bpntrb,bpntre,mb,kb,lb,&
                           prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,mb,kb,lb,prpty
      integer, dimension(:), intent(inout),target :: bindx,bpntrb,bpntre
      integer , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr, options, base , nnz
      logical :: COPY
      character :: message
      type(ispmat ),pointer :: dsp_data
      options = istat
      istat = -1 !if not changed later, routine has failed
      a = 0 
      call new_isp (nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      call accessdata(dsp_data,nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      dsp_data%FIDA = 'BSR'
      dsp_data%M = m
      dsp_data%K = n
      call set_descra(dsp_data%DESCRA,prpty,ierr)
      call get_descra(dsp_data%DESCRA,'b',message,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if      
      if (message.eq.'C') then
         base = C_BASE
      else !Assuming F base
         base = F_BASE
      end if
      call set_infoa(dsp_data%INFOA,'b',base,ierr)
      nnz = size(bindx) !no. of nonzero blocks
      call set_infoa(dsp_data%INFOA,'n',nnz,ierr)
      call set_infoa(dsp_data%INFOA,'d',lb,ierr) !row-dim of a block
      call set_infoa(dsp_data%INFOA,'e',lb,ierr) !col-dim of a block
      call set_infoa(dsp_data%INFOA,'f',mb,ierr) !row-dim in blocks
      call set_infoa(dsp_data%INFOA,'g',kb,ierr) !col-dim in blocks
      if((mb.ne.size(bpntrb)).or.(m.ne.mb*lb).or.(n.ne.kb*lb).or.&
        (minval(bindx).lt.base).or.(maxval(bindx).gt.kb-1+base).or.&
        (mb.ne.size(bpntre)).or.(nnz*lb*lb.ne.size(val))) then
        call usds(nmb,ierr)
        ierr = blas_error_param
        return
      end if
      if (options.gt.0) then
         ! decision rule whether or not to copy 
         COPY = .TRUE.
         if(COPY) then 
            options = -1 !copy
         else
            options = 0  !reference
         end if
      end if
      if (options.eq.0) then
         call set_infoa(dsp_data%INFOA,'c',REF_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! create reference to original matrix
         dsp_data%A => val
         dsp_data%IA1 => bindx
         dsp_data%pb => bpntrb
         dsp_data%pe => bpntre
         istat = 0
      else
         ! The additional required memory is DEALLOCATED later in USDS!
         call set_infoa(dsp_data%INFOA,'c',COP_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! copy original data
         allocate(dsp_data%A(size(val)),dsp_data%IA1(size(bindx)),&
                  dsp_data%pb(size(bpntrb)),dsp_data%pe(size(bpntre))&
                  ,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = bindx
         dsp_data%pb = bpntrb
         dsp_data%pe = bpntre
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine iuscr_bsr 
! **********************************************************************
! **********************************************************************
      subroutine suscr_bsr (m,n,val,bindx,bpntrb,bpntre,mb,kb,lb,&
                           prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,mb,kb,lb,prpty
      integer, dimension(:), intent(inout),target :: bindx,bpntrb,bpntre
      real(KIND=sp) , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr, options, base , nnz
      logical :: COPY
      character :: message
      type(sspmat ),pointer :: dsp_data
      options = istat
      istat = -1 !if not changed later, routine has failed
      a = 0 
      call new_ssp (nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      call accessdata(dsp_data,nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      dsp_data%FIDA = 'BSR'
      dsp_data%M = m
      dsp_data%K = n
      call set_descra(dsp_data%DESCRA,prpty,ierr)
      call get_descra(dsp_data%DESCRA,'b',message,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if      
      if (message.eq.'C') then
         base = C_BASE
      else !Assuming F base
         base = F_BASE
      end if
      call set_infoa(dsp_data%INFOA,'b',base,ierr)
      nnz = size(bindx) !no. of nonzero blocks
      call set_infoa(dsp_data%INFOA,'n',nnz,ierr)
      call set_infoa(dsp_data%INFOA,'d',lb,ierr) !row-dim of a block
      call set_infoa(dsp_data%INFOA,'e',lb,ierr) !col-dim of a block
      call set_infoa(dsp_data%INFOA,'f',mb,ierr) !row-dim in blocks
      call set_infoa(dsp_data%INFOA,'g',kb,ierr) !col-dim in blocks
      if((mb.ne.size(bpntrb)).or.(m.ne.mb*lb).or.(n.ne.kb*lb).or.&
        (minval(bindx).lt.base).or.(maxval(bindx).gt.kb-1+base).or.&
        (mb.ne.size(bpntre)).or.(nnz*lb*lb.ne.size(val))) then
        call usds(nmb,ierr)
        ierr = blas_error_param
        return
      end if
      if (options.gt.0) then
         ! decision rule whether or not to copy 
         COPY = .TRUE.
         if(COPY) then 
            options = -1 !copy
         else
            options = 0  !reference
         end if
      end if
      if (options.eq.0) then
         call set_infoa(dsp_data%INFOA,'c',REF_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! create reference to original matrix
         dsp_data%A => val
         dsp_data%IA1 => bindx
         dsp_data%pb => bpntrb
         dsp_data%pe => bpntre
         istat = 0
      else
         ! The additional required memory is DEALLOCATED later in USDS!
         call set_infoa(dsp_data%INFOA,'c',COP_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! copy original data
         allocate(dsp_data%A(size(val)),dsp_data%IA1(size(bindx)),&
                  dsp_data%pb(size(bpntrb)),dsp_data%pe(size(bpntre))&
                  ,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = bindx
         dsp_data%pb = bpntrb
         dsp_data%pe = bpntre
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine suscr_bsr 
! **********************************************************************
! **********************************************************************
      subroutine duscr_bsr (m,n,val,bindx,bpntrb,bpntre,mb,kb,lb,&
                           prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,mb,kb,lb,prpty
      integer, dimension(:), intent(inout),target :: bindx,bpntrb,bpntre
      real(KIND=dp) , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr, options, base , nnz
      logical :: COPY
      character :: message
      type(dspmat ),pointer :: dsp_data
      options = istat
      istat = -1 !if not changed later, routine has failed
      a = 0 
      call new_dsp (nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      call accessdata(dsp_data,nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      dsp_data%FIDA = 'BSR'
      dsp_data%M = m
      dsp_data%K = n
      call set_descra(dsp_data%DESCRA,prpty,ierr)
      call get_descra(dsp_data%DESCRA,'b',message,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if      
      if (message.eq.'C') then
         base = C_BASE
      else !Assuming F base
         base = F_BASE
      end if
      call set_infoa(dsp_data%INFOA,'b',base,ierr)
      nnz = size(bindx) !no. of nonzero blocks
      call set_infoa(dsp_data%INFOA,'n',nnz,ierr)
      call set_infoa(dsp_data%INFOA,'d',lb,ierr) !row-dim of a block
      call set_infoa(dsp_data%INFOA,'e',lb,ierr) !col-dim of a block
      call set_infoa(dsp_data%INFOA,'f',mb,ierr) !row-dim in blocks
      call set_infoa(dsp_data%INFOA,'g',kb,ierr) !col-dim in blocks
      if((mb.ne.size(bpntrb)).or.(m.ne.mb*lb).or.(n.ne.kb*lb).or.&
        (minval(bindx).lt.base).or.(maxval(bindx).gt.kb-1+base).or.&
        (mb.ne.size(bpntre)).or.(nnz*lb*lb.ne.size(val))) then
        call usds(nmb,ierr)
        ierr = blas_error_param
        return
      end if
      if (options.gt.0) then
         ! decision rule whether or not to copy 
         COPY = .TRUE.
         if(COPY) then 
            options = -1 !copy
         else
            options = 0  !reference
         end if
      end if
      if (options.eq.0) then
         call set_infoa(dsp_data%INFOA,'c',REF_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! create reference to original matrix
         dsp_data%A => val
         dsp_data%IA1 => bindx
         dsp_data%pb => bpntrb
         dsp_data%pe => bpntre
         istat = 0
      else
         ! The additional required memory is DEALLOCATED later in USDS!
         call set_infoa(dsp_data%INFOA,'c',COP_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! copy original data
         allocate(dsp_data%A(size(val)),dsp_data%IA1(size(bindx)),&
                  dsp_data%pb(size(bpntrb)),dsp_data%pe(size(bpntre))&
                  ,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = bindx
         dsp_data%pb = bpntrb
         dsp_data%pe = bpntre
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine duscr_bsr 
! **********************************************************************
! **********************************************************************
      subroutine cuscr_bsr (m,n,val,bindx,bpntrb,bpntre,mb,kb,lb,&
                           prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,mb,kb,lb,prpty
      integer, dimension(:), intent(inout),target :: bindx,bpntrb,bpntre
      complex(KIND=sp) , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr, options, base , nnz
      logical :: COPY
      character :: message
      type(cspmat ),pointer :: dsp_data
      options = istat
      istat = -1 !if not changed later, routine has failed
      a = 0 
      call new_csp (nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      call accessdata(dsp_data,nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      dsp_data%FIDA = 'BSR'
      dsp_data%M = m
      dsp_data%K = n
      call set_descra(dsp_data%DESCRA,prpty,ierr)
      call get_descra(dsp_data%DESCRA,'b',message,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if      
      if (message.eq.'C') then
         base = C_BASE
      else !Assuming F base
         base = F_BASE
      end if
      call set_infoa(dsp_data%INFOA,'b',base,ierr)
      nnz = size(bindx) !no. of nonzero blocks
      call set_infoa(dsp_data%INFOA,'n',nnz,ierr)
      call set_infoa(dsp_data%INFOA,'d',lb,ierr) !row-dim of a block
      call set_infoa(dsp_data%INFOA,'e',lb,ierr) !col-dim of a block
      call set_infoa(dsp_data%INFOA,'f',mb,ierr) !row-dim in blocks
      call set_infoa(dsp_data%INFOA,'g',kb,ierr) !col-dim in blocks
      if((mb.ne.size(bpntrb)).or.(m.ne.mb*lb).or.(n.ne.kb*lb).or.&
        (minval(bindx).lt.base).or.(maxval(bindx).gt.kb-1+base).or.&
        (mb.ne.size(bpntre)).or.(nnz*lb*lb.ne.size(val))) then
        call usds(nmb,ierr)
        ierr = blas_error_param
        return
      end if
      if (options.gt.0) then
         ! decision rule whether or not to copy 
         COPY = .TRUE.
         if(COPY) then 
            options = -1 !copy
         else
            options = 0  !reference
         end if
      end if
      if (options.eq.0) then
         call set_infoa(dsp_data%INFOA,'c',REF_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! create reference to original matrix
         dsp_data%A => val
         dsp_data%IA1 => bindx
         dsp_data%pb => bpntrb
         dsp_data%pe => bpntre
         istat = 0
      else
         ! The additional required memory is DEALLOCATED later in USDS!
         call set_infoa(dsp_data%INFOA,'c',COP_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! copy original data
         allocate(dsp_data%A(size(val)),dsp_data%IA1(size(bindx)),&
                  dsp_data%pb(size(bpntrb)),dsp_data%pe(size(bpntre))&
                  ,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = bindx
         dsp_data%pb = bpntrb
         dsp_data%pe = bpntre
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine cuscr_bsr 
! **********************************************************************
! **********************************************************************
      subroutine zuscr_bsr (m,n,val,bindx,bpntrb,bpntre,mb,kb,lb,&
                           prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,mb,kb,lb,prpty
      integer, dimension(:), intent(inout),target :: bindx,bpntrb,bpntre
      complex(KIND=dp) , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr, options, base , nnz
      logical :: COPY
      character :: message
      type(zspmat ),pointer :: dsp_data
      options = istat
      istat = -1 !if not changed later, routine has failed
      a = 0 
      call new_zsp (nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      call accessdata(dsp_data,nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if
      dsp_data%FIDA = 'BSR'
      dsp_data%M = m
      dsp_data%K = n
      call set_descra(dsp_data%DESCRA,prpty,ierr)
      call get_descra(dsp_data%DESCRA,'b',message,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if      
      if (message.eq.'C') then
         base = C_BASE
      else !Assuming F base
         base = F_BASE
      end if
      call set_infoa(dsp_data%INFOA,'b',base,ierr)
      nnz = size(bindx) !no. of nonzero blocks
      call set_infoa(dsp_data%INFOA,'n',nnz,ierr)
      call set_infoa(dsp_data%INFOA,'d',lb,ierr) !row-dim of a block
      call set_infoa(dsp_data%INFOA,'e',lb,ierr) !col-dim of a block
      call set_infoa(dsp_data%INFOA,'f',mb,ierr) !row-dim in blocks
      call set_infoa(dsp_data%INFOA,'g',kb,ierr) !col-dim in blocks
      if((mb.ne.size(bpntrb)).or.(m.ne.mb*lb).or.(n.ne.kb*lb).or.&
        (minval(bindx).lt.base).or.(maxval(bindx).gt.kb-1+base).or.&
        (mb.ne.size(bpntre)).or.(nnz*lb*lb.ne.size(val))) then
        call usds(nmb,ierr)
        ierr = blas_error_param
        return
      end if
      if (options.gt.0) then
         ! decision rule whether or not to copy 
         COPY = .TRUE.
         if(COPY) then 
            options = -1 !copy
         else
            options = 0  !reference
         end if
      end if
      if (options.eq.0) then
         call set_infoa(dsp_data%INFOA,'c',REF_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! create reference to original matrix
         dsp_data%A => val
         dsp_data%IA1 => bindx
         dsp_data%pb => bpntrb
         dsp_data%pe => bpntre
         istat = 0
      else
         ! The additional required memory is DEALLOCATED later in USDS!
         call set_infoa(dsp_data%INFOA,'c',COP_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! copy original data
         allocate(dsp_data%A(size(val)),dsp_data%IA1(size(bindx)),&
                  dsp_data%pb(size(bpntrb)),dsp_data%pe(size(bpntre))&
                  ,STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = bindx
         dsp_data%pb = bpntrb
         dsp_data%pe = bpntre
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine zuscr_bsr 
! **********************************************************************
! **********************************************************************
      end module mod_uscr_bsr

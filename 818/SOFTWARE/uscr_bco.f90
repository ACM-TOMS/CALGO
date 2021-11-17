module mod_uscr_bco
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : CREATION ROUTINE FOR MATRIX HANDLE FROM 'BCO'-FORMAT
! **********************************************************************
use representation_of_data
use properties
use mod_usds
implicit none
interface uscr_bco
  module procedure iuscr_bco
  module procedure suscr_bco
  module procedure duscr_bco
  module procedure cuscr_bco
  module procedure zuscr_bco
end interface
contains
! **********************************************************************
! **********************************************************************
      subroutine iuscr_bco (m,n,val,bindx,bjndx,bnnz,mb,kb,lb, &
                           prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,bnnz,mb,kb,lb,prpty
      integer, dimension(:), intent(inout),target :: bindx,bjndx
      integer , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr, options, base
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
      dsp_data%FIDA = 'BCO'
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
      call set_infoa(dsp_data%INFOA,'n',bnnz,ierr)
      call set_infoa(dsp_data%INFOA,'d',lb,ierr) !row-dim of a block
      call set_infoa(dsp_data%INFOA,'e',lb,ierr) !col-dim of a block
      call set_infoa(dsp_data%INFOA,'f',mb,ierr) !row-dim in blocks
      call set_infoa(dsp_data%INFOA,'g',kb,ierr) !col-dim in blocks
      if((bnnz.ne.size(bindx)).or.(bnnz.ne.size(bjndx)).or.&
        (m.ne.mb*lb).or.(n.ne.kb*lb).or.(bnnz*lb*lb.ne.size(val)).or.&
        (minval(bindx).lt.base).or.(maxval(bindx).gt.mb-1+base).or.&
        (minval(bjndx).lt.base).or.(maxval(bjndx).gt.kb-1+base)) then
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
         dsp_data%IA2 => bjndx
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
                  dsp_data%IA2(size(bjndx)),STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = bindx
         dsp_data%IA2 = bjndx
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine iuscr_bco 
! **********************************************************************
! **********************************************************************
      subroutine suscr_bco (m,n,val,bindx,bjndx,bnnz,mb,kb,lb, &
                           prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,bnnz,mb,kb,lb,prpty
      integer, dimension(:), intent(inout),target :: bindx,bjndx
      real(KIND=sp) , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr, options, base
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
      dsp_data%FIDA = 'BCO'
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
      call set_infoa(dsp_data%INFOA,'n',bnnz,ierr)
      call set_infoa(dsp_data%INFOA,'d',lb,ierr) !row-dim of a block
      call set_infoa(dsp_data%INFOA,'e',lb,ierr) !col-dim of a block
      call set_infoa(dsp_data%INFOA,'f',mb,ierr) !row-dim in blocks
      call set_infoa(dsp_data%INFOA,'g',kb,ierr) !col-dim in blocks
      if((bnnz.ne.size(bindx)).or.(bnnz.ne.size(bjndx)).or.&
        (m.ne.mb*lb).or.(n.ne.kb*lb).or.(bnnz*lb*lb.ne.size(val)).or.&
        (minval(bindx).lt.base).or.(maxval(bindx).gt.mb-1+base).or.&
        (minval(bjndx).lt.base).or.(maxval(bjndx).gt.kb-1+base)) then
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
         dsp_data%IA2 => bjndx
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
                  dsp_data%IA2(size(bjndx)),STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = bindx
         dsp_data%IA2 = bjndx
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine suscr_bco 
! **********************************************************************
! **********************************************************************
      subroutine duscr_bco (m,n,val,bindx,bjndx,bnnz,mb,kb,lb, &
                           prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,bnnz,mb,kb,lb,prpty
      integer, dimension(:), intent(inout),target :: bindx,bjndx
      real(KIND=dp) , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr, options, base
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
      dsp_data%FIDA = 'BCO'
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
      call set_infoa(dsp_data%INFOA,'n',bnnz,ierr)
      call set_infoa(dsp_data%INFOA,'d',lb,ierr) !row-dim of a block
      call set_infoa(dsp_data%INFOA,'e',lb,ierr) !col-dim of a block
      call set_infoa(dsp_data%INFOA,'f',mb,ierr) !row-dim in blocks
      call set_infoa(dsp_data%INFOA,'g',kb,ierr) !col-dim in blocks
      if((bnnz.ne.size(bindx)).or.(bnnz.ne.size(bjndx)).or.&
        (m.ne.mb*lb).or.(n.ne.kb*lb).or.(bnnz*lb*lb.ne.size(val)).or.&
        (minval(bindx).lt.base).or.(maxval(bindx).gt.mb-1+base).or.&
        (minval(bjndx).lt.base).or.(maxval(bjndx).gt.kb-1+base)) then
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
         dsp_data%IA2 => bjndx
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
                  dsp_data%IA2(size(bjndx)),STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = bindx
         dsp_data%IA2 = bjndx
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine duscr_bco 
! **********************************************************************
! **********************************************************************
      subroutine cuscr_bco (m,n,val,bindx,bjndx,bnnz,mb,kb,lb, &
                           prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,bnnz,mb,kb,lb,prpty
      integer, dimension(:), intent(inout),target :: bindx,bjndx
      complex(KIND=sp) , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr, options, base
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
      dsp_data%FIDA = 'BCO'
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
      call set_infoa(dsp_data%INFOA,'n',bnnz,ierr)
      call set_infoa(dsp_data%INFOA,'d',lb,ierr) !row-dim of a block
      call set_infoa(dsp_data%INFOA,'e',lb,ierr) !col-dim of a block
      call set_infoa(dsp_data%INFOA,'f',mb,ierr) !row-dim in blocks
      call set_infoa(dsp_data%INFOA,'g',kb,ierr) !col-dim in blocks
      if((bnnz.ne.size(bindx)).or.(bnnz.ne.size(bjndx)).or.&
        (m.ne.mb*lb).or.(n.ne.kb*lb).or.(bnnz*lb*lb.ne.size(val)).or.&
        (minval(bindx).lt.base).or.(maxval(bindx).gt.mb-1+base).or.&
        (minval(bjndx).lt.base).or.(maxval(bjndx).gt.kb-1+base)) then
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
         dsp_data%IA2 => bjndx
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
                  dsp_data%IA2(size(bjndx)),STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = bindx
         dsp_data%IA2 = bjndx
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine cuscr_bco 
! **********************************************************************
! **********************************************************************
      subroutine zuscr_bco (m,n,val,bindx,bjndx,bnnz,mb,kb,lb, &
                           prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,bnnz,mb,kb,lb,prpty
      integer, dimension(:), intent(inout),target :: bindx,bjndx
      complex(KIND=dp) , dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr, options, base
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
      dsp_data%FIDA = 'BCO'
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
      call set_infoa(dsp_data%INFOA,'n',bnnz,ierr)
      call set_infoa(dsp_data%INFOA,'d',lb,ierr) !row-dim of a block
      call set_infoa(dsp_data%INFOA,'e',lb,ierr) !col-dim of a block
      call set_infoa(dsp_data%INFOA,'f',mb,ierr) !row-dim in blocks
      call set_infoa(dsp_data%INFOA,'g',kb,ierr) !col-dim in blocks
      if((bnnz.ne.size(bindx)).or.(bnnz.ne.size(bjndx)).or.&
        (m.ne.mb*lb).or.(n.ne.kb*lb).or.(bnnz*lb*lb.ne.size(val)).or.&
        (minval(bindx).lt.base).or.(maxval(bindx).gt.mb-1+base).or.&
        (minval(bjndx).lt.base).or.(maxval(bjndx).gt.kb-1+base)) then
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
         dsp_data%IA2 => bjndx
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
                  dsp_data%IA2(size(bjndx)),STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = bindx
         dsp_data%IA2 = bjndx
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine zuscr_bco 
! **********************************************************************
! **********************************************************************
end module mod_uscr_bco

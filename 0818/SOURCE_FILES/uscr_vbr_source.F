      subroutine XUSCR_VBR(m,n,val,indx,bindx,rpntr,cpntr,bpntrb,&
                           bpntre,mb,kb,prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,mb,kb,prpty
      integer, dimension(:), intent(inout),target :: indx,bindx,&
                             rpntr,cpntr,bpntrb,bpntre
      DCOMPLEX, dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr, options, base , nnz
      logical :: COPY
      character :: message
      type(XSPMAT),pointer :: dsp_data

      options = istat
      istat = -1 !if not changed later, routine has failed
      a = 0 
      call NEW_XSP(nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_memalloc
         return
      end if
      call accessdata(dsp_data,nmb,ierr)
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if

      dsp_data%FIDA = 'VBR'
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
      call set_infoa(dsp_data%INFOA,'d',-1,ierr) !row-dim of block NOT fixed
      call set_infoa(dsp_data%INFOA,'e',-1,ierr) !col-dim of block NOT fixed
      call set_infoa(dsp_data%INFOA,'f',mb,ierr) !row-dim in blocks
      call set_infoa(dsp_data%INFOA,'g',kb,ierr) !col-dim in blocks

      if((mb.ne.size(bpntrb)).or.(mb.ne.size(bpntre)).or.&
        (size(val).ne.maxval(indx)-base).or.(minval(indx).lt.base).or.&
        (minval(bindx).lt.base).or.(maxval(bindx).gt.kb-1+base)) then
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
         dsp_data%IA2 => indx
         dsp_data%pb => bpntrb
         dsp_data%pe => bpntre
         dsp_data%bp1 => rpntr
         dsp_data%bp2 => cpntr
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
                  dsp_data%IA2(size(indx)),&
                  dsp_data%pb(size(bpntrb)),dsp_data%pe(size(bpntre)),&
                  dsp_data%bp1(size(rpntr)),dsp_data%bp2(size(cpntr)),&
                  STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = bindx
         dsp_data%IA2 = indx
         dsp_data%pb = bpntrb
         dsp_data%pe = bpntre
         dsp_data%bp1 = rpntr
         dsp_data%bp2 = cpntr
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine XUSCR_VBR

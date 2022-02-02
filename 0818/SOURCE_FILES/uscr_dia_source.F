      subroutine XUSCR_DIA(m,n,val,lda,idiag,ndiag,prpty,istat,a)
      implicit none
      integer, intent(in) :: m,n,lda,ndiag,prpty
      integer, dimension(:), intent(inout),target :: idiag
      DCOMPLEX, dimension(:), intent(inout),target :: val
      integer, intent(inout) :: istat
      integer, intent(out) :: a
      integer :: nmb,ierr,nnz,options,base
      logical :: COPY
      character :: message
      type(XSPMAT),pointer :: dsp_data

      options = istat
      istat = -1 !if not changed later, routine has failed

      a = 0 

      if((ndiag.ne.size(idiag)).or.(lda*ndiag.ne.size(val)).or.&
        (maxval(idiag).gt.n).or.(minval(idiag).lt.-m).or.&
        (lda.ne.min(m,n))) then
        ierr = blas_error_param
        return
      end if

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

      dsp_data%FIDA = 'DIA'
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
      if (ierr.ne.0) then
         ierr = blas_error_param
         return
      end if      
      nnz = count(val.ne.0.)
      if(nnz.le.lda*ndiag*0.5) then
         ! Warning Many zeros stored
      end if
      call set_infoa(dsp_data%INFOA,'n',nnz,ierr)
      call set_infoa(dsp_data%INFOA,'d',lda,ierr) !row-dim of val
      call set_infoa(dsp_data%INFOA,'e',ndiag,ierr) !col-dim of val

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
         dsp_data%IA1 => idiag
         istat = 0
      else
         ! The additional required memory is DEALLOCATED later in USDS!
         call set_infoa(dsp_data%INFOA,'c',COP_OF_SOURCE,ierr)
         if (ierr.ne.0) then
            ierr = blas_error_param
            return
         end if
         ! copy original data
         allocate(dsp_data%A(size(val)),dsp_data%IA1(size(idiag)),&
                             STAT=ierr)
         if(ierr.ne.0) then
            ierr = blas_error_memalloc
            return
         end if
         dsp_data%A = val
         dsp_data%IA1 = idiag
         istat = 1
      end if
      if(istat.ge.0)  a = nmb
      end subroutine XUSCR_DIA

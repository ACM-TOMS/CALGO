      module mod_info
! **********************************************************************
!     Author : C. Voemel
!     Date of last modification : 7.7.00
!     Description : FOR DEBUGGING ONLY !!! 
!                   "print" displays data for given handle number
! **********************************************************************
      use representation_of_data
      use properties
      implicit none
      contains
      subroutine print(nmb,ierr)
      implicit none
      intrinsic modulo
      integer, intent(in) :: nmb
      integer, intent(out) :: ierr
      type(ispmat),pointer :: isp_data
      type(sspmat),pointer :: ssp_data
      type(dspmat),pointer :: dsp_data
      type(cspmat),pointer :: csp_data
      type(zspmat),pointer :: zsp_data
      integer :: rest,base,copy,nnz,rowdim,coldim
      character :: style,diag,type,part,store
      rest = modulo(nmb,no_of_types)
      select case(rest)
      case(ISP_MATRIX)
! **********************************************************************
         call accessdata(isp_data ,nmb,ierr)
         if (ierr.ne.0) then
            write(*,*) '***********************************'
            write(*,*) 'No data for no. ',nmb,' available !' 
            write(*,*) '***********************************'
            return
         end if
         write(*,*) '***********************************'
         write(*,*) 'Matrix no. ', nmb
         write(*,*) 'number of rows : ', isp_data %M
         write(*,*) 'number of columns : ', isp_data %K
         write(*,*) 'Storage : ', isp_data %FIDA
         write(*,*) 'A : ', isp_data %A
         write(*,*) 'IA1 : ', isp_data %IA1
         write(*,*) 'IA2 : ', isp_data %IA2
         write(*,*) '***********************************'
         call get_descra(isp_data %DESCRA,'a',part,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Matrix part accessed : ',part
         end if
         call get_descra(isp_data %DESCRA,'b',style,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Index style : ',style
         end if
         call get_descra(isp_data %DESCRA,'d',diag,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Unity-diagonal : ',diag
         end if
         call get_descra(isp_data %DESCRA,'f',store,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Block-internal storage : ',store
         end if
         call get_descra(isp_data %DESCRA,'t',type,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Matrix type : ',type
         end if
         call get_infoa(isp_data %INFOA,'b',base,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Indices start at : ',base 
         end if
         call get_infoa(isp_data %INFOA,'c',copy,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Matrix is copy of original data : ',copy
         end if
         call get_infoa(isp_data %INFOA,'n',nnz,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'number of non-zero(-block)s : ',nnz
         end if
         call get_infoa(isp_data %INFOA,'d',rowdim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) '(Multi-dim arrays) row dim of block : ',rowdim
         end if
         call get_infoa(isp_data %INFOA,'e',coldim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) '(Multi-dim arrays) col dim of block : ',coldim
         end if
         call get_infoa(isp_data %INFOA,'f',rowdim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Block structure : row-dim in blocks : ',rowdim
         end if
         call get_infoa(isp_data %INFOA,'g',coldim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Block structure : col-dim in blocks : ',coldim
         end if
         write(*,*) '***********************************'
! **********************************************************************
      case(SSP_MATRIX)
! **********************************************************************
         call accessdata(ssp_data ,nmb,ierr)
         if (ierr.ne.0) then
            write(*,*) '***********************************'
            write(*,*) 'No data for no. ',nmb,' available !' 
            write(*,*) '***********************************'
            return
         end if
         write(*,*) '***********************************'
         write(*,*) 'Matrix no. ', nmb
         write(*,*) 'number of rows : ', ssp_data %M
         write(*,*) 'number of columns : ', ssp_data %K
         write(*,*) 'Storage : ', ssp_data %FIDA
         write(*,*) 'A : ', ssp_data %A
         write(*,*) 'IA1 : ', ssp_data %IA1
         write(*,*) 'IA2 : ', ssp_data %IA2
         write(*,*) '***********************************'
         call get_descra(ssp_data %DESCRA,'a',part,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Matrix part accessed : ',part
         end if
         call get_descra(ssp_data %DESCRA,'b',style,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Index style : ',style
         end if
         call get_descra(ssp_data %DESCRA,'d',diag,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Unity-diagonal : ',diag
         end if
         call get_descra(ssp_data %DESCRA,'f',store,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Block-internal storage : ',store
         end if
         call get_descra(ssp_data %DESCRA,'t',type,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Matrix type : ',type
         end if
         call get_infoa(ssp_data %INFOA,'b',base,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Indices start at : ',base 
         end if
         call get_infoa(ssp_data %INFOA,'c',copy,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Matrix is copy of original data : ',copy
         end if
         call get_infoa(ssp_data %INFOA,'n',nnz,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'number of non-zero(-block)s : ',nnz
         end if
         call get_infoa(ssp_data %INFOA,'d',rowdim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) '(Multi-dim arrays) row dim of block : ',rowdim
         end if
         call get_infoa(ssp_data %INFOA,'e',coldim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) '(Multi-dim arrays) col dim of block : ',coldim
         end if
         call get_infoa(ssp_data %INFOA,'f',rowdim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Block structure : row-dim in blocks : ',rowdim
         end if
         call get_infoa(ssp_data %INFOA,'g',coldim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Block structure : col-dim in blocks : ',coldim
         end if
         write(*,*) '***********************************'
! **********************************************************************
      case(DSP_MATRIX)
! **********************************************************************
         call accessdata(dsp_data ,nmb,ierr)
         if (ierr.ne.0) then
            write(*,*) '***********************************'
            write(*,*) 'No data for no. ',nmb,' available !' 
            write(*,*) '***********************************'
            return
         end if
         write(*,*) '***********************************'
         write(*,*) 'Matrix no. ', nmb
         write(*,*) 'number of rows : ', dsp_data %M
         write(*,*) 'number of columns : ', dsp_data %K
         write(*,*) 'Storage : ', dsp_data %FIDA
         write(*,*) 'A : ', dsp_data %A
         write(*,*) 'IA1 : ', dsp_data %IA1
         write(*,*) 'IA2 : ', dsp_data %IA2
         write(*,*) '***********************************'
         call get_descra(dsp_data %DESCRA,'a',part,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Matrix part accessed : ',part
         end if
         call get_descra(dsp_data %DESCRA,'b',style,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Index style : ',style
         end if
         call get_descra(dsp_data %DESCRA,'d',diag,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Unity-diagonal : ',diag
         end if
         call get_descra(dsp_data %DESCRA,'f',store,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Block-internal storage : ',store
         end if
         call get_descra(dsp_data %DESCRA,'t',type,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Matrix type : ',type
         end if
         call get_infoa(dsp_data %INFOA,'b',base,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Indices start at : ',base 
         end if
         call get_infoa(dsp_data %INFOA,'c',copy,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Matrix is copy of original data : ',copy
         end if
         call get_infoa(dsp_data %INFOA,'n',nnz,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'number of non-zero(-block)s : ',nnz
         end if
         call get_infoa(dsp_data %INFOA,'d',rowdim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) '(Multi-dim arrays) row dim of block : ',rowdim
         end if
         call get_infoa(dsp_data %INFOA,'e',coldim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) '(Multi-dim arrays) col dim of block : ',coldim
         end if
         call get_infoa(dsp_data %INFOA,'f',rowdim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Block structure : row-dim in blocks : ',rowdim
         end if
         call get_infoa(dsp_data %INFOA,'g',coldim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Block structure : col-dim in blocks : ',coldim
         end if
         write(*,*) '***********************************'
! **********************************************************************
      case(CSP_MATRIX)
! **********************************************************************
         call accessdata(csp_data ,nmb,ierr)
         if (ierr.ne.0) then
            write(*,*) '***********************************'
            write(*,*) 'No data for no. ',nmb,' available !' 
            write(*,*) '***********************************'
            return
         end if
         write(*,*) '***********************************'
         write(*,*) 'Matrix no. ', nmb
         write(*,*) 'number of rows : ', csp_data %M
         write(*,*) 'number of columns : ', csp_data %K
         write(*,*) 'Storage : ', csp_data %FIDA
         write(*,*) 'A : ', csp_data %A
         write(*,*) 'IA1 : ', csp_data %IA1
         write(*,*) 'IA2 : ', csp_data %IA2
         write(*,*) '***********************************'
         call get_descra(csp_data %DESCRA,'a',part,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Matrix part accessed : ',part
         end if
         call get_descra(csp_data %DESCRA,'b',style,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Index style : ',style
         end if
         call get_descra(csp_data %DESCRA,'d',diag,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Unity-diagonal : ',diag
         end if
         call get_descra(csp_data %DESCRA,'f',store,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Block-internal storage : ',store
         end if
         call get_descra(csp_data %DESCRA,'t',type,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Matrix type : ',type
         end if
         call get_infoa(csp_data %INFOA,'b',base,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Indices start at : ',base 
         end if
         call get_infoa(csp_data %INFOA,'c',copy,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Matrix is copy of original data : ',copy
         end if
         call get_infoa(csp_data %INFOA,'n',nnz,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'number of non-zero(-block)s : ',nnz
         end if
         call get_infoa(csp_data %INFOA,'d',rowdim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) '(Multi-dim arrays) row dim of block : ',rowdim
         end if
         call get_infoa(csp_data %INFOA,'e',coldim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) '(Multi-dim arrays) col dim of block : ',coldim
         end if
         call get_infoa(csp_data %INFOA,'f',rowdim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Block structure : row-dim in blocks : ',rowdim
         end if
         call get_infoa(csp_data %INFOA,'g',coldim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Block structure : col-dim in blocks : ',coldim
         end if
         write(*,*) '***********************************'
! **********************************************************************
      case(ZSP_MATRIX)
! **********************************************************************
         call accessdata(zsp_data ,nmb,ierr)
         if (ierr.ne.0) then
            write(*,*) '***********************************'
            write(*,*) 'No data for no. ',nmb,' available !' 
            write(*,*) '***********************************'
            return
         end if
         write(*,*) '***********************************'
         write(*,*) 'Matrix no. ', nmb
         write(*,*) 'number of rows : ', zsp_data %M
         write(*,*) 'number of columns : ', zsp_data %K
         write(*,*) 'Storage : ', zsp_data %FIDA
         write(*,*) 'A : ', zsp_data %A
         write(*,*) 'IA1 : ', zsp_data %IA1
         write(*,*) 'IA2 : ', zsp_data %IA2
         write(*,*) '***********************************'
         call get_descra(zsp_data %DESCRA,'a',part,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Matrix part accessed : ',part
         end if
         call get_descra(zsp_data %DESCRA,'b',style,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Index style : ',style
         end if
         call get_descra(zsp_data %DESCRA,'d',diag,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Unity-diagonal : ',diag
         end if
         call get_descra(zsp_data %DESCRA,'f',store,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Block-internal storage : ',store
         end if
         call get_descra(zsp_data %DESCRA,'t',type,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Matrix type : ',type
         end if
         call get_infoa(zsp_data %INFOA,'b',base,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Indices start at : ',base 
         end if
         call get_infoa(zsp_data %INFOA,'c',copy,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Matrix is copy of original data : ',copy
         end if
         call get_infoa(zsp_data %INFOA,'n',nnz,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'number of non-zero(-block)s : ',nnz
         end if
         call get_infoa(zsp_data %INFOA,'d',rowdim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) '(Multi-dim arrays) row dim of block : ',rowdim
         end if
         call get_infoa(zsp_data %INFOA,'e',coldim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) '(Multi-dim arrays) col dim of block : ',coldim
         end if
         call get_infoa(zsp_data %INFOA,'f',rowdim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Block structure : row-dim in blocks : ',rowdim
         end if
         call get_infoa(zsp_data %INFOA,'g',coldim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Block structure : col-dim in blocks : ',coldim
         end if
         write(*,*) '***********************************'
! **********************************************************************
      case default
         write(*,*) 'Wrong matrix type !'
         ierr = -1 
      end select
      end subroutine print
      end module mod_info

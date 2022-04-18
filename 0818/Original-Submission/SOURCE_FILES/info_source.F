         call accessdata(XSP_DATA,nmb,ierr)
         if (ierr.ne.0) then
            write(*,*) '***********************************'
            write(*,*) 'No data for no. ',nmb,' available !' 
            write(*,*) '***********************************'
            return
         end if
         write(*,*) '***********************************'
         write(*,*) 'Matrix no. ', nmb
         write(*,*) 'number of rows : ', XSP_DATA%M
         write(*,*) 'number of columns : ', XSP_DATA%K
         write(*,*) 'Storage : ', XSP_DATA%FIDA
         write(*,*) 'A : ', XSP_DATA%A
         write(*,*) 'IA1 : ', XSP_DATA%IA1
         write(*,*) 'IA2 : ', XSP_DATA%IA2
         write(*,*) '***********************************'
         call get_descra(XSP_DATA%DESCRA,'a',part,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Matrix part accessed : ',part
         end if
         call get_descra(XSP_DATA%DESCRA,'b',style,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Index style : ',style
         end if
         call get_descra(XSP_DATA%DESCRA,'d',diag,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Unity-diagonal : ',diag
         end if
         call get_descra(XSP_DATA%DESCRA,'f',store,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Block-internal storage : ',store
         end if
         call get_descra(XSP_DATA%DESCRA,'t',type,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Matrix type : ',type
         end if
         call get_infoa(XSP_DATA%INFOA,'b',base,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Indices start at : ',base 
         end if
         call get_infoa(XSP_DATA%INFOA,'c',copy,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Matrix is copy of original data : ',copy
         end if
         call get_infoa(XSP_DATA%INFOA,'n',nnz,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'number of non-zero(-block)s : ',nnz
         end if
         call get_infoa(XSP_DATA%INFOA,'d',rowdim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) '(Multi-dim arrays) row dim of block : ',rowdim
         end if
         call get_infoa(XSP_DATA%INFOA,'e',coldim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) '(Multi-dim arrays) col dim of block : ',coldim
         end if
         call get_infoa(XSP_DATA%INFOA,'f',rowdim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Block structure : row-dim in blocks : ',rowdim
         end if
         call get_infoa(XSP_DATA%INFOA,'g',coldim,ierr)
         if (ierr.ne.0) then
            write(*,*) 'No information avail. for that argument'
            return
         else
            write(*,*) 'Block structure : col-dim in blocks : ',coldim
         end if
         write(*,*) '***********************************'

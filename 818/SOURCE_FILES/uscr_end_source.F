
      istat=-1
      call access_matrix(XPMATRIX,a,istat)  
      if(istat.ne.0) return
      prpty=XPMATRIX%property  
      select case(XPMATRIX%format)
      case('block')
         call XUSCR_blockend(a,prpty,istat)
         if(istat.ne.0) return
      case('vblock')
         call XUSCR_varend (a,prpty,istat)
         if(istat.ne.0) return
      case('normal')
         call XUSCR_normend(a,prpty,istat)       
         if(istat.ne.0) return
      case default
         istat = blas_error_param
         return
      end select 
  

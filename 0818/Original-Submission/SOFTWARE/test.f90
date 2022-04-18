      program test
! **********************************************************************
!     Author : C. Voemel
!
!     Date of last modification : 7.7.00
!      
!     Description : SAMPLE ROUTINE FOR THE CERFACS TECHNICAL REPORT
! **********************************************************************

      ! include the SparseBLAS module
      use SparseBLAS 

      ! matrix handle and other parameters
      integer :: a, istat, prpty

      ! dense vectors
      double precision, dimension(5) :: x,y

      ! matrix data in COO format
      double precision, parameter, dimension(14) :: A_val_coo = &
       (/ 11,51,31,32,34,52,13,23,33,14,24,42,55,44 /)
      integer, parameter, dimension(14) :: A_indx_coo = &
       (/ 1,5,3,3,3,5,1,2,3,1,2,4,5,4 /)
      integer, parameter, dimension(14) :: A_jndx_coo = &
       (/ 1,1,1,2,4,2,3,3,3,4,4,2,5,4 /)
      prpty = blas_general + blas_one_base
      istat = -1 ! copy data when creating handle

      ! create matrix handle
      call uscr_coo(5,5,A_val_coo,A_indx_coo,A_jndx_coo,14, &
                          prpty,istat,a)
      if(istat.lt.0) then
         stop 'Error! No handle created!'
      end if

      ! calculate matrix-vector product y=A*x
      x = 1.0d0
      call usmv(a,x,y,istat)
      if (istat.ne.0) then
         stop 'Error! Can''t perform MV multiplication!'
      end if

      ! release matrix handle
      call usds(a, istat)
      if (istat.ne.0) then
         stop 'Error! Handle not released!'
      end if

      end program

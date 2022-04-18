program tester
use blas_sparse
use mod_test_parameters
real(kind=dp), dimension(:), allocatable :: x,y,z
real(kind=dp), dimension(:,:), allocatable :: dense_C, dense_B
complex(kind=dp), dimension(:), allocatable :: xc,yc,zc
complex(kind=dp), dimension(:,:), allocatable :: densec_C, densec_B
integer, dimension(:), allocatable :: x_i,y_i,z_i

integer :: i, prpty, a,d,ierr,b,c,dt,istat
integer :: n, nnz
integer,parameter::conj=1
integer :: USE_MATRIX = T_MATRIX
integer, parameter :: INDEX_BASE = blas_one_base
integer, parameter :: BLOCK_INTERN = blas_col_major 


 n = 5
 nnz = 14
  ! **********************************************************************     
       print *,'*********************************************************'
       if(conj.eq.0) then
          d=usdot(x_dot,indx_dot,y_dot)
       else
          d=usdot(x_dot,indx_dot,y_dot,conj)
       end if
       write(*,*) 'Testing USDOT'
       write(*,*) 'Errors of computed solution: ' 
       write(*,*) 'Error : ',abs(d - res_usdot)
       print *,'*********************************************************'
       if(alpha.eq.0) then
          call usaxpy(x_axpy,indx_axpy,y_axpy)
       else
          call usaxpy(x_axpy,indx_axpy,y_axpy,alpha)
       end if
       write(*,*) 'Testing USAXPY'
       write(*,*) 'Errors of computed solution: ' 
       write(*,*) 'Error : ',abs(y_axpy - res_usaxpy)
       print *,'*********************************************************'
       call usga(y_ga,x_ga,indx_ga)
       write(*,*) 'Testing USGA'
       write(*,*) 'Errors of computed solution: ' 
       write(*,*) 'Error : ',abs(x_ga - res_usga)
       print *,'*********************************************************'
       call usgz(y_gz,x_gz,indx_gz)
       write(*,*) 'Testing USGZ'
       write(*,*) 'Errors of computed solution: ' 
       write(*,*) 'Error : ',abs(x_gz - res_usgz_x)
       write(*,*) 'Error : ',abs(y_gz - res_usgz_y)
       print *,'*********************************************************'
       call ussc(x_sc,y_sc,indx_sc)
       write(*,*) 'Testing USSC'
       write(*,*) 'Errors of computed solution: ' 
       write(*,*) 'Error : ',abs(y_sc - res_ussc)
   
  ! **********************************************************************     
        print *,'*********************************************************'
         print *,'TEST WITH MATRIX A'
         prpty = blas_general + INDEX_BASE
         allocate(x(5),y(5))
         y=0.
         call duscr_begin(A_m,A_n,a,istat)
         call ussp(a,prpty,istat)
         do i=1,size(A_val)
            call uscr_insert_entry(a,A_val(i),A_indx(i),A_jndx(i),istat)
         end do
         call uscr_end(a,istat )
         do i=1,size(x)
            x(i) = dble(i)
         end do
         allocate(z(size(x)))
         z = x
         allocate(dense_C(size(y),3),dense_B(size(x),3))
         do i = 1,3
            dense_B(:,i) = x
            dense_C(:,i) = 0.
         end do
         write(*,*) '* Test of MV multiplication *'
         call usmv(a,x,y,istat)
         write(*,*) 'Error : ',abs(y-Ax)
         write(*,*) '* Test of MM multiplication *'
         call usmm(a,dense_B,dense_C,istat)
         write(*,*) 'Error : ',(abs(dense_C(:,i)-Ax),i=1,3)
         write(*,*) '* Deleting matrix handle    *'
         deallocate(x,y,dense_C,dense_B,z)
         call usds(a, istat)
         if (istat.ne.0) then
            write(*,*) 'Deallocation failure'
            stop
         end if


        print *,'*********************************************************'
         print *,'TEST WITH MATRIX A_SU'

         prpty = blas_symmetric + INDEX_BASE + blas_upper
         allocate(x(5),y(5))
         y=0.
         call duscr_begin(A_SU_m,A_SU_n,c,istat) 
         call ussp(c,prpty,istat)
         call uscr_insert_clique(c,A_SU,A_SU_indx,A_SU_jndx,istat)
         call uscr_end(c,istat )
        
         do i=1,size(x)
            x(i) = dble(i)
         end do
         allocate(z(size(x)))
         z = x
         allocate(dense_C(size(y),3),dense_B(size(x),3))
         do i = 1,3
            dense_B(:,i) = x
            dense_C(:,i) = 0.
         end do
         write(*,*) '* Test of MV multiplication *'
         call usmv(c,x,y,istat)
         write(*,*) 'Error : ',abs(y-A_SUx)
         write(*,*) '* Test of MM multiplication *'
         call usmm(c,dense_B,dense_C,istat)
         write(*,*) 'Error : ',(abs(dense_C(:,i)-A_SUx),i=1,3)
         write(*,*) '* Deleting matrix handle    *'
         deallocate(x,y,dense_C,dense_B,z)
         call usds(c, istat)
         if (istat.ne.0) then
            write(*,*) 'Deallocation failure'
            stop
         end if
        print *,'*********************************************************'
         print *,'TEST WITH MATRIX A_SU'
         prpty = blas_symmetric + INDEX_BASE + blas_lower
         allocate(x(5),y(5))
         y=0.
         call duscr_begin(A_SL_m,A_SL_n,a,istat)
         call ussp(a,prpty,istat)
         call  uscr_insert_col(a,1,A_SL_1,A_SL_1_INDX,istat)
         call  uscr_insert_col(a,2,A_SL_2,A_SL_2_INDX,istat)
         call  uscr_insert_col(a,3,A_SL_3,A_SL_3_INDX,istat)
         call  uscr_insert_col(a,4,A_SL_4,A_SL_4_INDX,istat)
         call  uscr_insert_col(a,5,A_SL_5,A_SL_5_INDX,istat)
         call uscr_end(a,istat )
         do i=1,size(x)
            x(i) = dble(i)
         end do
         allocate(z(size(x)))
         z = x
         allocate(dense_C(size(y),3),dense_B(size(x),3))
         do i = 1,3
            dense_B(:,i) = x
            dense_C(:,i) = 0.
         end do
         write(*,*) '* Test of MV multiplication *'
         call usmv(a,x,y,istat)
         write(*,*) 'Error : ',abs(y-A_SLx)
         write(*,*) '* Test of MM multiplication *'
         call usmm(a,dense_B,dense_C,istat)
         write(*,*) 'Error : ',(abs(dense_C(:,i)-A_SLx),i=1,3)
         write(*,*) '* Deleting matrix handle    *'
         deallocate(x,y,dense_C,dense_B,z)
         call usds(a, istat)
         if (istat.ne.0) then
            write(*,*) 'Deallocation failure'
            stop
         end if

  ! **********************************************************************     
        print *,'*********************************************************'
         print *,'TEST WITH MATRIX B'
         prpty = blas_general + INDEX_BASE
         allocate(x(5),y(3))
         y=0.
         call duscr_begin(B_m,B_n,b,istat)
         call ussp(b,prpty,istat)
         call uscr_insert_row(b,1,B_1,B_1_JNDX,istat)
         call uscr_insert_row(b,2,B_2,B_2_JNDX,istat)
         call uscr_insert_row(b,3,B_3,B_3_JNDX,istat)
         call uscr_end(b,istat )
         !call print(b,istat) 
         do i=1,size(x)
            x(i) = dble(i)
         end do
         allocate(z(size(x)))
         z = x
         allocate(dense_C(size(y),3),dense_B(size(x),3))
         do i = 1,3
            dense_B(:,i) = x
            dense_C(:,i) = 0.
         end do
         write(*,*) '* Test of MV multiplication *'
         call usmv(b,x,y,istat)
         write(*,*) 'Error : ',abs(y-Bx)
         write(*,*) '* Test of MM multiplication *'
         call usmm(b,dense_B,dense_C,istat)
         write(*,*) 'Error : ',(abs(dense_C(:,i)-Bx),i=1,3)
         write(*,*) '* Deleting matrix handle    *'
         deallocate(x,y,dense_C,dense_B,z)
         call usds(b, istat)
         if (istat.ne.0) then
            write(*,*) 'Deallocation failure'
            stop
         end if
  ! **********************************************************************     
        print *,'*********************************************************'
         print *,'TEST WITH MATRIX C'
         USE_MATRIX = O_MATRIX
         prpty = blas_general + INDEX_BASE
         allocate(xc(3),yc(3))    
         call zuscr_begin(C_m,C_n,a,istat)
         call ussp(a,prpty,istat)
         call  uscr_insert_entries(a,VAL_C,INDX_C,JNDX_C,istat)
         call uscr_end(a,istat )
         do i=1,size(xc)
            xc(i) = (1.,1.)
         end do
         allocate(zc(size(xc)))     
         zc = xc
         yc = (0.,0.)
         allocate(densec_C(size(yc),3),densec_B(size(xc),3))
  
         do i = 1,3
            densec_B(:,i) = xc
            densec_C(:,i) = (0.,0.)
         end do
         write(*,*) '* Test of MV multiplication *'
         select case(USE_MATRIX)
         case (O_MATRIX)
            call usmv(a,xc,yc,istat)
         case (T_MATRIX)
            call usmv(a,xc,yc,istat, transa = TRANSP_MATRIX)
         case (H_MATRIX)
            call usmv(a,xc,yc,istat, transa = HERMIT_MATRIX)
         case default
            stop 'No valid choice'
         end select
         if (istat.ne.0) then
            write(*,*) 'Can''t  perform MV multiplication'
            stop
         end if
         write(*,*) 'Errors of computed solution: ' 
         if(USE_MATRIX.eq.O_MATRIX) then
            write(*,*) 'Error : ',abs(yc-Cx)
         else if(USE_MATRIX.eq.T_MATRIX) then 
            write(*,*) 'Error : ',abs(yc-CTx)
         else if(USE_MATRIX.eq.H_MATRIX) then 
            write(*,*) 'Error : ',abs(yc-CHx)
         end if
         write(*,*) '* Test of MM multiplication *'
         select case(USE_MATRIX)
         case (O_MATRIX)
            call usmm(a,densec_B,densec_C,istat)
         case (T_MATRIX)
            call usmm(a,densec_B,densec_C,istat, transa =&
                 & TRANSP_MATRIX)
         case (H_MATRIX)
            call usmm(a,densec_B,densec_C,istat, transa =&
                 & HERMIT_MATRIX)
         case default
            stop 'No valid choice'
         end select
         if (istat.ne.0) then
            write(*,*) 'Can''t  perform MM multiplication'
            stop
         end if
         write(*,*) 'Errors of computed solution: ' 
         if(USE_MATRIX.eq.O_MATRIX) then
            write(*,*) 'Error : ',(abs(densec_C(:,i)-Cx),i=1,3)
         else if(USE_MATRIX.eq.T_MATRIX) then 
            write(*,*) 'Error : ',(abs(densec_C(:,i)-CTx),i=1,3)
         else if(USE_MATRIX.eq.H_MATRIX) then 
            write(*,*) 'Error : ',(abs(densec_C(:,i)-CHx),i=1,3)
         end if
         
         write(*,*) '*****************************'
         write(*,*) '* Deleting matrix handle    *'
         call usds(a,istat)
         deallocate(xc,yc,zc,densec_C,densec_B)
  ! **********************************************************************     
        print *,'*********************************************************'
         print *,'TEST WITH MATRIX D'
         USE_MATRIX = O_MATRIX
         prpty = blas_general + INDEX_BASE
         allocate(xc(3),yc(3))    
         call zuscr_begin(d_m,d_n,a,istat)
         call ussp(a,prpty,istat)
         call  uscr_insert_entries(a,VAL_d,INDX_d,JNDX_d,istat)
         call uscr_end(a,istat )
         do i=1,size(xc)
            xc(i) = (1.,1.)
         end do
      allocate(zc(size(xc)))     
      zc = xc
      yc = (0.,0.)
    allocate(densec_C(size(yc),3),densec_B(size(xc),3))
  
      do i = 1,3
         densec_B(:,i) = xc
         densec_C(:,i) = (0.,0.)
      end do
      write(*,*) '* Test of MV multiplication *'
      select case(USE_MATRIX)
      case (O_MATRIX)
         call usmv(a,xc,yc,istat)
      case (T_MATRIX)
         call usmv(a,xc,yc,istat, transa = TRANSP_MATRIX)
      case (H_MATRIX)
         call usmv(a,xc,yc,istat, transa = HERMIT_MATRIX)
      case default
         stop 'No valid choice'
      end select
      if (istat.ne.0) then
         write(*,*) 'Can''t  perform MV multiplication'
         stop
      end if
      write(*,*) 'Errors of computed solution: ' 
      if(USE_MATRIX.eq.O_MATRIX) then
         write(*,*) 'Error : ',abs(yc-dx)
      else if(USE_MATRIX.eq.T_MATRIX) then 
         write(*,*) 'Error : ',abs(yc-dTx)
      else if(USE_MATRIX.eq.H_MATRIX) then 
         write(*,*) 'Error : ',abs(yc-dHx)
      end if
      write(*,*) '* Test of MM multiplication *'
      select case(USE_MATRIX)

      case (O_MATRIX)
         call usmm(a,densec_B,densec_C,istat)
      case (T_MATRIX)
         call usmm(a,densec_B,densec_C,istat, transa = TRANSP_MATRIX)
      case (H_MATRIX)
         call usmm(a,densec_B,densec_C,istat, transa = HERMIT_MATRIX)
      case default
         stop 'No valid choice'
      end select
      if (istat.ne.0) then
         write(*,*) 'Can''t  perform MM multiplication'
         stop
      end if
      write(*,*) 'Errors of computed solution: ' 
      if(USE_MATRIX.eq.O_MATRIX) then
         write(*,*) 'Error : ',(abs(densec_C(:,i)-dx),i=1,3)
      else if(USE_MATRIX.eq.T_MATRIX) then 
         write(*,*) 'Error : ',(abs(densec_C(:,i)-dTx),i=1,3)
      else if(USE_MATRIX.eq.H_MATRIX) then 
         write(*,*) 'Error : ',(abs(densec_C(:,i)-dHx),i=1,3)
      end if
      write(*,*) '*****************************'
      write(*,*) '* Deleting matrix handle    *'
      call usds(a,istat)
      deallocate(xc,yc,zc,densec_C,densec_B)

      print *,'*****************************************************&
           &****'
      print *,'TEST WITH MATRIX DT'
      USE_MATRIX = O_MATRIX
      prpty = blas_upper_triangular + INDEX_BASE
      allocate(xc(3),yc(3))    
      call zuscr_begin(DT_m,DT_n,dt,istat)
      call ussp(dt,prpty,istat)
      call  uscr_insert_entries(dt,VAL_DT,INDX_DT,JNDX_DT,istat)
      call uscr_end(dt,istat )
      do i=1,size(xc)
         xc(i) = (1.,1.)
      end do
      allocate(zc(size(xc)))     
      zc = xc
      yc = (0.,0.)
      allocate(densec_C(size(yc),3),densec_B(size(xc),3))
      
      do i = 1,3
         densec_B(:,i) = xc
         densec_C(:,i) = (0.,0.)
      end do
      write(*,*) '* Test of MV multiplication *'
      select case(USE_MATRIX)
      case (O_MATRIX)
         call usmv(dt,xc,yc,istat)
      case (T_MATRIX)
         call usmv(dt,xc,yc,istat, transa = TRANSP_MATRIX)
      case (H_MATRIX)
         call usmv(dt,xc,yc,istat, transa = HERMIT_MATRIX)
      case default
         stop 'No valid choice'
      end select
      if (istat.ne.0) then
         write(*,*) 'Can''t  perform MV multiplication'
         stop
      end if
      write(*,*) 'Errors of computed solution: ' 
      if(USE_MATRIX.eq.O_MATRIX) then
         write(*,*) 'Error : ',abs(yc-DTTx)
      else if(USE_MATRIX.eq.T_MATRIX) then 
         write(*,*) 'Error : ',abs(yc-DTTx)
      else if(USE_MATRIX.eq.H_MATRIX) then 
         write(*,*) 'Error : ',abs(yc-DTHx)
      end if
      write(*,*) '* Test of MM multiplication *'
      select case(USE_MATRIX)

      case (O_MATRIX)
         call usmm(dt,densec_B,densec_C,istat)
      case (T_MATRIX)
         call usmm(dt,densec_B,densec_C,istat, transa = TRANSP_MATRIX)
      case (H_MATRIX)
         call usmm(dt,densec_B,densec_C,istat, transa = HERMIT_MATRIX)
      case default
         stop 'No valid choice'
      end select
      if (istat.ne.0) then
         write(*,*) 'Can''t  perform MM multiplication'
         stop
      end if
      write(*,*) 'Errors of computed solution: ' 
      if(USE_MATRIX.eq.O_MATRIX) then
         write(*,*) 'Error : ',(abs(densec_C(:,i)-DTTx),i=1,3)
      else if(USE_MATRIX.eq.T_MATRIX) then 
         write(*,*) 'Error : ',(abs(densec_C(:,i)-DTTx),i=1,3)
      else if(USE_MATRIX.eq.H_MATRIX) then 
         write(*,*) 'Error : ',(abs(densec_C(:,i)-DTHx),i=1,3)
      end if
    
      write(*,*) '* Test of tri. vec. solver  *'
      select case(USE_MATRIX)
         case (O_MATRIX)
            call ussv(dt,yc,istat)
         case (T_MATRIX)
            call ussv(dt,yc,istat, transa = TRANSP_MATRIX)
         case (H_MATRIX)
            call ussv(dt,yc,istat, transa = HERMIT_MATRIX)
         case default
            stop 'No valid choice'
         end select
         if (istat.ne.0) then
            write(*,*) 'Can''t  perform triangular solve'
            stop
         else
            write(*,*) 'Errors of computed solution: ' 
            write(*,*) 'Error : ',abs(yc-xc)
         end if
         write(*,*) '* Testing tri. mat. solver  *'
         select case(USE_MATRIX)
         case (O_MATRIX)
            call ussm(dt,densec_C,istat)
         case (T_MATRIX)
            call ussm(dt,densec_C,istat, transa = TRANSP_MATRIX)
         case (H_MATRIX)
            call ussm(dt,densec_C,istat, transa = HERMIT_MATRIX)
         case default
            stop 'No valid choice'
         end select
         if (istat.ne.0) then
            write(*,*) 'Can''t  perform triangular solve'
            stop
         else
            write(*,*) 'Errors of computed solution: ' 
            write(*,*) 'Error : ',abs(densec_C-densec_B)
         end if

         write(*,*) '*****************************'
         write(*,*) '* Deleting matrix handle    *'
         call usds(dt,istat)
         deallocate(xc,yc,zc,densec_C,densec_B)

        print *,'*********************************************************'
         print *,'TEST WITH MATRIX DH'
         USE_MATRIX = O_MATRIX
         prpty = blas_hermitian + blas_upper + INDEX_BASE
         allocate(xc(3),yc(3))    
         call zuscr_begin(DT_m,DT_n,dt,istat)
         call ussp(dt,prpty,istat)
         call  uscr_insert_entries(dt,VAL_DT,INDX_DT,JNDX_DT,istat)
         call uscr_end(dt,istat )
         do i=1,size(xc)
            xc(i) = (1.,1.)
         end do
         allocate(zc(size(xc)))     
         zc = xc
         yc = (0.,0.)
         allocate(densec_C(size(yc),3),densec_B(size(xc),3))
         
         do i = 1,3
            densec_B(:,i) = xc
            densec_C(:,i) = (0.,0.)
         end do
         write(*,*) '* Test of MV multiplication *'
         select case(USE_MATRIX)
      case (O_MATRIX)
         call usmv(dt,xc,yc,istat)
      case (T_MATRIX)
         call usmv(dt,xc,yc,istat, transa = TRANSP_MATRIX)
      case (H_MATRIX)
         call usmv(dt,xc,yc,istat, transa = HERMIT_MATRIX)
      case default
         stop 'No valid choice'
      end select
      if (istat.ne.0) then
         write(*,*) 'Can''t  perform MV multiplication'
         stop
      end if
      write(*,*) 'Errors of computed solution: ' 
      if(USE_MATRIX.eq.O_MATRIX) then
         write(*,*) 'Error : ',abs(yc-dx)
      else if(USE_MATRIX.eq.T_MATRIX) then 
         write(*,*) 'Error : ',abs(yc-DTHx)
      else if(USE_MATRIX.eq.H_MATRIX) then 
         write(*,*) 'Error : ',abs(yc-DTHx)
      end if
      write(*,*) '* Test of MM multiplication *'
      select case(USE_MATRIX)
      case (O_MATRIX)
         call usmm(dt,densec_B,densec_C,istat)
      case (T_MATRIX)
         call usmm(dt,densec_B,densec_C,istat, transa = TRANSP_MATRIX)
      case (H_MATRIX)
         call usmm(dt,densec_B,densec_C,istat, transa = HERMIT_MATRIX)
      case default
         stop 'No valid choice'
      end select
      if (istat.ne.0) then
         write(*,*) 'Can''t  perform MM multiplication'
         stop
      end if
      write(*,*) 'Errors of computed solution: ' 
      if(USE_MATRIX.eq.O_MATRIX) then
         write(*,*) 'Error : ',(abs(densec_C(:,i)-dx),i=1,3)
      else if(USE_MATRIX.eq.T_MATRIX) then 
         write(*,*) 'Error : ',(abs(densec_C(:,i)-DTHx),i=1,3)
      else if(USE_MATRIX.eq.H_MATRIX) then 
         write(*,*) 'Error : ',(abs(densec_C(:,i)-DTHx),i=1,3)
      end if
  write(*,*) '*****************************'
      write(*,*) '* Deleting matrix handle    *'
      call usds(dt,istat)
 deallocate(xc,yc,zc,densec_C,densec_B)
  ! **********************************************************************     

        print *,'*********************************************************'
         print *,'TEST WITH MATRIX I'
         prpty = blas_general + INDEX_BASE
         allocate(x_i(5),y_i(5))
         y_i=0
         call iuscr_begin(I_m,I_n,a,istat)
         call ussp(a,prpty,istat)
         do i=1,size(I_val)
            call uscr_insert_entry(a,I_val(i),I_indx(i),I_jndx(i),istat)
         end do
         call uscr_end(a,istat )
         do i=1,size(x_i)
            x_i(i) = int(i)
         end do
         allocate(z_i(size(x_i)))
         z_i = x_i
         write(*,*) '* Test of MV multiplication *'
         call usmv(a,x_i,y_i,istat)
         write(*,*) 'Error : ',abs(y_i-Ix)
         write(*,*) '* Deleting matrix handle    *'
         deallocate(x_i,y_i,z_i)
         call usds(a, istat)
         if (istat.ne.0) then
            write(*,*) 'Deallocation failure'
            stop
         end if

  ! **********************************************************************     
        print *,'*********************************************************'
         print *,'TEST WITH MATRIX J'  
         prpty = blas_upper_triangular + INDEX_BASE
         allocate(x_i(5),y_i(5))
         y_i=0
         call iuscr_begin(J_m,J_n,a,istat)
         call uscr_insert_entries(a,J_VAL,J_indx,J_jndx,istat)
         call ussp(a,prpty,istat)
         call uscr_end(a,istat )
         do i=1,size(x_i)
            x_i(i) = int(i)
         end do
         allocate(z_i(size(x_i)))
         z_i = x_i
         write(*,*) '* Test of MV multiplication *'
         call usmv(a,x_i,y_i,istat)
         write(*,*) 'Error : ',abs(y_i-Jx)
         write(*,*) '* Test of tri. vec. solver  *'
         call ussv(a,y_i,istat)
         if (istat.ne.0) then
            write(*,*) 'Can''t  perform triangular solve'
            stop
         else
            write(*,*) 'Errors of computed solution: ' 
            write(*,*) 'Error : ',abs(y_i-x_i)
         end if
         write(*,*) '* Deleting matrix handle    *'
         deallocate(x_i,y_i,z_i)
         call usds(a, istat)
         if (istat.ne.0) then
            write(*,*) 'Deallocation failure'
            stop
         end if

  ! **********************************************************************     
       print *,'*********************************************************'
        print *,'TEST WITH MATRIX M'
        prpty = blas_upper_triangular + INDEX_BASE + BLOCK_INTERN
        allocate(x(4),y(4))
         y=0.
        call duscr_block_begin(M_Mb,M_Nb,M_k,M_l,a,istat)
        call ussp(a,prpty,istat)
        call uscr_insert_block(a,M11,1,1,istat)
        call uscr_insert_block(a,M12,1,2,istat)
        call uscr_insert_block(a,M22,2,2,istat)
        call uscr_end(a,istat )
        do i=1,size(x)
            x(i) = dble(i)
         end do
         allocate(z(size(x)))
         z = x
         allocate(dense_C(size(y),3),dense_B(size(x),3))
         do i = 1,3
            dense_B(:,i) = x
            dense_C(:,i) = 0.
         end do
         write(*,*) '* Test of MV multiplication *'
         call usmv(a,x,y,istat)
         write(*,*) 'Error : ',abs(y-Mx)
         write(*,*) '* Test of MM multiplication *'
         call usmm(a,dense_B,dense_C,istat)
         write(*,*) 'Error : ',(abs(dense_C(:,i)-Mx),i=1,3)
         write(*,*) '* Deleting matrix handle    *'
         deallocate(x,y,dense_C,dense_B,z)
         call usds(a, istat)
         if (istat.ne.0) then
            write(*,*) 'Deallocation failure'
            stop
         end if
  ! **********************************************************************     
        print *,'*********************************************************'
         print *,'TEST WITH MATRIX N'
        prpty = blas_lower_triangular + INDEX_BASE + BLOCK_INTERN
        allocate(x(4),y(4))
        y=0.
        call duscr_block_begin(N_Mb,N_Nb,N_k,N_l,a,istat)
        call ussp(a,prpty,istat)
        call uscr_insert_block(a,N11,1,1,istat)
        call uscr_insert_block(a,N21,2,1,istat)
        call uscr_insert_block(a,N22,2,2,istat)
        call uscr_end(a,istat )
        do i=1,size(x)
           x(i) = dble(i)
        end do
        allocate(z(size(x)))
        z = x
        allocate(dense_C(size(y),3),dense_B(size(x),3))
        do i = 1,3
           dense_B(:,i) = x
           dense_C(:,i) = 0.
        end do
        write(*,*) '* Test of MV multiplication *'
        call usmv(a,x,y,istat)
        write(*,*) 'Error : ',abs(y-Nx)
        write(*,*) '* Test of MM multiplication *'
        call usmm(a,dense_B,dense_C,istat)
        write(*,*) 'Error : ',(abs(dense_C(:,i)-Nx),i=1,3)
        write(*,*) '* Deleting matrix handle    *'
        deallocate(x,y,dense_C,dense_B,z)
        call usds(a, istat)
        if (istat.ne.0) then
           write(*,*) 'Deallocation failure'
           stop
        end if
  ! **********************************************************************     
        print *,'*********************************************************'
         print *,'TEST WITH MATRIX T'  
         prpty = blas_upper_triangular + INDEX_BASE
         allocate(x(5),y(5))
         y=0.
         call duscr_begin(T_m,T_n,a,istat)
         call uscr_insert_entries(a,T_VAL,T_indx,T_jndx,istat)
         call ussp(a,prpty,istat)
         call uscr_end(a,istat )
         do i=1,size(x)
            x(i) = dble(i)
         end do
         allocate(z(size(x)))
         z = x
         allocate(dense_C(size(y),3),dense_B(size(x),3))
         do i = 1,3
            dense_B(:,i) = x
            dense_C(:,i) = 0.
         end do
         write(*,*) '* Test of MV multiplication *'
         call usmv(a,x,y,istat)
         write(*,*) 'Error : ',abs(y-Tx)
         write(*,*) '* Test of MM multiplication *'
         call usmm(a,dense_B,dense_C,istat)
         write(*,*) 'Error : ',(abs(dense_C(:,i)-Tx),i=1,3)
         write(*,*) '* Test of tri. vec. solver  *'
         call ussv(a,y,istat)
         if (istat.ne.0) then
            write(*,*) 'Can''t  perform triangular solve'
            stop
         else
            write(*,*) 'Errors of computed solution: ' 
            write(*,*) 'Error : ',abs(y-x)
         end if
         write(*,*) '* Testing tri. mat. solver  *'
         call ussm(a,dense_C,istat)
         if (istat.ne.0) then
            write(*,*) 'Can''t  perform triangular solve'
            stop
         else
            write(*,*) 'Errors of computed solution: ' 
            write(*,*) 'Error : ',abs(dense_C-dense_B)
         end if
         write(*,*) '* Deleting matrix handle    *'
         deallocate(x,y,dense_C,dense_B,z)
         call usds(a, istat)
         if (istat.ne.0) then
            write(*,*) 'Deallocation failure'
            stop
         end if
  ! **********************************************************************     
        print *,'*********************************************************'
         print *,'TEST WITH MATRIX U'
         prpty = blas_lower_triangular + INDEX_BASE
         allocate(x(5),y(5))
         y=0.
         call duscr_begin(U_m,U_n,a,istat)
         call ussp(a,prpty,istat)
         call uscr_insert_entries(a,U_VAL,U_indx,U_jndx,istat)   
         call uscr_end(a,istat )
         do i=1,size(x)
            x(i) = dble(i)
         end do
         allocate(z(size(x)))
         z = x
         allocate(dense_C(size(y),3),dense_B(size(x),3))
         do i = 1,3
            dense_B(:,i) = x
            dense_C(:,i) = 0.
         end do
         write(*,*) '* Test of MV multiplication *'
         call usmv(a,x,y,istat)
         write(*,*) 'Error : ',abs(y-Ux)
         write(*,*) '* Test of MM multiplication *'
         call usmm(a,dense_B,dense_C,istat)
         write(*,*) 'Error : ',(abs(dense_C(:,i)-Ux),i=1,3)
         write(*,*) '* Test of tri. vec. solver  *'
         call ussv(a,y,istat)
         if (istat.ne.0) then
            write(*,*) 'Can''t  perform triangular solve'
            stop
         else
            write(*,*) 'Errors of computed solution: ' 
            write(*,*) 'Error : ',abs(y-x)
         end if
         write(*,*) '* Testing tri. mat. solver  *'
         call ussm(a,dense_C,istat)
         if (istat.ne.0) then
            write(*,*) 'Can''t  perform triangular solve'
            stop
         else
            write(*,*) 'Errors of computed solution: ' 
            write(*,*) 'Error : ',abs(dense_C-dense_B)
         end if
         write(*,*) '* Deleting matrix handle    *'
         deallocate(x,y,dense_C,dense_B,z)
         call usds(a, istat)
         if (istat.ne.0) then
            write(*,*) 'Deallocation failure'
            stop
         end if
  ! **********************************************************************     
        print *,'*********************************************************'
         print *,'TEST WITH MATRIX X'
         prpty = blas_general + INDEX_BASE + BLOCK_INTERN
         allocate(x(8),y(6))
         y=0.
         call duscr_block_begin(X_Mb,X_Nb,X_k,X_l,a,istat)
         call ussp(a,prpty,istat)
         call uscr_insert_block(a,X11,1,1,istat)
         call uscr_insert_block(a,X31,3,1,istat)
         call uscr_insert_block(a,X22,2,2,istat)
         call uscr_insert_block(a,X13,1,3,istat)
         call uscr_insert_block(a,X23,2,3,istat)
         call uscr_end(a,istat )
         do i=1,size(x)
            x(i) = dble(i)
         end do
         allocate(z(size(x)))
         z = x
         allocate(dense_C(size(y),3),dense_B(size(x),3))
         do i = 1,3
            dense_B(:,i) = x
            dense_C(:,i) = 0.
         end do
         write(*,*) '* Test of MV multiplication *'
         call usmv(a,x,y,istat)
         write(*,*) 'Error : ',abs(y-Xx)
         write(*,*) '* Test of MM multiplication *'
         call usmm(a,dense_B,dense_C,istat)
         write(*,*) 'Error : ',(abs(dense_C(:,i)-Xx),i=1,3)
         write(*,*) '* Deleting matrix handle    *'
         deallocate(x,y,dense_C,dense_B,z)
         call usds(a, istat)
         if (istat.ne.0) then
            write(*,*) 'Deallocation failure'
            stop
         end if
  ! **********************************************************************     
        print *,'*********************************************************'
         print *,'TEST WITH MATRIX Y'
         prpty = blas_general + INDEX_BASE + BLOCK_INTERN
         allocate(x(6),y(6))
         y=0.
         call duscr_block_begin(Y_Mb,Y_Nb,Y_k,Y_l,a,istat)
         call ussp(a,prpty,istat)
         call uscr_insert_block(a,Y11,1,1,istat)
         call uscr_insert_block(a,Y12,1,2,istat)
         call uscr_insert_block(a,Y22,2,2,istat)
         call uscr_insert_block(a,Y21,2,1,istat)
         call uscr_end(a,istat )
         do i=1,size(x)
            x(i) = dble(i)
         end do
         allocate(z(size(x)))
         z = x
         allocate(dense_C(size(y),3),dense_B(size(x),3))
         do i = 1,3
            dense_B(:,i) = x
            dense_C(:,i) = 0.
         end do
         write(*,*) '* Test of MV multiplication *'
         call usmv(a,x,y,istat)
         write(*,*) 'Error : ',abs(y-Yx)
         write(*,*) '* Test of MM multiplication *'
         call usmm(a,dense_B,dense_C,istat)
         write(*,*) 'Error : ',(abs(dense_C(:,i)-Yx),i=1,3)
         write(*,*) '* Deleting matrix handle    *'
         deallocate(x,y,dense_C,dense_B,z)
         call usds(a, istat)
         if (istat.ne.0) then
            write(*,*) 'Deallocation failure'
            stop
         end if
        print *,'*********************************************************'
         print *,'TEST WITH MATRIX Y_SU'
         prpty = blas_symmetric + INDEX_BASE+blas_upper+BLOCK_INTERN
         allocate(x(6),y(6))
         y=0.
         call duscr_block_begin(Y_SU_Mb,Y_SU_Nb,Y_SU_k,Y_SU_l,a,istat)
         call ussp(a,prpty,istat)
         call uscr_insert_block(a,Y_SU11,1,1,istat)
         call uscr_insert_block(a,Y_SU12,1,2,istat)
         call uscr_insert_block(a,Y_SU22,2,2,istat)
         call uscr_insert_block(a,Y_SU21,2,1,istat)
         call uscr_end(a,istat )
         do i=1,size(x)
            x(i) = dble(i)
         end do
         allocate(z(size(x)))
         z = x
         allocate(dense_C(size(y),3),dense_B(size(x),3))
         do i = 1,3
            dense_B(:,i) = x
            dense_C(:,i) = 0.
         end do
         write(*,*) '* Test of MV multiplication *'
         call usmv(a,x,y,istat)
         write(*,*) 'Error : ',abs(y-Y_SUx)
         write(*,*) '* Test of MM multiplication *'
         call usmm(a,dense_B,dense_C,istat)
         write(*,*) 'Error : ',(abs(dense_C(:,i)-Y_SUx),i=1,3)
         write(*,*) '* Deleting matrix handle    *'
         deallocate(x,y,dense_C,dense_B,z)
         call usds(a, istat)
         if (istat.ne.0) then
            write(*,*) 'Deallocation failure'
            stop
         end if
        print *,'*********************************************************'
         print *,'TEST WITH MATRIX Y_SL'
         prpty = blas_symmetric + INDEX_BASE+blas_lower+BLOCK_INTERN
         allocate(x(6),y(6))
         y=0.
         call duscr_block_begin(Y_SL_Mb,Y_SL_Nb,Y_SL_k,Y_SL_l,a,istat)
         call ussp(a,prpty,istat)
         call uscr_insert_block(a,Y_SL11,1,1,istat)
         call uscr_insert_block(a,Y_SL12,1,2,istat)
         call uscr_insert_block(a,Y_SL22,2,2,istat)
         call uscr_insert_block(a,Y_SL21,2,1,istat)
         call uscr_end(a,istat )
         do i=1,size(x)
            x(i) = dble(i)
         end do
         allocate(z(size(x)))
         z = x
         allocate(dense_C(size(y),3),dense_B(size(x),3))
         do i = 1,3
            dense_B(:,i) = x
            dense_C(:,i) = 0.
         end do
         write(*,*) '* Test of MV multiplication *'
         call usmv(a,x,y,istat)
         write(*,*) 'Error : ',abs(y-Y_SLx)
         write(*,*) '* Test of MM multiplication *'
         call usmm(a,dense_B,dense_C,istat)
         write(*,*) 'Error : ',(abs(dense_C(:,i)-Y_SLx),i=1,3)
         write(*,*) '* Deleting matrix handle    *'
         deallocate(x,y,dense_C,dense_B,z)
         call usds(a, istat)
         if (istat.ne.0) then
            write(*,*) 'Deallocation failure'
            stop
         end if
  ! **********************************************************************     
        print *,'*********************************************************'
         print *,'TEST WITH MATRIX Z'
         prpty = blas_general + INDEX_BASE + BLOCK_INTERN
         allocate(x(11),y(11))
         y=0.
         call duscr_variable_block_begin(Z_Mb,Z_Nb,Z_kk,Z_ll,a,istat)
         call ussp(a,prpty,istat)
         call uscr_insert_block(a,Z11,1,1,istat)
         call uscr_insert_block(a,Z13,1,3,istat)
         call uscr_insert_block(a,Z22,2,2,istat)
         call uscr_insert_block(a,Z32,3,2,istat)
         call uscr_insert_block(a,Z23,2,3,istat)
         call uscr_insert_block(a,Z31,3,1,istat)
         call uscr_insert_block(a,Z33,3,3,istat)
         call uscr_insert_block(a,Z55,5,5,istat)
         call uscr_insert_block(a,Z44,4,4,istat)
         call uscr_insert_block(a,Z51,5,1,istat)
         call uscr_insert_block(a,Z15,1,5,istat)
         call uscr_insert_block(a,Z34,3,4,istat)
         call uscr_insert_block(a,Z43,4,3,istat)
         call uscr_end(a,istat )
         do i=1,size(x)
            x(i) = dble(i)
         end do
         allocate(z(size(x)))
         z = x
         allocate(dense_C(size(y),3),dense_B(size(x),3))
         do i = 1,3
            dense_B(:,i) = x
            dense_C(:,i) = 0.
         end do
         write(*,*) '* Test of MV multiplication *'
         call usmv(a,x,y,istat)
         write(*,*) 'Error : ',abs(y-Zx)
         write(*,*) '* Test of MM multiplication *'
         call usmm(a,dense_B,dense_C,istat)
         write(*,*) 'Error : ',(abs(dense_C(:,i)-Zx),i=1,3)
         write(*,*) '* Deleting matrix handle    *'
         deallocate(x,y,dense_C,dense_B,z)
         call usds(a, istat)
         if (istat.ne.0) then
            write(*,*) 'Deallocation failure'
            stop
         end if
  ! **********************************************************************     

      write(*,*) '*****************************'
      write(*,*) '* REGULAR END OF PROGRAM    *' 
      write(*,*) '*****************************'
      ierr=0
            
end program tester
     

      module mod_usaxpy
      use blas_sparse_namedconstants
      interface usaxpy
      module procedure iusaxpy
      module procedure susaxpy
      module procedure dusaxpy
      module procedure cusaxpy
      module procedure zusaxpy
      end interface
      contains
! **********************************************************************
! **********************************************************************
 subroutine iusaxpy (x,indx,y,alpha)
    integer ,dimension(:),intent(in) ::x
    integer ,dimension(:),intent(inout) ::y
    integer,dimension(:),intent(in) ::indx
    integer ,intent(in) ,optional ::alpha
    integer  :: i,t
    t=size(indx)
    if(t.gt.0) then             
       if(present(alpha)) then                          
          do i=1,t
             y(indx(i))=y(indx(i))+x(i)*alpha
          end do
       else
          do i=1,t
             y(indx(i))=y(indx(i))+x(i)
          end do
       end if
    end if
  end subroutine iusaxpy 
! **********************************************************************
! **********************************************************************
 subroutine susaxpy (x,indx,y,alpha)
    real(KIND=sp) ,dimension(:),intent(in) ::x
    real(KIND=sp) ,dimension(:),intent(inout) ::y
    integer,dimension(:),intent(in) ::indx
    real(KIND=sp) ,intent(in) ,optional ::alpha
    integer  :: i,t
    t=size(indx)
    if(t.gt.0) then             
       if(present(alpha)) then                          
          do i=1,t
             y(indx(i))=y(indx(i))+x(i)*alpha
          end do
       else
          do i=1,t
             y(indx(i))=y(indx(i))+x(i)
          end do
       end if
    end if
  end subroutine susaxpy 
! **********************************************************************
! **********************************************************************
 subroutine dusaxpy (x,indx,y,alpha)
    real(KIND=dp) ,dimension(:),intent(in) ::x
    real(KIND=dp) ,dimension(:),intent(inout) ::y
    integer,dimension(:),intent(in) ::indx
    real(KIND=dp) ,intent(in) ,optional ::alpha
    integer  :: i,t
    t=size(indx)
    if(t.gt.0) then             
       if(present(alpha)) then                          
          do i=1,t
             y(indx(i))=y(indx(i))+x(i)*alpha
          end do
       else
          do i=1,t
             y(indx(i))=y(indx(i))+x(i)
          end do
       end if
    end if
  end subroutine dusaxpy 
! **********************************************************************
! **********************************************************************
 subroutine cusaxpy (x,indx,y,alpha)
    complex(KIND=sp) ,dimension(:),intent(in) ::x
    complex(KIND=sp) ,dimension(:),intent(inout) ::y
    integer,dimension(:),intent(in) ::indx
    complex(KIND=sp) ,intent(in) ,optional ::alpha
    integer  :: i,t
    t=size(indx)
    if(t.gt.0) then             
       if(present(alpha)) then                          
          do i=1,t
             y(indx(i))=y(indx(i))+x(i)*alpha
          end do
       else
          do i=1,t
             y(indx(i))=y(indx(i))+x(i)
          end do
       end if
    end if
  end subroutine cusaxpy 
! **********************************************************************
! **********************************************************************
 subroutine zusaxpy (x,indx,y,alpha)
    complex(KIND=dp) ,dimension(:),intent(in) ::x
    complex(KIND=dp) ,dimension(:),intent(inout) ::y
    integer,dimension(:),intent(in) ::indx
    complex(KIND=dp) ,intent(in) ,optional ::alpha
    integer  :: i,t
    t=size(indx)
    if(t.gt.0) then             
       if(present(alpha)) then                          
          do i=1,t
             y(indx(i))=y(indx(i))+x(i)*alpha
          end do
       else
          do i=1,t
             y(indx(i))=y(indx(i))+x(i)
          end do
       end if
    end if
  end subroutine zusaxpy 
! **********************************************************************
! **********************************************************************
      end module mod_usaxpy

 subroutine XUSAXPY(x,indx,y,alpha)
    DCOMPLEX,dimension(:),intent(in) ::x
    DCOMPLEX,dimension(:),intent(inout) ::y
    integer,dimension(:),intent(in) ::indx
    DCOMPLEX,intent(in) ,optional ::alpha
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
  end subroutine XUSAXPY

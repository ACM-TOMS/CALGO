
  subroutine XUSSC(x,y,indx)
 DCOMPLEX,dimension(:),intent(in) ::x
 DCOMPLEX ,dimension(:),intent(inout) ::y
  integer,dimension(:),intent(in)::indx
  integer  :: i,t
     
  t=size(indx)
  if(t.gt.0) then      
     do i=1,t
        y(indx(i))=  x(i)
     end do
  end if
  
  end subroutine XUSSC

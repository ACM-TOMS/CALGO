 subroutine XUSGA(y,x,indx)
    DCOMPLEX ,dimension(:),intent(inout) ::x
    DCOMPLEX ,dimension(:),intent(in) ::y
    integer,dimension(:),intent(in)::indx
    integer  ::t,i 
    t=size(x) 
    if(t.gt.0) then
       do i=1,t
          x(i)=y(indx(i))        
       end do
    end if
  end subroutine XUSGA

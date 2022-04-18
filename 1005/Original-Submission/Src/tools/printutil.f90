MODULE PRINTUTIL
CONTAINS
  SUBROUTINE PRINTMAT(S, A)
    character(*), intent(in) :: s
    character(20) :: blk
    integer :: i
    double precision, intent(in) :: A(:,:)
    blk = ' '
    print '(A,10F11.4)', s, A(1,:)
    do i=2,size(A,1)
      print '(A,10F11.4)', blk(1:len(s)), A(i,:)
    enddo
    print *
  END SUBROUTINE PRINTMAT

  SUBROUTINE PRINTVEC(S, X)
    character(*), intent(in) :: s
    character(*), parameter :: blk = '                    '!
    integer :: i
    double precision, intent(in) :: x(:)
    print '(A,F11.4)', s, x(1)
    print '(A,F11.4)', (blk(1:len(s)), x(i), i=2, size(x))
    print *
  END SUBROUTINE PRINTVEC
END MODULE PRINTUTIL

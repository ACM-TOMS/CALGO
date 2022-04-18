program test_dlyap
  !use dispmodule
  implicit none
  integer, parameter :: dble = kind(0d0)
  real(dble), allocatable :: A(:,:), Sig(:,:), S(:,:), C(:,:)
  real(dble) :: d
  integer :: n, rep, i, itr
  character(2) :: transp
  transp = 'nt'
  do n=1,5
    allocate(A(n,n), Sig(n,n), S(n,n), C(n,n))
    do rep=1,5
      do itr = 1,2
        call random_number(A)
        call dscal(n*n, 0.1d0, A, 1)
        do i=1,n ! make Sig symmetric
          call random_number(Sig(i,1:i))
          Sig(1:i-1,i) = Sig(i,1:i-1)
        enddo
        call dlyap_slicot(transp(itr:itr), n, A, Sig, S)
        d = 0d0
        do i=1,n
          d = max(d, maxval(abs(Sig(1:i,i) - Sig(i,1:i))))
        enddo
        if (d > 0d0) then
          print *,'Error: Non-symmetric solution from dlyap'
          stop 1
        end if
        if (itr==1) then ! Compute residual, Sig - (S - A*S*A')
          call dsymm('right', 'low', n, n,    1d0, S, n, A, n,  0d0, C, n)
          call dgemm('notr', 'tran', n, n, n, 1d0, C, n, A, n, -1d0, S, n)
        else ! ... or Sig - (S - A'*S*A)
          call dsymm('left',  'low', n, n,    1d0, S, n, A, n,  0d0, C, n)
          call dgemm('tran', 'notr', n, n, n, 1d0, C, n, A, n, -1d0, S, n)
        endif
        S = S + Sig
        if (maxval(abs(S)) > 1d-8) then
          print *,'Error: Wrong solution from dlyap'
          print *,'n=',n,', rep=',rep,', itr=',itr
          stop 1
        endif
      enddo
    enddo
    deallocate(A, Sig, S, C)
  enddo
end program test_dlyap

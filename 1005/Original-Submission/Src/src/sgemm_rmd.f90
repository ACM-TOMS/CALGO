SUBROUTINE SGEMM_RMD(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, ldc, Aa, Ba, Ca, sel)
  implicit none
  real :: alpha, beta
  integer :: m,n,k, lda,ldb,ldc, i
  character :: transa, transb
  real :: A(lda,*), B(ldb,*), Aa(lda,*), Ba(ldb,*), Ca(ldc,*)
  character(3) :: sel

  ! PURPOSE
  !    Calculate the reverse mode derivative of SGEMM from BLAS.
  !
  ! ARGUMENTS
  !    If SGEMM was called with the arguments
  !    
  !       transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc
  !    
  !    then the corresponding call to SGEMM_RMD should begin with the arguments
  !
  !       transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, ldc
  !
  !    with the same values. Note that C is omitted. All these arguments will
  !    remain unchanged on exit. In addition the following arguments should
  !    be provided:
  !
  !    Aa   (input, output, real matrix of the same dimensions as A)
  !         Aa += the adjoint of A due to the SGEMM call.
  !
  !    Ba   (input, output, real matrix of the same dimensions as B)
  !         Ba += the adjoint of B due to the SGEMM call.
  !
  !    Ca   (input, output, real matrix of the same dimensions as C)
  !         On entry: the adjoint of the C produced by SGEMM
  !         On exit: the adjoint of the C supplied to SGEMM
  !
  !    sel  (input, character*3)
  !         Used to select which adjoints to update/compute:
  !            sel(1:1) = '1' if Aa should be updated, else sel(1:1) = '0'
  !            sel(2:2) = '1' if Ba should be updated, else sel(2:2) = '0'
  !            sel(3:3) = '1' if Ca should be computed, else sel(3:3) = '0'
  !         For example, to update only Aa, set sel = '100'.
  !       
  ! NOTE
  !    To compute C := alpha*A*A + beta*C one may call sgemm with a repeated
  !    argument, e.g.
  !       call sgemm('N', 'N', n, n, n, 1.0, A, n, A, n, 1.0, C, n).
  !    The correct adjoint of A, Aa = Ca*A' + A'*Ca, can then be computed
  !    with:
  !       call sgemm_rmd('N','N',n,n,n,1.0,A,n,A,n,1.0,n,Aa,dummy,dummy,'100')
  !       call sgemm_rmd('N','N',n,n,n,1.0,A,n,A,n,1.0,n,dummy,Aa,dummy,'010')
  !       
  ! OPERATIONS
  !          SGEMM('N', 'N'...)        SGEMM('T', 'N'...)
  !    BLAS: C = alpha*A*B + beta*C    C = alpha*A'*B + beta*C
  !    RMD:  Aa += alpha*Ca*B'         Aa += alpha*B*Ca'
  !          Ba += alpha*A'*Ca         Ba += alpha*A*Ca
  !          Ca := beta*Ca             Ca := beta*Ca
  !
  !          SGEMM('N', 'T'...)        SGEMM('T', 'T'...)
  !    BLAS: C = alpha*A*B' + beta*C   C = alpha*A'*B' + beta*C
  !    RMD:  Aa += alpha*Ca*B          Aa += alpha*B'*Ca'
  !          Ba += alpha*Ca'*A         Ba += alpha*Ca'*A'
  !          Ca := beta*Ca             Ca := beta*Ca
  !    The Ca on the right hand sides of the equals signs is its value on entry

  ! Local variables
  logical sela, selb, selc
  character :: ntransa, ntransb, trno(2) = ['T', 'N']

  sela = sel(1:1) == '1'
  selb = sel(2:2) == '1'
  selc = sel(3:3) == '1'
  ntransa = trno(1 + (index('NnTt', transa)-1)/2) ! N-->T, T-->N
  ntransb = trno(1 + (index('NnTt', transb)-1)/2) ! N-->T, T-->N

  ! BLAS operation
  !   C = alpha*op(A)*op(B) + beta*C

  if(sela) then
    ! Calculate Aa
    if(transa == 'n' .or. transa == 'N') then
      ! A is not transposed, i.e. op(A) = A
      ! RMD operation
      !   Aa += alpha*Ca*(op(B))'
      call sgemm('n',ntransb, m,k,n, alpha,Ca,ldc,B,ldb, 1.0,Aa,lda)
    else
      ! A is transposed, i.e. op(A) = A'
      ! RMD operation
      !   Aa += alpha*(op(B))*Ca'
      call sgemm(transb,'t', k,m,n, alpha,B,ldb,Ca,ldc, 1.0,Aa,lda)
    end if
  end if
    
  if(selb) then
    ! Calculate Ba
    if(transb == 'n' .or. transb == 'N') then
      ! B is not transposed, i.e. op(B) = B
      ! RMD operation
      !   Ba +=alpha*(op(A))'*Ca
      call sgemm(ntransa,'n', k,n,m, alpha,A,lda,Ca,ldc, 1.0,Ba,ldb)
    else 
      ! B is transposed, i.e. op(B) = B'
      ! RMD operation
      !   Ba += alpha*Ca'*(op(A))
      call sgemm('t',transa, n,k,m, alpha,Ca,ldc,A,lda, 1.0,Ba,ldb)
    end if
  end if
    
  if(selc) then
    ! Calculate Ca
    ! Note: This must be done last, since Aa and Ba rely on the old 
    ! value of Ca.
    ! RMD operation
    !   Ca = beta*Ca
    if (m == ldc) then
      call sscal(m*n, beta, Ca, 1)
    else
      do i = 1,n
        call sscal(m, beta, Ca(1,i), 1)
      enddo
    end if
    !call sgemm(transa,transb, m,n,k, 0.0,A,lda,B,ldb, beta,Ca,ldc)
  end if

END SUBROUTINE SGEMM_RMD

! This module implements basic operations for overlapped block-diagonal
! matrices that consist of upper quasi-triangular blocks.
!
! References:
!
!   _. _______,
!   On a sep-inverse estimate and its application to a modified
!   block Schur--Parlett algorithm,
!   submitted to ACM transactions on mathematical software, 2008.
!
! This module is intended for internal use only.

module ovlpmatrix

  use floattypes
  use matrixpwrtag

  implicit none

  public

contains


!    Name : ovlp_uniqsort
!         :
! Purpose : This subroutine eliminates all diagonal blocks that are completely
!         : included in some of the other diagonal blocks. After this
!         : elimination, the diagonal blocks are sorted according to their
!         : position along the diagonal.
!         :
!   Input : - ``m'' is the number of diagonal blocks.
!         : - diag(i,1) is the row/column number of the head of the i-th
!         :   diagonal block, and diag(i,2) is the row/column number of
!         :   the tail of the i-th diagonal block.
!         :
!  Output : - diag(*,*) stores the information on the diagonal blocks
!         :   after the above elimination and sorting.
!         : - ``n'' is the number of diagonal blocks after the above operation.

pure subroutine ovlp_uniqsort(m,diag,n)
  integer, intent(in   ) :: m                       ! the # of input  blocks
  integer, intent(inout) :: diag(:,:)
  integer, intent(out  ) :: n                       ! the # of output blocks

  integer :: wrk(m,size(diag,2))
  integer :: i, j, ptr, count, tmp(size(diag,2))
  integer :: included

  wrk = 0; n = 1
  do i=1,m                                          ! eliminate included ones
    included = 0

    do j=1,m           
      if (i == j) cycle
      if (                                        &
           ( diag(i,1) >  diag(j,1) .and.         &
             diag(i,2) <= diag(j,2)       ) .or.  &
           ( diag(i,1) >= diag(j,1) .and.         &
             diag(i,2) <  diag(j,2)       )       &
      ) then
        included = 1; exit
      end if
    end do

    if (included == 0) then
      wrk(n,:) = diag(i,:)
      n = n + 1
    end if
  end do

  n = n - 1

  do                                                ! bubble sort
    count = 0
    do i=1,n-1
      if (wrk(i,1) > wrk(i+1,1)) then
        tmp(:)     = wrk(i+1,:)
        wrk(i+1,:) = wrk(i  ,:)
        wrk(i  ,:) = tmp(:)
        count = count + 1
      end if
    end do
    if (count == 0) exit
  end do

  ptr = 1                                           ! uniq
  do while ( ptr < n )
    if (wrk(ptr,1) == wrk(ptr+1,1) .and. wrk(ptr,2) == wrk(ptr+1,2)) then
      do i=ptr+1,n-1; wrk(i,:) = wrk(i+1,:); end do
      n = n - 1
    else
      ptr = ptr + 1
    end if
  end do

  diag = 0; diag(1:n,:) = wrk(1:n,:)
end subroutine ovlp_uniqsort


!    Name : ovlp_blockdecomp
!         :
! Purpose : This subroutine introduces another partition to the matrix in order
!         : to apply block algorithms to the overlapped block-diagonal region.
!         :
!   Input : - ``blks'' specifies the size of a block in the new partition.
!         : - ``n'' is the size of the whole matrix.
!         : - diag(i,1) is the row/column number of the head of the i-th
!         :   diagonal block, and diag(i,2) is the row/column number of
!         :   the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is nonzero.
!         :
!  Output : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from succeeding computations.
!         :
!         : - ``nblk'' is the number of partitions in the row/column direction.

pure subroutine ovlp_blockdcmp(blks,n,diag,subd,head,tail,nonz,nblk)
  integer, intent(in ) :: blks, n, diag(:,:), subd(:)
  integer, intent(out) :: head(:), tail(:), nonz(:,:), nblk

  integer :: h(n), t(n), blksize
  integer :: i, j, k, m, acts

  ! For a small matrix, a new partition is not introduced.

  blksize = max(4, blks)
  if (n <= blksize) then
    head(1) = 1; tail(1) = n; nonz(1,1) = 1; nblk = 1
    return
  end if

  ! A new block partition is defined so that a two-by-two diagonal
  ! block is not separated.

  nblk = n / blksize 
  if (mod(n,blksize) /= 0) nblk = nblk + 1
  acts = n / nblk
    
  m = 0
  do i=1,nblk-1
    h(i) = m + 1 
    m = m + acts - subd(m + acts)
    t(i) = m
  end do
  h(nblk) = m + 1; t(nblk) = n

  ! The block mask, ``nonz'', is defined.

  head = 0; tail = 0; nonz = 0
  head(1:nblk) = h(1:nblk); tail(1:nblk) = t(1:nblk)

  nonz = 0
  do i=1,nblk
    do j=1,nblk
      do k=1,size(diag,1)
        if ( diag(k,2) >= head(i) .and. &
             diag(k,1) <= tail(i) .and. &
             diag(k,2) >= head(j) .and. &
             diag(k,1) <= tail(j)       ) then
          nonz(i,j) = 1; exit
        end if
      end do
    end do
  end do
end subroutine ovlp_blockdcmp


!    Name : rs_ovlp_times
!         :
! Purpose : This function computes a block-wise product of two overlapped
!         : block-diagonal matrices without repeating equivalent operations.
!         : Two relevant matrices must have the same size and the same block
!         : structure.
!         :
!   Input : - The square arrays, ``A'' and ``B'', are two overlapped
!         :   block-diagonal matrices. Both matrices must have the same size
!         :   and the same block structure.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The result is the product of A and B, which is equivalent
!         : to the block-wise product, ( A_{ii}B_{ii}; 1<=i<=size(diag,1) ).

pure function rs_ovlp_times(diag,subd,head,tail,nonz,a,b)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  real   (kind=sp), intent(in) :: a(:,:), b(:,:)
  real   (kind=sp)             :: rs_ovlp_times(size(a,1),size(b,2))

  integer :: i, j, k, hi, ti, hj, tj, hk, tk

  rs_ovlp_times = 0.0_sp

  do i=1,size(nonz,1)
    hi = head(i); ti = tail(i)
    do j=i,size(nonz,1)
      if (nonz(i,j) == 0) cycle
      hj = head(j); tj = tail(j)
      do k=i,j
        if (nonz(i,k) == 0 .or. nonz(k,j) == 0) cycle
        hk = head(k); tk = tail(k)
        rs_ovlp_times(hi:ti,hj:tj) = &
        rs_ovlp_times(hi:ti,hj:tj) + (a(hi:ti,hk:tk) .sqtimes. b(hk:tk,hj:tj))
      end do
    end do
  end do

  ! Clear unnecessary entries.

  call rs_ovlp_clear(diag,subd,head,tail,nonz,rs_ovlp_times)
end function rs_ovlp_times


!    Name : rs_ovlp_clear
!         :
! Purpose : This subroutine clears all entries outside the coverage of
!         : the overlapped diagonal blocks.
!         :
!   Input : - The square array ``A'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : ``A'' stores the result, where all entries are zero outside the
!         : the overlapped block-diagonal region.

pure subroutine rs_ovlp_clear(diag,subd,head,tail,nonz,a)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  real   (kind=sp), intent(inout) :: a(:,:)

  integer :: i, j, m, n, band(size(a,1))

  i = head(1); i = tail(1); i = nonz(1,1)     ! to avoid compiler's warning

  m = size(diag,1); n = size(a,1)
  do i=1,m
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do
  do j=1,n-2
    if (subd(j) == 0) a(j+1,j) = 0.0_sp
    a(j+2:n,j) = 0.0_sp
  end do
  if (1 < n) then
    if (subd(n-1) == 0) a(n,n-1) = 0.0_sp
  end if

  do i=1,n-1
    a(i,band(i)+1:n) = 0.0_sp
  end do
end subroutine rs_ovlp_clear


!    Name : rs_ludcmp
! Purpose : This is a formal subroutine to perform the LU factorization.
!    Note : Nothing is done in this function.

pure subroutine rs_ovlp_ludcmp(diag,subd,head,tail,nonz,a,p)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  real   (kind=sp), intent(inout) :: a(:,:)
  integer,          intent(out  ) :: p(:)
  integer :: i

  p(1) = diag(1,1); p(1) = subd(1)                ! to avoid compilers warning
  p(1) = head(1); p(1) = tail(1); p(1) = nonz(1,1)
  p(1) = int(a(1,1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine rs_ovlp_ludcmp


!    Name : rs_ovlp_ainvbf
!         :
! Purpose : This function computes the block-wise product of an inverse matrix
!         : and a matrix without repeating equivalent operations. The relevant
!         : two matrices are overlapped block-diagonal matrices that have the
!         : same size and the same structure.
!         :
!   Input : - ``A'' is an overlapped block-diagonal matrix.
!         : - ``p'' is an array of integer, which is not used.
!         : - ``B'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The result is the product of A^{-1} and B, which is equivalent to
!         : the block-wise product, ( A_{ii}^{-1}B_{ii}; 1<=i<=size(diag,1) ).

pure function rs_ovlp_ainvbf(diag,subd,head,tail,nonz,a,p,b)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  real   (kind=sp), intent(in) :: a(:,:), b(:,:)
  integer,          intent(in) :: p(:)
  real   (kind=sp)             :: rs_ovlp_ainvbf(size(a,1),size(b,2))

  call rs_ovlp_ainvb(diag,subd,head,tail,nonz,a,p,b,rs_ovlp_ainvbf)
end function rs_ovlp_ainvbf


!    Name : rs_ovlp_ainvb
!         :
! Purpose : This subroutine computes the block-wise product of an inverse matrix
!         : and a matrix without repeating equivalent operations. The relevant
!         : two matrices are overlapped block-diagonal matrices that have the
!         : same size and the same block structure.
!         :
!   Input : - ``A'' is an overlapped block-diagonal matrix.
!         : - ``p'' is an array of integer, which is not used.
!         : - ``B'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The matrix ``V'' stores the product of A^{-1} and B, which is
!         : equivalent to the block-wise product,
!         :    ( A_{ii}^{-1}B_{ii}; 1<=i<=size(diag,1) ).

pure subroutine rs_ovlp_ainvb(diag,subd,head,tail,nonz,a,p,b,v)
  integer,          intent(in ) :: diag(:,:), subd(:)
  integer,          intent(in ) :: head(:), tail(:), nonz(:,:)
  real   (kind=sp), intent(in ) :: a(:,:), b(:,:)
  integer,          intent(in ) :: p(:)
  real   (kind=sp), intent(out) :: v(:,:)

  real   (kind=sp) :: w(1:size(a,1),1:size(a,2))
  integer          :: i, j, k, n, hi, ti, hj, tj, hk, tk

  n = size(nonz,1)                                 ! the number of blocks
  v = 0.0_sp; w = 0.0_sp                           ! solution and working area

  do i=n,1,-1
    hi = head(i); ti = tail(i)
    do j=i,n
      if (nonz(i,j) == 0) cycle
      hj = head(j); tj = tail(j)
      do k=i+1,j
        if (nonz(i,k) == 0 .or. nonz(k,j) == 0) cycle
        hk = head(k); tk = tail(k)
        v(hi:ti,hj:tj) = &
        v(hi:ti,hj:tj) - (a(hi:ti,hk:tk) .sqtimes. v(hk:tk,hj:tj))
      end do
      v(hi:ti,hj:tj) = v(hi:ti,hj:tj) + b(hi:ti,hj:tj)
      call rstr_ainvbblock(                                                   &
        a(hi:ti,hi:ti), p(hi:ti), v(hi:ti,hj:tj), w(hi:ti,hj:tj)              &
      )
      v(hi:ti,hj:tj) = w(hi:ti,hj:tj)
    end do
  end do

  ! Clear unnecessary entries.

  call rs_ovlp_clear(diag,subd,head,tail,nonz,v)   ! clear off-diagonal entries
end subroutine rs_ovlp_ainvb


!    Name : cs_ovlp_times
!         :
! Purpose : This function computes a block-wise product of two overlapped
!         : block-diagonal matrices without repeating equivalent operations.
!         : Two relevant matrices must have the same size and the same block
!         : structure.
!         :
!   Input : - The square arrays, ``A'' and ``B'', are two overlapped
!         :   block-diagonal matrices. Both matrices must have the same size
!         :   and the same block structure.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The result is the product of A and B, which is equivalent
!         : to the block-wise product, ( A_{ii}B_{ii}; 1<=i<=size(diag,1) ).

pure function cs_ovlp_times(diag,subd,head,tail,nonz,a,b)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  complex(kind=sp), intent(in) :: a(:,:), b(:,:)
  complex(kind=sp)             :: cs_ovlp_times(size(a,1),size(b,2))

  integer :: i, j, k, hi, ti, hj, tj, hk, tk

  cs_ovlp_times = cmplx(0.0_sp,0.0_sp,sp)

  do i=1,size(nonz,1)
    hi = head(i); ti = tail(i)
    do j=i,size(nonz,1)
      if (nonz(i,j) == 0) cycle
      hj = head(j); tj = tail(j)
      do k=i,j
        if (nonz(i,k) == 0 .or. nonz(k,j) == 0) cycle
        hk = head(k); tk = tail(k)
        cs_ovlp_times(hi:ti,hj:tj) = &
        cs_ovlp_times(hi:ti,hj:tj) + (a(hi:ti,hk:tk) .sqtimes. b(hk:tk,hj:tj))
      end do
    end do
  end do

  ! Clear unnecessary entries.

  call cs_ovlp_clear(diag,subd,head,tail,nonz,cs_ovlp_times)
end function cs_ovlp_times


!    Name : cs_ovlp_clear
!         :
! Purpose : This subroutine clears all entries outside the coverage of
!         : the overlapped diagonal blocks.
!         :
!   Input : - The square array ``A'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : ``A'' stores the result, where all entries are zero outside the
!         : the overlapped block-diagonal region.

pure subroutine cs_ovlp_clear(diag,subd,head,tail,nonz,a)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  complex(kind=sp), intent(inout) :: a(:,:)

  integer :: i, j, m, n, band(size(a,1))

  i = head(1); i = tail(1); i = nonz(1,1)     ! to avoid compiler's warning

  m = size(diag,1); n = size(a,1)
  do i=1,m
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do
  do j=1,n-2
    if (subd(j) == 0) a(j+1,j) = cmplx(0.0_sp,0.0_sp,sp)
    a(j+2:n,j) = cmplx(0.0_sp,0.0_sp,sp)
  end do
  if (1 < n) then
    if (subd(n-1) == 0) a(n,n-1) = cmplx(0.0_sp,0.0_sp,sp)
  end if

  do i=1,n-1
    a(i,band(i)+1:n) = cmplx(0.0_sp,0.0_sp,sp)
  end do
end subroutine cs_ovlp_clear


!    Name : cs_ludcmp
! Purpose : This is a formal subroutine to perform the LU factorization.
!    Note : Nothing is done in this function.

pure subroutine cs_ovlp_ludcmp(diag,subd,head,tail,nonz,a,p)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  complex(kind=sp), intent(inout) :: a(:,:)
  integer,          intent(out  ) :: p(:)
  integer :: i

  p(1) = diag(1,1); p(1) = subd(1)                ! to avoid compilers warning
  p(1) = head(1); p(1) = tail(1); p(1) = nonz(1,1)
  p(1) = int(a(1,1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine cs_ovlp_ludcmp


!    Name : cs_ovlp_ainvbf
!         :
! Purpose : This function computes the block-wise product of an inverse matrix
!         : and a matrix without repeating equivalent operations. The relevant
!         : two matrices are overlapped block-diagonal matrices that have the
!         : same size and the same structure.
!         :
!   Input : - ``A'' is an overlapped block-diagonal matrix.
!         : - ``p'' is an array of integer, which is not used.
!         : - ``B'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The result is the product of A^{-1} and B, which is equivalent to
!         : the block-wise product, ( A_{ii}^{-1}B_{ii}; 1<=i<=size(diag,1) ).

pure function cs_ovlp_ainvbf(diag,subd,head,tail,nonz,a,p,b)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  complex(kind=sp), intent(in) :: a(:,:), b(:,:)
  integer,          intent(in) :: p(:)
  complex(kind=sp)             :: cs_ovlp_ainvbf(size(a,1),size(b,2))

  call cs_ovlp_ainvb(diag,subd,head,tail,nonz,a,p,b,cs_ovlp_ainvbf)
end function cs_ovlp_ainvbf


!    Name : cs_ovlp_ainvb
!         :
! Purpose : This subroutine computes the block-wise product of an inverse matrix
!         : and a matrix without repeating equivalent operations. The relevant
!         : two matrices are overlapped block-diagonal matrices that have the
!         : same size and the same block structure.
!         :
!   Input : - ``A'' is an overlapped block-diagonal matrix.
!         : - ``p'' is an array of integer, which is not used.
!         : - ``B'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The matrix ``V'' stores the product of A^{-1} and B, which is
!         : equivalent to the block-wise product,
!         :    ( A_{ii}^{-1}B_{ii}; 1<=i<=size(diag,1) ).

pure subroutine cs_ovlp_ainvb(diag,subd,head,tail,nonz,a,p,b,v)
  integer,          intent(in ) :: diag(:,:), subd(:)
  integer,          intent(in ) :: head(:), tail(:), nonz(:,:)
  complex(kind=sp), intent(in ) :: a(:,:), b(:,:)
  integer,          intent(in ) :: p(:)
  complex(kind=sp), intent(out) :: v(:,:)

  complex(kind=sp) :: w(1:size(a,1),1:size(a,2))
  integer          :: i, j, k, n, hi, ti, hj, tj, hk, tk

  n = size(nonz,1)                                 ! the number of blocks
  v = cmplx(0.0_sp,0.0_sp,sp); w = cmplx(0.0_sp,0.0_sp,sp)                           ! solution and working area

  do i=n,1,-1
    hi = head(i); ti = tail(i)
    do j=i,n
      if (nonz(i,j) == 0) cycle
      hj = head(j); tj = tail(j)
      do k=i+1,j
        if (nonz(i,k) == 0 .or. nonz(k,j) == 0) cycle
        hk = head(k); tk = tail(k)
        v(hi:ti,hj:tj) = &
        v(hi:ti,hj:tj) - (a(hi:ti,hk:tk) .sqtimes. v(hk:tk,hj:tj))
      end do
      v(hi:ti,hj:tj) = v(hi:ti,hj:tj) + b(hi:ti,hj:tj)
      call cstr_ainvbblock(                                                   &
        a(hi:ti,hi:ti), p(hi:ti), v(hi:ti,hj:tj), w(hi:ti,hj:tj)              &
      )
      v(hi:ti,hj:tj) = w(hi:ti,hj:tj)
    end do
  end do

  ! Clear unnecessary entries.

  call cs_ovlp_clear(diag,subd,head,tail,nonz,v)   ! clear off-diagonal entries
end subroutine cs_ovlp_ainvb


!    Name : rw_ovlp_times
!         :
! Purpose : This function computes a block-wise product of two overlapped
!         : block-diagonal matrices without repeating equivalent operations.
!         : Two relevant matrices must have the same size and the same block
!         : structure.
!         :
!   Input : - The square arrays, ``A'' and ``B'', are two overlapped
!         :   block-diagonal matrices. Both matrices must have the same size
!         :   and the same block structure.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The result is the product of A and B, which is equivalent
!         : to the block-wise product, ( A_{ii}B_{ii}; 1<=i<=size(diag,1) ).

pure function rw_ovlp_times(diag,subd,head,tail,nonz,a,b)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  real   (kind=wp), intent(in) :: a(:,:), b(:,:)
  real   (kind=wp)             :: rw_ovlp_times(size(a,1),size(b,2))

  integer :: i, j, k, hi, ti, hj, tj, hk, tk

  rw_ovlp_times = 0.0_wp

  do i=1,size(nonz,1)
    hi = head(i); ti = tail(i)
    do j=i,size(nonz,1)
      if (nonz(i,j) == 0) cycle
      hj = head(j); tj = tail(j)
      do k=i,j
        if (nonz(i,k) == 0 .or. nonz(k,j) == 0) cycle
        hk = head(k); tk = tail(k)
        rw_ovlp_times(hi:ti,hj:tj) = &
        rw_ovlp_times(hi:ti,hj:tj) + (a(hi:ti,hk:tk) .sqtimes. b(hk:tk,hj:tj))
      end do
    end do
  end do

  ! Clear unnecessary entries.

  call rw_ovlp_clear(diag,subd,head,tail,nonz,rw_ovlp_times)
end function rw_ovlp_times


!    Name : rw_ovlp_clear
!         :
! Purpose : This subroutine clears all entries outside the coverage of
!         : the overlapped diagonal blocks.
!         :
!   Input : - The square array ``A'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : ``A'' stores the result, where all entries are zero outside the
!         : the overlapped block-diagonal region.

pure subroutine rw_ovlp_clear(diag,subd,head,tail,nonz,a)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  real   (kind=wp), intent(inout) :: a(:,:)

  integer :: i, j, m, n, band(size(a,1))

  i = head(1); i = tail(1); i = nonz(1,1)     ! to avoid compiler's warning

  m = size(diag,1); n = size(a,1)
  do i=1,m
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do
  do j=1,n-2
    if (subd(j) == 0) a(j+1,j) = 0.0_wp
    a(j+2:n,j) = 0.0_wp
  end do
  if (1 < n) then
    if (subd(n-1) == 0) a(n,n-1) = 0.0_wp
  end if

  do i=1,n-1
    a(i,band(i)+1:n) = 0.0_wp
  end do
end subroutine rw_ovlp_clear


!    Name : rw_ludcmp
! Purpose : This is a formal subroutine to perform the LU factorization.
!    Note : Nothing is done in this function.

pure subroutine rw_ovlp_ludcmp(diag,subd,head,tail,nonz,a,p)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  real   (kind=wp), intent(inout) :: a(:,:)
  integer,          intent(out  ) :: p(:)
  integer :: i

  p(1) = diag(1,1); p(1) = subd(1)                ! to avoid compilers warning
  p(1) = head(1); p(1) = tail(1); p(1) = nonz(1,1)
  p(1) = int(a(1,1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine rw_ovlp_ludcmp


!    Name : rw_ovlp_ainvbf
!         :
! Purpose : This function computes the block-wise product of an inverse matrix
!         : and a matrix without repeating equivalent operations. The relevant
!         : two matrices are overlapped block-diagonal matrices that have the
!         : same size and the same structure.
!         :
!   Input : - ``A'' is an overlapped block-diagonal matrix.
!         : - ``p'' is an array of integer, which is not used.
!         : - ``B'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The result is the product of A^{-1} and B, which is equivalent to
!         : the block-wise product, ( A_{ii}^{-1}B_{ii}; 1<=i<=size(diag,1) ).

pure function rw_ovlp_ainvbf(diag,subd,head,tail,nonz,a,p,b)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  real   (kind=wp), intent(in) :: a(:,:), b(:,:)
  integer,          intent(in) :: p(:)
  real   (kind=wp)             :: rw_ovlp_ainvbf(size(a,1),size(b,2))

  call rw_ovlp_ainvb(diag,subd,head,tail,nonz,a,p,b,rw_ovlp_ainvbf)
end function rw_ovlp_ainvbf


!    Name : rw_ovlp_ainvb
!         :
! Purpose : This subroutine computes the block-wise product of an inverse matrix
!         : and a matrix without repeating equivalent operations. The relevant
!         : two matrices are overlapped block-diagonal matrices that have the
!         : same size and the same block structure.
!         :
!   Input : - ``A'' is an overlapped block-diagonal matrix.
!         : - ``p'' is an array of integer, which is not used.
!         : - ``B'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The matrix ``V'' stores the product of A^{-1} and B, which is
!         : equivalent to the block-wise product,
!         :    ( A_{ii}^{-1}B_{ii}; 1<=i<=size(diag,1) ).

pure subroutine rw_ovlp_ainvb(diag,subd,head,tail,nonz,a,p,b,v)
  integer,          intent(in ) :: diag(:,:), subd(:)
  integer,          intent(in ) :: head(:), tail(:), nonz(:,:)
  real   (kind=wp), intent(in ) :: a(:,:), b(:,:)
  integer,          intent(in ) :: p(:)
  real   (kind=wp), intent(out) :: v(:,:)

  real   (kind=wp) :: w(1:size(a,1),1:size(a,2))
  integer          :: i, j, k, n, hi, ti, hj, tj, hk, tk

  n = size(nonz,1)                                 ! the number of blocks
  v = 0.0_wp; w = 0.0_wp                           ! solution and working area

  do i=n,1,-1
    hi = head(i); ti = tail(i)
    do j=i,n
      if (nonz(i,j) == 0) cycle
      hj = head(j); tj = tail(j)
      do k=i+1,j
        if (nonz(i,k) == 0 .or. nonz(k,j) == 0) cycle
        hk = head(k); tk = tail(k)
        v(hi:ti,hj:tj) = &
        v(hi:ti,hj:tj) - (a(hi:ti,hk:tk) .sqtimes. v(hk:tk,hj:tj))
      end do
      v(hi:ti,hj:tj) = v(hi:ti,hj:tj) + b(hi:ti,hj:tj)
      call rwtr_ainvbblock(                                                   &
        a(hi:ti,hi:ti), p(hi:ti), v(hi:ti,hj:tj), w(hi:ti,hj:tj)              &
      )
      v(hi:ti,hj:tj) = w(hi:ti,hj:tj)
    end do
  end do

  ! Clear unnecessary entries.

  call rw_ovlp_clear(diag,subd,head,tail,nonz,v)   ! clear off-diagonal entries
end subroutine rw_ovlp_ainvb


!    Name : cw_ovlp_times
!         :
! Purpose : This function computes a block-wise product of two overlapped
!         : block-diagonal matrices without repeating equivalent operations.
!         : Two relevant matrices must have the same size and the same block
!         : structure.
!         :
!   Input : - The square arrays, ``A'' and ``B'', are two overlapped
!         :   block-diagonal matrices. Both matrices must have the same size
!         :   and the same block structure.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The result is the product of A and B, which is equivalent
!         : to the block-wise product, ( A_{ii}B_{ii}; 1<=i<=size(diag,1) ).

pure function cw_ovlp_times(diag,subd,head,tail,nonz,a,b)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  complex(kind=wp), intent(in) :: a(:,:), b(:,:)
  complex(kind=wp)             :: cw_ovlp_times(size(a,1),size(b,2))

  integer :: i, j, k, hi, ti, hj, tj, hk, tk

  cw_ovlp_times = cmplx(0.0_wp,0.0_wp,wp)

  do i=1,size(nonz,1)
    hi = head(i); ti = tail(i)
    do j=i,size(nonz,1)
      if (nonz(i,j) == 0) cycle
      hj = head(j); tj = tail(j)
      do k=i,j
        if (nonz(i,k) == 0 .or. nonz(k,j) == 0) cycle
        hk = head(k); tk = tail(k)
        cw_ovlp_times(hi:ti,hj:tj) = &
        cw_ovlp_times(hi:ti,hj:tj) + (a(hi:ti,hk:tk) .sqtimes. b(hk:tk,hj:tj))
      end do
    end do
  end do

  ! Clear unnecessary entries.

  call cw_ovlp_clear(diag,subd,head,tail,nonz,cw_ovlp_times)
end function cw_ovlp_times


!    Name : cw_ovlp_clear
!         :
! Purpose : This subroutine clears all entries outside the coverage of
!         : the overlapped diagonal blocks.
!         :
!   Input : - The square array ``A'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : ``A'' stores the result, where all entries are zero outside the
!         : the overlapped block-diagonal region.

pure subroutine cw_ovlp_clear(diag,subd,head,tail,nonz,a)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  complex(kind=wp), intent(inout) :: a(:,:)

  integer :: i, j, m, n, band(size(a,1))

  i = head(1); i = tail(1); i = nonz(1,1)     ! to avoid compiler's warning

  m = size(diag,1); n = size(a,1)
  do i=1,m
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do
  do j=1,n-2
    if (subd(j) == 0) a(j+1,j) = cmplx(0.0_wp,0.0_wp,wp)
    a(j+2:n,j) = cmplx(0.0_wp,0.0_wp,wp)
  end do
  if (1 < n) then
    if (subd(n-1) == 0) a(n,n-1) = cmplx(0.0_wp,0.0_wp,wp)
  end if

  do i=1,n-1
    a(i,band(i)+1:n) = cmplx(0.0_wp,0.0_wp,wp)
  end do
end subroutine cw_ovlp_clear


!    Name : cw_ludcmp
! Purpose : This is a formal subroutine to perform the LU factorization.
!    Note : Nothing is done in this function.

pure subroutine cw_ovlp_ludcmp(diag,subd,head,tail,nonz,a,p)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  complex(kind=wp), intent(inout) :: a(:,:)
  integer,          intent(out  ) :: p(:)
  integer :: i

  p(1) = diag(1,1); p(1) = subd(1)                ! to avoid compilers warning
  p(1) = head(1); p(1) = tail(1); p(1) = nonz(1,1)
  p(1) = int(a(1,1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine cw_ovlp_ludcmp


!    Name : cw_ovlp_ainvbf
!         :
! Purpose : This function computes the block-wise product of an inverse matrix
!         : and a matrix without repeating equivalent operations. The relevant
!         : two matrices are overlapped block-diagonal matrices that have the
!         : same size and the same structure.
!         :
!   Input : - ``A'' is an overlapped block-diagonal matrix.
!         : - ``p'' is an array of integer, which is not used.
!         : - ``B'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The result is the product of A^{-1} and B, which is equivalent to
!         : the block-wise product, ( A_{ii}^{-1}B_{ii}; 1<=i<=size(diag,1) ).

pure function cw_ovlp_ainvbf(diag,subd,head,tail,nonz,a,p,b)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  complex(kind=wp), intent(in) :: a(:,:), b(:,:)
  integer,          intent(in) :: p(:)
  complex(kind=wp)             :: cw_ovlp_ainvbf(size(a,1),size(b,2))

  call cw_ovlp_ainvb(diag,subd,head,tail,nonz,a,p,b,cw_ovlp_ainvbf)
end function cw_ovlp_ainvbf


!    Name : cw_ovlp_ainvb
!         :
! Purpose : This subroutine computes the block-wise product of an inverse matrix
!         : and a matrix without repeating equivalent operations. The relevant
!         : two matrices are overlapped block-diagonal matrices that have the
!         : same size and the same block structure.
!         :
!   Input : - ``A'' is an overlapped block-diagonal matrix.
!         : - ``p'' is an array of integer, which is not used.
!         : - ``B'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The matrix ``V'' stores the product of A^{-1} and B, which is
!         : equivalent to the block-wise product,
!         :    ( A_{ii}^{-1}B_{ii}; 1<=i<=size(diag,1) ).

pure subroutine cw_ovlp_ainvb(diag,subd,head,tail,nonz,a,p,b,v)
  integer,          intent(in ) :: diag(:,:), subd(:)
  integer,          intent(in ) :: head(:), tail(:), nonz(:,:)
  complex(kind=wp), intent(in ) :: a(:,:), b(:,:)
  integer,          intent(in ) :: p(:)
  complex(kind=wp), intent(out) :: v(:,:)

  complex(kind=wp) :: w(1:size(a,1),1:size(a,2))
  integer          :: i, j, k, n, hi, ti, hj, tj, hk, tk

  n = size(nonz,1)                                 ! the number of blocks
  v = cmplx(0.0_wp,0.0_wp,wp); w = cmplx(0.0_wp,0.0_wp,wp)                           ! solution and working area

  do i=n,1,-1
    hi = head(i); ti = tail(i)
    do j=i,n
      if (nonz(i,j) == 0) cycle
      hj = head(j); tj = tail(j)
      do k=i+1,j
        if (nonz(i,k) == 0 .or. nonz(k,j) == 0) cycle
        hk = head(k); tk = tail(k)
        v(hi:ti,hj:tj) = &
        v(hi:ti,hj:tj) - (a(hi:ti,hk:tk) .sqtimes. v(hk:tk,hj:tj))
      end do
      v(hi:ti,hj:tj) = v(hi:ti,hj:tj) + b(hi:ti,hj:tj)
      call cwtr_ainvbblock(                                                   &
        a(hi:ti,hi:ti), p(hi:ti), v(hi:ti,hj:tj), w(hi:ti,hj:tj)              &
      )
      v(hi:ti,hj:tj) = w(hi:ti,hj:tj)
    end do
  end do

  ! Clear unnecessary entries.

  call cw_ovlp_clear(diag,subd,head,tail,nonz,v)   ! clear off-diagonal entries
end subroutine cw_ovlp_ainvb

#ifdef __USE_TPREC

!    Name : rt_ovlp_times
!         :
! Purpose : This function computes a block-wise product of two overlapped
!         : block-diagonal matrices without repeating equivalent operations.
!         : Two relevant matrices must have the same size and the same block
!         : structure.
!         :
!   Input : - The square arrays, ``A'' and ``B'', are two overlapped
!         :   block-diagonal matrices. Both matrices must have the same size
!         :   and the same block structure.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The result is the product of A and B, which is equivalent
!         : to the block-wise product, ( A_{ii}B_{ii}; 1<=i<=size(diag,1) ).

pure function rt_ovlp_times(diag,subd,head,tail,nonz,a,b)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  real   (kind=tp), intent(in) :: a(:,:), b(:,:)
  real   (kind=tp)             :: rt_ovlp_times(size(a,1),size(b,2))

  integer :: i, j, k, hi, ti, hj, tj, hk, tk

  rt_ovlp_times = 0.0_tp

  do i=1,size(nonz,1)
    hi = head(i); ti = tail(i)
    do j=i,size(nonz,1)
      if (nonz(i,j) == 0) cycle
      hj = head(j); tj = tail(j)
      do k=i,j
        if (nonz(i,k) == 0 .or. nonz(k,j) == 0) cycle
        hk = head(k); tk = tail(k)
        rt_ovlp_times(hi:ti,hj:tj) = &
        rt_ovlp_times(hi:ti,hj:tj) + (a(hi:ti,hk:tk) .sqtimes. b(hk:tk,hj:tj))
      end do
    end do
  end do

  ! Clear unnecessary entries.

  call rt_ovlp_clear(diag,subd,head,tail,nonz,rt_ovlp_times)
end function rt_ovlp_times


!    Name : rt_ovlp_clear
!         :
! Purpose : This subroutine clears all entries outside the coverage of
!         : the overlapped diagonal blocks.
!         :
!   Input : - The square array ``A'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : ``A'' stores the result, where all entries are zero outside the
!         : the overlapped block-diagonal region.

pure subroutine rt_ovlp_clear(diag,subd,head,tail,nonz,a)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  real   (kind=tp), intent(inout) :: a(:,:)

  integer :: i, j, m, n, band(size(a,1))

  i = head(1); i = tail(1); i = nonz(1,1)     ! to avoid compiler's warning

  m = size(diag,1); n = size(a,1)
  do i=1,m
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do
  do j=1,n-2
    if (subd(j) == 0) a(j+1,j) = 0.0_tp
    a(j+2:n,j) = 0.0_tp
  end do
  if (1 < n) then
    if (subd(n-1) == 0) a(n,n-1) = 0.0_tp
  end if

  do i=1,n-1
    a(i,band(i)+1:n) = 0.0_tp
  end do
end subroutine rt_ovlp_clear


!    Name : rt_ludcmp
! Purpose : This is a formal subroutine to perform the LU factorization.
!    Note : Nothing is done in this function.

pure subroutine rt_ovlp_ludcmp(diag,subd,head,tail,nonz,a,p)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  real   (kind=tp), intent(inout) :: a(:,:)
  integer,          intent(out  ) :: p(:)
  integer :: i

  p(1) = diag(1,1); p(1) = subd(1)                ! to avoid compilers warning
  p(1) = head(1); p(1) = tail(1); p(1) = nonz(1,1)
  p(1) = int(a(1,1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine rt_ovlp_ludcmp


!    Name : rt_ovlp_ainvbf
!         :
! Purpose : This function computes the block-wise product of an inverse matrix
!         : and a matrix without repeating equivalent operations. The relevant
!         : two matrices are overlapped block-diagonal matrices that have the
!         : same size and the same structure.
!         :
!   Input : - ``A'' is an overlapped block-diagonal matrix.
!         : - ``p'' is an array of integer, which is not used.
!         : - ``B'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The result is the product of A^{-1} and B, which is equivalent to
!         : the block-wise product, ( A_{ii}^{-1}B_{ii}; 1<=i<=size(diag,1) ).

pure function rt_ovlp_ainvbf(diag,subd,head,tail,nonz,a,p,b)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  real   (kind=tp), intent(in) :: a(:,:), b(:,:)
  integer,          intent(in) :: p(:)
  real   (kind=tp)             :: rt_ovlp_ainvbf(size(a,1),size(b,2))

  call rt_ovlp_ainvb(diag,subd,head,tail,nonz,a,p,b,rt_ovlp_ainvbf)
end function rt_ovlp_ainvbf


!    Name : rt_ovlp_ainvb
!         :
! Purpose : This subroutine computes the block-wise product of an inverse matrix
!         : and a matrix without repeating equivalent operations. The relevant
!         : two matrices are overlapped block-diagonal matrices that have the
!         : same size and the same block structure.
!         :
!   Input : - ``A'' is an overlapped block-diagonal matrix.
!         : - ``p'' is an array of integer, which is not used.
!         : - ``B'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The matrix ``V'' stores the product of A^{-1} and B, which is
!         : equivalent to the block-wise product,
!         :    ( A_{ii}^{-1}B_{ii}; 1<=i<=size(diag,1) ).

pure subroutine rt_ovlp_ainvb(diag,subd,head,tail,nonz,a,p,b,v)
  integer,          intent(in ) :: diag(:,:), subd(:)
  integer,          intent(in ) :: head(:), tail(:), nonz(:,:)
  real   (kind=tp), intent(in ) :: a(:,:), b(:,:)
  integer,          intent(in ) :: p(:)
  real   (kind=tp), intent(out) :: v(:,:)

  real   (kind=tp) :: w(1:size(a,1),1:size(a,2))
  integer          :: i, j, k, n, hi, ti, hj, tj, hk, tk

  n = size(nonz,1)                                 ! the number of blocks
  v = 0.0_tp; w = 0.0_tp                           ! solution and working area

  do i=n,1,-1
    hi = head(i); ti = tail(i)
    do j=i,n
      if (nonz(i,j) == 0) cycle
      hj = head(j); tj = tail(j)
      do k=i+1,j
        if (nonz(i,k) == 0 .or. nonz(k,j) == 0) cycle
        hk = head(k); tk = tail(k)
        v(hi:ti,hj:tj) = &
        v(hi:ti,hj:tj) - (a(hi:ti,hk:tk) .sqtimes. v(hk:tk,hj:tj))
      end do
      v(hi:ti,hj:tj) = v(hi:ti,hj:tj) + b(hi:ti,hj:tj)
      call rttr_ainvbblock(                                                   &
        a(hi:ti,hi:ti), p(hi:ti), v(hi:ti,hj:tj), w(hi:ti,hj:tj)              &
      )
      v(hi:ti,hj:tj) = w(hi:ti,hj:tj)
    end do
  end do

  ! Clear unnecessary entries.

  call rt_ovlp_clear(diag,subd,head,tail,nonz,v)   ! clear off-diagonal entries
end subroutine rt_ovlp_ainvb


!    Name : ct_ovlp_times
!         :
! Purpose : This function computes a block-wise product of two overlapped
!         : block-diagonal matrices without repeating equivalent operations.
!         : Two relevant matrices must have the same size and the same block
!         : structure.
!         :
!   Input : - The square arrays, ``A'' and ``B'', are two overlapped
!         :   block-diagonal matrices. Both matrices must have the same size
!         :   and the same block structure.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The result is the product of A and B, which is equivalent
!         : to the block-wise product, ( A_{ii}B_{ii}; 1<=i<=size(diag,1) ).

pure function ct_ovlp_times(diag,subd,head,tail,nonz,a,b)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  complex(kind=tp), intent(in) :: a(:,:), b(:,:)
  complex(kind=tp)             :: ct_ovlp_times(size(a,1),size(b,2))

  integer :: i, j, k, hi, ti, hj, tj, hk, tk

  ct_ovlp_times = cmplx(0.0_tp,0.0_tp,tp)

  do i=1,size(nonz,1)
    hi = head(i); ti = tail(i)
    do j=i,size(nonz,1)
      if (nonz(i,j) == 0) cycle
      hj = head(j); tj = tail(j)
      do k=i,j
        if (nonz(i,k) == 0 .or. nonz(k,j) == 0) cycle
        hk = head(k); tk = tail(k)
        ct_ovlp_times(hi:ti,hj:tj) = &
        ct_ovlp_times(hi:ti,hj:tj) + (a(hi:ti,hk:tk) .sqtimes. b(hk:tk,hj:tj))
      end do
    end do
  end do

  ! Clear unnecessary entries.

  call ct_ovlp_clear(diag,subd,head,tail,nonz,ct_ovlp_times)
end function ct_ovlp_times


!    Name : ct_ovlp_clear
!         :
! Purpose : This subroutine clears all entries outside the coverage of
!         : the overlapped diagonal blocks.
!         :
!   Input : - The square array ``A'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : ``A'' stores the result, where all entries are zero outside the
!         : the overlapped block-diagonal region.

pure subroutine ct_ovlp_clear(diag,subd,head,tail,nonz,a)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  complex(kind=tp), intent(inout) :: a(:,:)

  integer :: i, j, m, n, band(size(a,1))

  i = head(1); i = tail(1); i = nonz(1,1)     ! to avoid compiler's warning

  m = size(diag,1); n = size(a,1)
  do i=1,m
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do
  do j=1,n-2
    if (subd(j) == 0) a(j+1,j) = cmplx(0.0_tp,0.0_tp,tp)
    a(j+2:n,j) = cmplx(0.0_tp,0.0_tp,tp)
  end do
  if (1 < n) then
    if (subd(n-1) == 0) a(n,n-1) = cmplx(0.0_tp,0.0_tp,tp)
  end if

  do i=1,n-1
    a(i,band(i)+1:n) = cmplx(0.0_tp,0.0_tp,tp)
  end do
end subroutine ct_ovlp_clear


!    Name : ct_ludcmp
! Purpose : This is a formal subroutine to perform the LU factorization.
!    Note : Nothing is done in this function.

pure subroutine ct_ovlp_ludcmp(diag,subd,head,tail,nonz,a,p)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  complex(kind=tp), intent(inout) :: a(:,:)
  integer,          intent(out  ) :: p(:)
  integer :: i

  p(1) = diag(1,1); p(1) = subd(1)                ! to avoid compilers warning
  p(1) = head(1); p(1) = tail(1); p(1) = nonz(1,1)
  p(1) = int(a(1,1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine ct_ovlp_ludcmp


!    Name : ct_ovlp_ainvbf
!         :
! Purpose : This function computes the block-wise product of an inverse matrix
!         : and a matrix without repeating equivalent operations. The relevant
!         : two matrices are overlapped block-diagonal matrices that have the
!         : same size and the same structure.
!         :
!   Input : - ``A'' is an overlapped block-diagonal matrix.
!         : - ``p'' is an array of integer, which is not used.
!         : - ``B'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The result is the product of A^{-1} and B, which is equivalent to
!         : the block-wise product, ( A_{ii}^{-1}B_{ii}; 1<=i<=size(diag,1) ).

pure function ct_ovlp_ainvbf(diag,subd,head,tail,nonz,a,p,b)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  complex(kind=tp), intent(in) :: a(:,:), b(:,:)
  integer,          intent(in) :: p(:)
  complex(kind=tp)             :: ct_ovlp_ainvbf(size(a,1),size(b,2))

  call ct_ovlp_ainvb(diag,subd,head,tail,nonz,a,p,b,ct_ovlp_ainvbf)
end function ct_ovlp_ainvbf


!    Name : ct_ovlp_ainvb
!         :
! Purpose : This subroutine computes the block-wise product of an inverse matrix
!         : and a matrix without repeating equivalent operations. The relevant
!         : two matrices are overlapped block-diagonal matrices that have the
!         : same size and the same block structure.
!         :
!   Input : - ``A'' is an overlapped block-diagonal matrix.
!         : - ``p'' is an array of integer, which is not used.
!         : - ``B'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The matrix ``V'' stores the product of A^{-1} and B, which is
!         : equivalent to the block-wise product,
!         :    ( A_{ii}^{-1}B_{ii}; 1<=i<=size(diag,1) ).

pure subroutine ct_ovlp_ainvb(diag,subd,head,tail,nonz,a,p,b,v)
  integer,          intent(in ) :: diag(:,:), subd(:)
  integer,          intent(in ) :: head(:), tail(:), nonz(:,:)
  complex(kind=tp), intent(in ) :: a(:,:), b(:,:)
  integer,          intent(in ) :: p(:)
  complex(kind=tp), intent(out) :: v(:,:)

  complex(kind=tp) :: w(1:size(a,1),1:size(a,2))
  integer          :: i, j, k, n, hi, ti, hj, tj, hk, tk

  n = size(nonz,1)                                 ! the number of blocks
  v = cmplx(0.0_tp,0.0_tp,tp); w = cmplx(0.0_tp,0.0_tp,tp)                           ! solution and working area

  do i=n,1,-1
    hi = head(i); ti = tail(i)
    do j=i,n
      if (nonz(i,j) == 0) cycle
      hj = head(j); tj = tail(j)
      do k=i+1,j
        if (nonz(i,k) == 0 .or. nonz(k,j) == 0) cycle
        hk = head(k); tk = tail(k)
        v(hi:ti,hj:tj) = &
        v(hi:ti,hj:tj) - (a(hi:ti,hk:tk) .sqtimes. v(hk:tk,hj:tj))
      end do
      v(hi:ti,hj:tj) = v(hi:ti,hj:tj) + b(hi:ti,hj:tj)
      call cttr_ainvbblock(                                                   &
        a(hi:ti,hi:ti), p(hi:ti), v(hi:ti,hj:tj), w(hi:ti,hj:tj)              &
      )
      v(hi:ti,hj:tj) = w(hi:ti,hj:tj)
    end do
  end do

  ! Clear unnecessary entries.

  call ct_ovlp_clear(diag,subd,head,tail,nonz,v)   ! clear off-diagonal entries
end subroutine ct_ovlp_ainvb

#endif
#ifdef __USE_QPREC

!    Name : rq_ovlp_times
!         :
! Purpose : This function computes a block-wise product of two overlapped
!         : block-diagonal matrices without repeating equivalent operations.
!         : Two relevant matrices must have the same size and the same block
!         : structure.
!         :
!   Input : - The square arrays, ``A'' and ``B'', are two overlapped
!         :   block-diagonal matrices. Both matrices must have the same size
!         :   and the same block structure.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The result is the product of A and B, which is equivalent
!         : to the block-wise product, ( A_{ii}B_{ii}; 1<=i<=size(diag,1) ).

pure function rq_ovlp_times(diag,subd,head,tail,nonz,a,b)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  real   (kind=qp), intent(in) :: a(:,:), b(:,:)
  real   (kind=qp)             :: rq_ovlp_times(size(a,1),size(b,2))

  integer :: i, j, k, hi, ti, hj, tj, hk, tk

  rq_ovlp_times = 0.0_qp

  do i=1,size(nonz,1)
    hi = head(i); ti = tail(i)
    do j=i,size(nonz,1)
      if (nonz(i,j) == 0) cycle
      hj = head(j); tj = tail(j)
      do k=i,j
        if (nonz(i,k) == 0 .or. nonz(k,j) == 0) cycle
        hk = head(k); tk = tail(k)
        rq_ovlp_times(hi:ti,hj:tj) = &
        rq_ovlp_times(hi:ti,hj:tj) + (a(hi:ti,hk:tk) .sqtimes. b(hk:tk,hj:tj))
      end do
    end do
  end do

  ! Clear unnecessary entries.

  call rq_ovlp_clear(diag,subd,head,tail,nonz,rq_ovlp_times)
end function rq_ovlp_times


!    Name : rq_ovlp_clear
!         :
! Purpose : This subroutine clears all entries outside the coverage of
!         : the overlapped diagonal blocks.
!         :
!   Input : - The square array ``A'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : ``A'' stores the result, where all entries are zero outside the
!         : the overlapped block-diagonal region.

pure subroutine rq_ovlp_clear(diag,subd,head,tail,nonz,a)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  real   (kind=qp), intent(inout) :: a(:,:)

  integer :: i, j, m, n, band(size(a,1))

  i = head(1); i = tail(1); i = nonz(1,1)     ! to avoid compiler's warning

  m = size(diag,1); n = size(a,1)
  do i=1,m
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do
  do j=1,n-2
    if (subd(j) == 0) a(j+1,j) = 0.0_qp
    a(j+2:n,j) = 0.0_qp
  end do
  if (1 < n) then
    if (subd(n-1) == 0) a(n,n-1) = 0.0_qp
  end if

  do i=1,n-1
    a(i,band(i)+1:n) = 0.0_qp
  end do
end subroutine rq_ovlp_clear


!    Name : rq_ludcmp
! Purpose : This is a formal subroutine to perform the LU factorization.
!    Note : Nothing is done in this function.

pure subroutine rq_ovlp_ludcmp(diag,subd,head,tail,nonz,a,p)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  real   (kind=qp), intent(inout) :: a(:,:)
  integer,          intent(out  ) :: p(:)
  integer :: i

  p(1) = diag(1,1); p(1) = subd(1)                ! to avoid compilers warning
  p(1) = head(1); p(1) = tail(1); p(1) = nonz(1,1)
  p(1) = int(a(1,1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine rq_ovlp_ludcmp


!    Name : rq_ovlp_ainvbf
!         :
! Purpose : This function computes the block-wise product of an inverse matrix
!         : and a matrix without repeating equivalent operations. The relevant
!         : two matrices are overlapped block-diagonal matrices that have the
!         : same size and the same structure.
!         :
!   Input : - ``A'' is an overlapped block-diagonal matrix.
!         : - ``p'' is an array of integer, which is not used.
!         : - ``B'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The result is the product of A^{-1} and B, which is equivalent to
!         : the block-wise product, ( A_{ii}^{-1}B_{ii}; 1<=i<=size(diag,1) ).

pure function rq_ovlp_ainvbf(diag,subd,head,tail,nonz,a,p,b)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  real   (kind=qp), intent(in) :: a(:,:), b(:,:)
  integer,          intent(in) :: p(:)
  real   (kind=qp)             :: rq_ovlp_ainvbf(size(a,1),size(b,2))

  call rq_ovlp_ainvb(diag,subd,head,tail,nonz,a,p,b,rq_ovlp_ainvbf)
end function rq_ovlp_ainvbf


!    Name : rq_ovlp_ainvb
!         :
! Purpose : This subroutine computes the block-wise product of an inverse matrix
!         : and a matrix without repeating equivalent operations. The relevant
!         : two matrices are overlapped block-diagonal matrices that have the
!         : same size and the same block structure.
!         :
!   Input : - ``A'' is an overlapped block-diagonal matrix.
!         : - ``p'' is an array of integer, which is not used.
!         : - ``B'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The matrix ``V'' stores the product of A^{-1} and B, which is
!         : equivalent to the block-wise product,
!         :    ( A_{ii}^{-1}B_{ii}; 1<=i<=size(diag,1) ).

pure subroutine rq_ovlp_ainvb(diag,subd,head,tail,nonz,a,p,b,v)
  integer,          intent(in ) :: diag(:,:), subd(:)
  integer,          intent(in ) :: head(:), tail(:), nonz(:,:)
  real   (kind=qp), intent(in ) :: a(:,:), b(:,:)
  integer,          intent(in ) :: p(:)
  real   (kind=qp), intent(out) :: v(:,:)

  real   (kind=qp) :: w(1:size(a,1),1:size(a,2))
  integer          :: i, j, k, n, hi, ti, hj, tj, hk, tk

  n = size(nonz,1)                                 ! the number of blocks
  v = 0.0_qp; w = 0.0_qp                           ! solution and working area

  do i=n,1,-1
    hi = head(i); ti = tail(i)
    do j=i,n
      if (nonz(i,j) == 0) cycle
      hj = head(j); tj = tail(j)
      do k=i+1,j
        if (nonz(i,k) == 0 .or. nonz(k,j) == 0) cycle
        hk = head(k); tk = tail(k)
        v(hi:ti,hj:tj) = &
        v(hi:ti,hj:tj) - (a(hi:ti,hk:tk) .sqtimes. v(hk:tk,hj:tj))
      end do
      v(hi:ti,hj:tj) = v(hi:ti,hj:tj) + b(hi:ti,hj:tj)
      call rqtr_ainvbblock(                                                   &
        a(hi:ti,hi:ti), p(hi:ti), v(hi:ti,hj:tj), w(hi:ti,hj:tj)              &
      )
      v(hi:ti,hj:tj) = w(hi:ti,hj:tj)
    end do
  end do

  ! Clear unnecessary entries.

  call rq_ovlp_clear(diag,subd,head,tail,nonz,v)   ! clear off-diagonal entries
end subroutine rq_ovlp_ainvb


!    Name : cq_ovlp_times
!         :
! Purpose : This function computes a block-wise product of two overlapped
!         : block-diagonal matrices without repeating equivalent operations.
!         : Two relevant matrices must have the same size and the same block
!         : structure.
!         :
!   Input : - The square arrays, ``A'' and ``B'', are two overlapped
!         :   block-diagonal matrices. Both matrices must have the same size
!         :   and the same block structure.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The result is the product of A and B, which is equivalent
!         : to the block-wise product, ( A_{ii}B_{ii}; 1<=i<=size(diag,1) ).

pure function cq_ovlp_times(diag,subd,head,tail,nonz,a,b)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  complex(kind=qp), intent(in) :: a(:,:), b(:,:)
  complex(kind=qp)             :: cq_ovlp_times(size(a,1),size(b,2))

  integer :: i, j, k, hi, ti, hj, tj, hk, tk

  cq_ovlp_times = cmplx(0.0_qp,0.0_qp,qp)

  do i=1,size(nonz,1)
    hi = head(i); ti = tail(i)
    do j=i,size(nonz,1)
      if (nonz(i,j) == 0) cycle
      hj = head(j); tj = tail(j)
      do k=i,j
        if (nonz(i,k) == 0 .or. nonz(k,j) == 0) cycle
        hk = head(k); tk = tail(k)
        cq_ovlp_times(hi:ti,hj:tj) = &
        cq_ovlp_times(hi:ti,hj:tj) + (a(hi:ti,hk:tk) .sqtimes. b(hk:tk,hj:tj))
      end do
    end do
  end do

  ! Clear unnecessary entries.

  call cq_ovlp_clear(diag,subd,head,tail,nonz,cq_ovlp_times)
end function cq_ovlp_times


!    Name : cq_ovlp_clear
!         :
! Purpose : This subroutine clears all entries outside the coverage of
!         : the overlapped diagonal blocks.
!         :
!   Input : - The square array ``A'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry is not zero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : ``A'' stores the result, where all entries are zero outside the
!         : the overlapped block-diagonal region.

pure subroutine cq_ovlp_clear(diag,subd,head,tail,nonz,a)
  integer,          intent(in   ) :: diag(:,:), subd(:)
  integer,          intent(in   ) :: head(:), tail(:), nonz(:,:)
  complex(kind=qp), intent(inout) :: a(:,:)

  integer :: i, j, m, n, band(size(a,1))

  i = head(1); i = tail(1); i = nonz(1,1)     ! to avoid compiler's warning

  m = size(diag,1); n = size(a,1)
  do i=1,m
    band(diag(i,1):diag(i,2)) = diag(i,2)
  end do
  do j=1,n-2
    if (subd(j) == 0) a(j+1,j) = cmplx(0.0_qp,0.0_qp,qp)
    a(j+2:n,j) = cmplx(0.0_qp,0.0_qp,qp)
  end do
  if (1 < n) then
    if (subd(n-1) == 0) a(n,n-1) = cmplx(0.0_qp,0.0_qp,qp)
  end if

  do i=1,n-1
    a(i,band(i)+1:n) = cmplx(0.0_qp,0.0_qp,qp)
  end do
end subroutine cq_ovlp_clear


!    Name : cq_ludcmp
! Purpose : This is a formal subroutine to perform the LU factorization.
!    Note : Nothing is done in this function.

pure subroutine cq_ovlp_ludcmp(diag,subd,head,tail,nonz,a,p)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  complex(kind=qp), intent(inout) :: a(:,:)
  integer,          intent(out  ) :: p(:)
  integer :: i

  p(1) = diag(1,1); p(1) = subd(1)                ! to avoid compilers warning
  p(1) = head(1); p(1) = tail(1); p(1) = nonz(1,1)
  p(1) = int(a(1,1)*0)
  p = (/(i,i=1,size(p))/)
end subroutine cq_ovlp_ludcmp


!    Name : cq_ovlp_ainvbf
!         :
! Purpose : This function computes the block-wise product of an inverse matrix
!         : and a matrix without repeating equivalent operations. The relevant
!         : two matrices are overlapped block-diagonal matrices that have the
!         : same size and the same structure.
!         :
!   Input : - ``A'' is an overlapped block-diagonal matrix.
!         : - ``p'' is an array of integer, which is not used.
!         : - ``B'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The result is the product of A^{-1} and B, which is equivalent to
!         : the block-wise product, ( A_{ii}^{-1}B_{ii}; 1<=i<=size(diag,1) ).

pure function cq_ovlp_ainvbf(diag,subd,head,tail,nonz,a,p,b)
  integer,          intent(in) :: diag(:,:), subd(:)
  integer,          intent(in) :: head(:), tail(:), nonz(:,:)
  complex(kind=qp), intent(in) :: a(:,:), b(:,:)
  integer,          intent(in) :: p(:)
  complex(kind=qp)             :: cq_ovlp_ainvbf(size(a,1),size(b,2))

  call cq_ovlp_ainvb(diag,subd,head,tail,nonz,a,p,b,cq_ovlp_ainvbf)
end function cq_ovlp_ainvbf


!    Name : cq_ovlp_ainvb
!         :
! Purpose : This subroutine computes the block-wise product of an inverse matrix
!         : and a matrix without repeating equivalent operations. The relevant
!         : two matrices are overlapped block-diagonal matrices that have the
!         : same size and the same block structure.
!         :
!   Input : - ``A'' is an overlapped block-diagonal matrix.
!         : - ``p'' is an array of integer, which is not used.
!         : - ``B'' is an overlapped block-diagonal matrix.
!         :
!         : - diag(i,1) stores the row/column number of the head
!         :   of the i-th diagonal block, and diag(i,2) stores the
!         :   row/column number of the tail of the i-th diagonal block.
!         : - When subd(i)==0, the (i+1,i)-th entry of A and B is zero, while
!         :   if   subd(i)==1, the (i+1,i)-th entry of A and B is nonzero.
!         :
!         : The following three arguments, head, tail, and nonz, define
!         : another partition for the matrix. This partition is introduced
!         : to apply a block algorithm.
!         :
!         : - head(i) is the row/column number of the head of i-th partition.
!         : - tail(i) is the row/column number of the tail of i-th partition.
!         : - When nonz(i,j)==0, the (i,j)-th partition is out of the coverage
!         :   of the diagonal blocks defined by diag(:,:). When nonz(i,j)==1,
!         :   the intersection of the diagonal blocks and (i,j)-th partition
!         :   is non-empty. In the former case, the (i,j)-th partition is
!         :   omitted from the computation.
!         :
!  Output : The matrix ``V'' stores the product of A^{-1} and B, which is
!         : equivalent to the block-wise product,
!         :    ( A_{ii}^{-1}B_{ii}; 1<=i<=size(diag,1) ).

pure subroutine cq_ovlp_ainvb(diag,subd,head,tail,nonz,a,p,b,v)
  integer,          intent(in ) :: diag(:,:), subd(:)
  integer,          intent(in ) :: head(:), tail(:), nonz(:,:)
  complex(kind=qp), intent(in ) :: a(:,:), b(:,:)
  integer,          intent(in ) :: p(:)
  complex(kind=qp), intent(out) :: v(:,:)

  complex(kind=qp) :: w(1:size(a,1),1:size(a,2))
  integer          :: i, j, k, n, hi, ti, hj, tj, hk, tk

  n = size(nonz,1)                                 ! the number of blocks
  v = cmplx(0.0_qp,0.0_qp,qp); w = cmplx(0.0_qp,0.0_qp,qp)                           ! solution and working area

  do i=n,1,-1
    hi = head(i); ti = tail(i)
    do j=i,n
      if (nonz(i,j) == 0) cycle
      hj = head(j); tj = tail(j)
      do k=i+1,j
        if (nonz(i,k) == 0 .or. nonz(k,j) == 0) cycle
        hk = head(k); tk = tail(k)
        v(hi:ti,hj:tj) = &
        v(hi:ti,hj:tj) - (a(hi:ti,hk:tk) .sqtimes. v(hk:tk,hj:tj))
      end do
      v(hi:ti,hj:tj) = v(hi:ti,hj:tj) + b(hi:ti,hj:tj)
      call cqtr_ainvbblock(                                                   &
        a(hi:ti,hi:ti), p(hi:ti), v(hi:ti,hj:tj), w(hi:ti,hj:tj)              &
      )
      v(hi:ti,hj:tj) = w(hi:ti,hj:tj)
    end do
  end do

  ! Clear unnecessary entries.

  call cq_ovlp_clear(diag,subd,head,tail,nonz,v)   ! clear off-diagonal entries
end subroutine cq_ovlp_ainvb

#endif

end module ovlpmatrix


! Eigentest is a package that produces real test matrices with known
! eigensystems.  A test matrix, called an eigenmat, is generated in a
! factored form, in which the user can specify the eigenvalues
! (including complex conjugate eigenvalues and Jordan blocks) and has
! some control over the condition of the eigenvalues and eigenvectors.
! An eigenmat $A$ of order $n$ requires only $O(n)$ storage for its
! representation.  Auxiliary programs permit the computation of
! 
!  (A - s*I)*b, (A - s*I)'*b, inv(A - s*I)*b, and inv(A -s*I)'*b
! 
! in $O(n)$ operations.  Thus eigenmats are suitable for testing
! algorithms based on Krylov sequences, as well as others based on
! matrix vector products.  A special routine computes specified
! eigenvectors of an eigenmat and the condition of its eigenvalue.
! 
! For more details see
!
!    1. Eigentest: A Test Maatrix Generator for Large Scale
!       Eigenproblem
!
!    2. Eigentest: The Fortran 95 User's Guide
!
! Pdf files are included with this distribution.


! Contents
! 
!    subroutine EigenmatAlloc
!       Allocates storage for an eigenmat.
! 
!    subroutine EigenmatFree
!       Frees the storage for an eigenmat.
! 
!    subroutine EigenmatProd
!       Computes the above products.
! 
!    subroutine HsvdProd
!       Computes products of special matrices in the factorization.
!    
!    subroutine Hprod
!       Computes products with a Householder transformation.
! 
!    subroutine EigenmatVecs
!       Computes eigenvectors of an eigenmat.
! 
! Coded by Roger Lee and Pete Stewart


module eigentest

implicit none

integer, parameter :: wp = kind(0.0d0)

!  An hsvd matrix has the form
!
!      X = diag(X1, X2, ..., Xnblocks),
!
!  where
!
!      Xi = (I - ui*ui^T)*Si*(I - vi*vi^T).
!
!  with Sk diagonal having positive diagonal entries.  The vectors
!  uk, vk, and the diagonals of Sk are stacked in order in the arrays
!  X%u, X%v, and X%sig.  abs(X%bs(i)) points to the beginning of the ith
!  block and abs(X.bs(X.nblocks-1)) = n+1.  If X.bs(i+1)<0,
!  the ith block is an identity matrix.

   type hsvdmat
      integer               :: n       ! The order of the matrix.
      integer               :: nblocks ! The number of blocks in the hsvdmat.
      integer, allocatable  :: bs(:)   ! abs(bs(i)) is the index of the
                                       ! start of the i-th block.
                                       ! abs(bs(nblocks+1)) = n+1.
                                       ! If bs(i+1)<0, the i-th block
                                       ! is an identity.
      real(wp), allocatable :: u(:)    ! The vectors generating the left
                                       ! Householder transformations.
      real(wp), allocatable :: v(:)    ! The vectors generating the right
                                       ! Householder transformations.
      real(wp), allocatable :: sig(:)  ! The diagonals of the Si.

   end type hsvdmat


!  An eigenmat has the form
!
!      A = Y*Z*(L - shift*I)*Z^-1*Y^-1.
!
!  Here L is a block triangular with 1x1 blocks containing the
!  real eigenvalues of A and 2x2 blocks conaining the complex
!  eignevalues of A in the form
!
!               |  mu  nu |
!               |         |.
!               | -nu  mu |
!
!  Y and Z are hsvd transformations, with Y having only one block

   type eigenmat
      integer               :: n       ! The order of the matrix.
      real(wp), allocatable :: eig(:)  ! Array containing the
                                       ! eigenvalues of A
                                       ! or the supediagonals
                                       ! of a Jordan block.
      integer, allocatable  :: type(:) ! The type of the entry
                                       ! in eig (see below).
      type(hsvdmat)          :: Y, Z   ! The outer and inner hsvd
                                       ! transformations.
   end type eigenmat

!  The contents of eig and type are determined as follows.
!
!  Real eigenvalue
!     type(i) = 1    eig(i) = the eigenvalue.
!
!  Complex conjugate pair
!     type(i) = 2    eig(i) = the real part.
!     type(i) = 3    eig(i) = the complex part.
!
!  Jordan block of order k
!     type(i) = -k     eig(i)   = the eigenvalue.
!     type(i+j) = -1   eig(i+j) = the jth superdiagonal
!                i = 1,...,k-1.

contains

   subroutine EigenmatAlloc(A, n, nblocks, yident, zident)
      type(eigenmat), intent(inout) :: A
      integer, intent(in)           :: n
      integer, intent(in)           :: nblocks
      logical                       :: yident, zident

!     EigenmatAlloc allocates memory for an eigenmat A.  In addition,
!     it initializes A%n, A%Z%n, A%Z%nblocks, A%Z%bs(1), A%Z%bs(nblocks+1),
!     A%Y%n, A%Y%nblocks, A%Y%bs.
!
!     A        The eigenmat to be initialized.
!     n        The order of A.
!     nblocks  The number of blocks in the hsvd matrix Z.
!     yident   If yident is .true., Y is declared to be an identity
!              matrix and only Y%nblock and Y%bs are initialized.
!     zident   If zident is .true., Z is declared to be an identity
!              matrix and only Z%nblock and Z%bs are initialized.
!
!     EigenmatAlloc allocates memory for the arrays in A, A%Y, and A%Z.
!     In addition it initializes A%n, A%Z%n, A%Z%nblocks, A%Z%bs(1),
!     A%Z%bs(nblocks+1), A%Y%n, A%y%nblocks, and A%Y%bs.  All other
!     arrays are initialized to zero.


      A%n = n

      ! Allocate the eigenvalues.

      allocate(A%eig(n))
      allocate(A%type(n))

      ! Allocoate the outer hsvdmat Y.

      A%Y%n = n
      A%Y%nblocks = 1
      allocate(A%Y%bs(2))
      if (yident) then ! Y is an identity.
         A%Y%bs(1) = 1
         A%Y%bs(2) = -(n+1)
      else
         A%Y%bs(1) = 1;
         A%Y%bs(2) = n+1;
         allocate(A%Y%u(n))
         allocate(A%Y%v(n))
         allocate(A%Y%sig(n))
      end if

      ! Allocate the inner hsvdmat Z.

      A%Z%n = n
      if (zident) then ! Z is an identity.
         allocate(A%Z%bs(2))
         A%Z%nblocks = 1
         A%Z%bs(1) = 1
         A%Z%bs(2) = -(n+1)
      else
         allocate(A%Z%bs(nblocks+1))
         A%Z%nblocks = nblocks
         A%Z%bs(1) = 1
         A%Z%bs(nblocks+1) = n+1
         allocate(A%Z%u(n))
         allocate(A%Z%v(n))
         allocate(A%Z%sig(n))
      end if

   end subroutine EigenmatAlloc

   subroutine EigenmatFree(A)
      type(Eigenmat), intent(inout) :: A

!     EigenmatFree deallocates the storage of the eigenmat A.

      if (allocated(A%eig))   deallocate(A%eig)
      if (allocated(A%type))  deallocate(A%type)
      if (allocated(A%Y%bs))  deallocate(A%Y%bs)
      if (allocated(A%Y%u))   deallocate(A%Y%u)
      if (allocated(A%Y%v))   deallocate(A%Y%v)
      if (allocated(A%Y%sig)) deallocate(A%Y%sig)
      if (allocated(A%Z%bs))  deallocate(A%Z%bs)
      if (allocated(A%Z%u))   deallocate(A%Z%u)
      if (allocated(A%Z%v))   deallocate(A%Z%v)
      if (allocated(A%Z%sig)) deallocate(A%Z%sig)

   end subroutine EigenmatFree

   subroutine EigenmatProd(A, ncols, B, ldb, C, ldc, shift, job)
      type(Eigenmat), intent(in) :: A
      integer, intent(in)        :: ncols, ldb, ldc
      real(wp), intent(in)       :: B(ldb,*)
      real(wp), intent(inout)    :: C(ldc,*)
      real(wp), intent(in)       :: shift
      character(*), intent(in)   :: job

!     EigmatProd computes product, perhaps shifted and inverted, of the
!     eigenmat A and a matrix B, placing the result in C.
!
!     A       The eigenmat
!     ncols   Number of columns in the matrix B.
!     B       The array containing the matrix B.
!     ldb     The leading dimension of B.
!     C       The array array containing the matrix C.
!     ldc     The leading dimension of C.
!     shift   A shift.
!     job     A string specifying the operation to be performed.
!
!             "ab"    C = (A - shift*I)*B
!             "atb"   C = (A - shift*I)^T*B
!             "aib"   C = (A - shift*I)^-1*B
!             "aitb"  C = (A - shift*I)^-T*B

      integer   :: i, n, nc, j, endj
      real(wp)  :: mu, nu, mult, u11, u12, u22
      real(wp)  :: temp(1, ncols)

      

      n = A%n    ! Get the order of the matrix.
      nc = ncols ! An abbreviation for ncols.

      ! Copy B into C.

      C(1:n, 1:ncols) = B(1:n, 1:ncols) 

      select case (job)

      case ('ab')

         ! Compute C = A*B  = (Y*Z*L*Z^-1*Y^-1)*C.

         ! Compute Z^-1*Y^-1*C.

         call HsvdProd(A%Y, ncols, C, ldc, 'aib') 
         call HsvdProd(A%Z, ncols, C, ldc, 'aib')

         ! Compute L*C.  The index i points to successive blocks
         ! of L.

         i = 1
         do while (i<=n)

            if (A%type(i) .eq. 1) then

               ! 1x1 block.

               C(i,1:nc) = (A%eig(i)-shift)*C(i,1:nc)
               i = i + 1

            elseif (i .eq. n) then

               stop "Error in EigenmatProd: 2x2 block starts at bs(n)."

            elseif ((A%type(i) .eq. 2) .and. (A%type(i+1) .eq. 3)) then

               ! 2x2 block.

               mu = A%eig(i)-shift
               nu = A%eig(i+1)
               temp(1,:) = mu*C(i,1:nc) + nu*C(i+1,1:nc)
               C(i+1,1:nc) = -nu*C(i,1:nc) + mu*C(i+1,1:nc)
               C(i,1:nc) = temp(1,:)
               i = i + 2

            elseif (A%type(i) .lt. -1) then
               
               ! A Jordan block
               
               endj = i+abs(A%type(i))-1
               
               if (endj .gt. A%n) then
                  stop "Error in EigematProd: Jordan block is too large"
               end if

               if (any(A%type(i+1:endj).ne.-1)) then
                  stop "Error in EigematProd: Illegal type."
               end if
               
               mu = A%eig(i)-shift

               do j = i, endj-1
                  C(j,1:nc) = mu*C(j,1:nc) + A%eig(j+1)*C(j+1,1:nc)
               end do
               C(endj,1:nc) = mu*C(endj,1:nc)
               i = endj + 1

            else
               stop "Error in EigematProd: Illegal type."
            end if

         end do

         ! Compute Y*Z*C.

         call HsvdProd(A%Z, ncols, C, ldc, 'ab')
         call Hsvdprod(A%Y, ncols, C, ldc, 'ab')

      case ('atb')
         
         ! Compute C = A^T*B = (Y^-T*Z^-T*L^T*Z^T*Y^T)*C.

         ! Compute Z^T*Y^T*C.

         call HsvdProd(A%Y, ncols, C, ldc, 'atb')
         call Hsvdprod(A%Z, ncols, C, ldc, 'atb')

         ! Compute L*C.  The index i points to successive blocks
         ! of L.

         i = 1
         do while(i<=n)

            ! 1x1 block.

            if (A%type(i) .eq. 1) then
               C(i,1:nc) = (A%eig(i)-shift)*C(i,1:nc)
               i = i + 1

            elseif (i .eq. n) then

               stop "Error in EigenmatProd: 2x2 block starts at bs(n)."

            elseif ((A%type(i) .eq. 2) .and. (A%type(i+1) .eq. 3)) then

               ! 2x2 block.

               mu = A%eig(i)-shift
               nu = -A%eig(i+1)

               temp(1,:) = mu*C(i,1:nc) + nu*C(i+1,1:nc)
               C(i+1,1:nc) = -nu*C(i,1:nc) + mu*C(i+1,1:nc)
               C(i,1:nc) = temp(1,:)
               i = i + 2                  

            elseif (A%type(i) .lt. -1) then
               
               ! A Jordan block
               
               endj = i+abs(A%type(i))-1

               if (endj .gt. A%n) then
                  stop "Error in EigematProd: Jordan block is too large."
               end if

               if (any(A%type(i+1:endj).ne.-1)) then
                  stop "Error in EigematProd: Illegal type."
               end if
               
               mu = A%eig(i)-shift

               do j = endj, i+1, -1
                  C(j,1:nc) = mu*C(j,1:nc) + A%eig(j)*C(j-1,1:nc)
               end do
               C(i,1:nc) = mu*C(i,1:nc)

               i = endj + 1
            
            else
               stop "Error in Eigematprod: Illegal type."
            end if
         end do

         ! Compute C = Y^-T*Z.^-T*C

         call HsvdProd(A%Z, ncols, C, ldc, 'aitb')
         call HsvdProd(A%Y, ncols, C, ldc, 'aitb')  

      case ('aib')

         ! Compute C = A^-1*C = (Y*Z*L^-1*Y^-1*Z^-1)*C.

         ! Compute Z^-1*Y^-1*C.

         call HsvdProd(A%Y, ncols, C, ldc, 'aib')  
         call HsvdProd(A%Z, ncols, C, ldc, 'aib')

         ! Compute L^-1*C.  The index i points to successive blocks
         ! of L.
         i = 1
         do while (i<=A%n)
            if (A%type(i) .eq. 1) then

               ! 1x1 block.

               mu = 1.0_wp/(A%eig(i)-shift)
               C(i,1:nc) = mu*C(i,1:nc)
               i = i + 1

            elseif (i .eq. n) then

               stop "Error in EigenmatProd: 2x2 block starts at bs(n)."

            elseif ((A%type(i) .eq. 2) .and. (A%type(i+1) .eq. 3)) then

               ! 2x2 block.  Use Gaussian elimination with
               ! partial pivoting to apply its inverse.

               mu = A%eig(i)-shift
               nu = A%eig(i+1)

               if (abs(mu) .ge. abs(nu)) then

                  ! Do not pivot.

                  u11 = mu
                  u12 = nu
                  mult = -nu/mu
                  u22 = mu - mult*nu

                  C(i+1,1:nc) = (C(i+1,1:nc) - mult*C(i,1:nc))/u22;
                  C(i,1:nc)   = (C(i,1:nc) - u12*C(i+1,1:nc))/u11;

               else

                  ! Pivot.

                  u11 = -nu
                  u12 = mu
                  mult = -mu/nu
                  u22 = nu - mult*mu

                  temp(1,:) = C(i+1,1:nc)
                  C(i+1,1:nc) = (C(i,1:nc) - mult*temp(1,:))/u22;
                  C(i,1:nc)   = (temp(1,:) - u12*C(i+1,1:nc))/u11;

               end if

               i = i + 2

            elseif (A%type(i) .lt. -1) then
               
               ! A Jordan block
               
               endj = i+abs(A%type(i))-1
               
               if (endj .gt. A%n) then
                  stop "Error in EigematProd: Jordan block too large."
               end if
               
               if (any(A%type(i+1:endj).ne.-1)) then
                  stop "Error in EigematProd: Illegal type."
               end if
               
               mu = 1.0D0/(A%eig(i)-shift)

               C(endj,1:nc) = C(endj,1:nc)*mu
               do j=endj-1, i, -1
                  C(j,1:nc) = (C(j,1:nc)-(C(j+1,1:nc)*A%eig(j+1)))*mu
               end do

               i = endj + 1

            else
               stop "Error in Eigematprod: Illegal type."
            end if
         end do

         ! Compute C = Y*Z*C.

         call HsvdProd(A%Z, ncols, C, ldc, 'ab')
         call HsvdProd(A%Y, ncols, C, ldc, 'ab')

      case ('aitb')

         ! Compute C = A^-T*C = (Y^-T*Z^-T*L^-T*Z^T*Y^T)*C.

         ! Compute C = Z^T*Y^T*C.

         call HsvdProd(A%Y, ncols, C, ldc, 'atb' )
         call HsvdProd(A%Z, ncols, C, ldc, 'atb')

         ! Compute L^-T*C.  The index i points to successive blocks
         ! of L.

         i = 1
         do while(i<=A%n)
            if (A%type(i) .eq. 1) then

               ! 1x1 block.

               mu = 1.0_wp/(A%eig(i)-shift)
               C(i,1:nc) = mu*c(i,1:nc)
               i = i + 1

            elseif (i .eq. n) then

               stop "Error in EigenmatProd: 2x2 block starts at bs(n)."

            elseif ((A%type(i) .eq. 2) .and. (A%type(i+1) .eq. 3)) then

               ! 2x2 block.  Use Gaussian eliminaton with
               ! partial pivoting to apply its inverse transpose.

               mu = A%eig(i)-shift
               nu = -A%eig(i+1)

               if (abs(mu) .ge. abs(nu)) then

                  ! Do not pivot.

                  u11 = mu
                  u12 = nu
                  mult = -nu/mu
                  u22 = mu - mult*nu

                  C(i+1,1:nc) = (C(i+1,1:nc) - mult*C(i,1:nc))/u22;
                  C(i,1:nc)   = (C(i,1:nc) - u12*C(i+1,1:nc))/u11;

               else

                  ! Pivot

                  u11 = -nu
                  u12 = mu
                  mult = -mu/nu
                  u22 = nu - mult*mu

                  temp(1,:) = C(i+1,1:nc)
                  C(i+1,1:nc) = (C(i,1:nc) - mult*temp(1,:))/u22;
                  C(i,1:nc)   = (temp(1,:) - u12*c(i+1,1:nc))/u11;

               end if

               i = i + 2

            elseif (A%type(i) .lt. -1) then
               
               ! A Jordan block
               
               endj = i+abs(A%type(i))-1
               
               if (endj .gt. A%n) then
                  stop "Error in EigematProd: Jordan block too large."
               end if

               if (any(A%type(i+1:endj) .ne. -1)) then
                  stop "Error in EigematProd: Illegal type."
               end if
               
               mu = 1.0D0/(A%eig(i)-shift)

               C(i,1:nc) = C(i,1:nc)*mu
               do j = i+1,endj
                  C(j,1:nc) = (C(j,1:nc)-C(j-1,1:nc)*A%eig(j))*mu
               end do

               i = endj + 1

            else
               stop "Error in EigematProd: Illegal type."
            end if
         end do

        ! Compute C = Y^-T*Z^-T*C.

         call HsvdProd(A%Z, ncols, C, ldc, 'aitb')
         call HsvdProd(A%Y, ncols, C, ldc, 'aitb')

      case default

         stop "Error in EigementProd: Illegal operation."

      end select

   end subroutine EigenmatProd


   subroutine HsvdProd(X, ncols, B, ldb, job)
      type(hsvdmat), intent(in) :: X
      integer, intent(in)       :: ncols, ldb
      real(wp), intent(inout)   :: B(ldb, *)
      character(*), intent(in)  :: job

!     HsvdProd computes the product of computes the product of an hsvdmat
!     X and a matrix B, overwriting B with the product.
! 
!     X       Pointer to the hsvdmat.
!     ncols   The number of columns in B.
!     B       The array B.
!     ldb     The leading dimension of B.
!     job     A string specifying the operation to be performed.
! 
!             "ab"   B <- X*B
!             "atb"  B <- X^T*B
!             "aib"  B <- X^-1*B
!             "aitb" B <- X^-T*B

      integer :: i, i1, i2, k, nblocks, nc

      nblocks = X%nblocks
      nc = ncols           ! Abbreviation for ncols.

      select case (job)
      case ('ab')

         ! Compute B= X*B .

         ! Loop over the blocks.

         do k = 1, nblocks 

            if (X%bs(k+1) > 0) then
               i1 = abs(X%bs(k))
               i2 = abs(X%bs(k+1))-1

               ! B = (I - vk*vk^T)*B.

               call Hprod(i1, i2, X%v, nc, B, ldb)

               ! B = Sk*B

               do i = i1, i2
                  B(i,1:nc) = X%sig(i)*B(i,1:nc)
               end do

               ! B = (I - uk*uk^T)*B.

               call Hprod(i1, i2, X%u, nc, B, ldb)
            end if
         end do

      case ('atb')

         ! Compute B = X^T*B.

         do k = 1, nblocks

            if (X%bs(k+1) > 0) then
               i1 = abs(X%bs(k))
               i2 = abs(X%bs(k+1))-1
               
               ! B = (I - uk*uk^T)*B.

               call Hprod(i1, i2, X%u, nc, B, ldb)

               ! B = Sk*B.

               do i = i1, i2
                  B(i,1:nc) = X%sig(i)*B(i,1:nc)
               end do

               ! B = (I - vk*vk^T)*B.

               call Hprod(i1, i2, X%v, nc, B, ldb)
              
            end if
         end do

      case ('aib')

         ! Compute B = X^-1*B

         do k = 1, nblocks

            if (X%bs(k+1) > 0) then
               i1 = abs(X%bs(k))
               i2 = abs(X%bs(k+1))-1
               
               ! B = (I - uk*uk^T)*B.

               call Hprod(i1, i2, X%u, nc, B, ldb)

               ! B = Sk^-1*B.

               do i = i1, i2
                  B(i,1:nc) = B(i,1:nc)/X%sig(i)
               end do

               ! B = (I - vk*vk^T)*B.

               call Hprod(i1, i2, X%v, nc, B, ldb)
            end if
         end do

      case ('aitb')
         
         ! Compute B = X^-T*B.

         do k = 1, nblocks

            if (X%bs(k+1) > 0) then
               i1 = abs(X%bs(k))
               i2 = abs(X%bs(k+1))-1
               
               ! B = (I - vk*vk^T)*B.

               call Hprod(i1, i2, X%v, nc, B, ldb)

               ! B = Sk^-1*B.

               do i = i1, i2
                  B(i,1:nc) = B(i,1:nc)/X%sig(i)
               end do
               
               ! B = (I - vk*vk^T)*B.

               call Hprod(i1, i2, X%u, nc, B, ldb)
               
            end if
         end do
         
      case default

         stop "Error in HsvdProd: Illegal operation."

      end select

   end subroutine HsvdProd

   subroutine Hprod(i1, i2, w, ncols, B, ldb)
      real(wp), intent(in) :: w(:)
      integer, intent(in)  ::  i1, i2, ncols, ldb 
      real(wp), intent(inout) :: B(ldb,*)

               
!     hprod computes the product of a block Householder transformation
!     in a hsvdmat with a matrix contained in B.  Specifically,
!
!        B[i1:i2,1:ncols] = (I - w[i1:i2]*w[i1:i2]^T)*B[i1:i2,1:ncols].
!
!     i1    Points to the beginning of the block.
!     i2    Points to the end of the block.
!     ncol  The number of columns of B
!     w     The vector from the hsvd.
!     B     The array containing the matrix B.

      real(wp) :: dot
      integer  :: i,j

      do j = 1,ncols
         dot = 0.0_wp
         do i = i1,i2
            dot = dot + w(i)*B(i,j)
         end do
         do i = i1,i2
            B(i,j) = B(i,j) - dot*w(i)
         end do
      end do

   end subroutine Hprod

   subroutine EigenmatVecs(A, eignum, eig, x, y, cond, job)
      type(Eigenmat), intent(in) :: A
      integer, intent(in)        :: eignum
      complex(wp), intent(out)   :: eig
      complex(wp), intent(inout) :: x(:)
      complex(wp), intent(inout) :: y(:)
      real(wp), intent(out)      :: cond
      character, intent(in)      :: job

!     EigenmatVecs computes an eigenvalue of A and the corresponding
!     left or right eigenvectors or principal vectors
!     as specified job.  If both left and rights
!     are computed,  EigenmatVecs also returns the condition
!     number of the eigenvalue  (or -1 if the eigenvalue belongs
!     to a Jordan block).  The eigenvectors are scaled
!     to have Euclidean norm 1.
!
!     A        The eigenmat whose vectors are to be computed.
!     eignum   The position in A%eig of the eigenvalue.
!     eig      The eigenvalue.
!     x(:)     The right eigenvector or principal vector.
!     y(:)     The left eigenvector or principal vector.
!     cond     The condition number of the eigenvalue.
!              (or -1, if the eigenvalue belongs to a
!              Jordan block).
!     job      A string specifying what to compute.
!
!              "r"  Compute the right eigenvector.
!              "l"  Compute the left eigenvector.
!              "b"  Compute both and the condition number.
!              (Note: For Jordan blocks, principal vectors
!              are computed and -1 is returned for the
!              condition number.)

      integer  :: n, i, found
      real(wp) :: xr(A%n), xi(A%n), yr(A%n), yi(A%n)
      real(wp) :: ZERO = 0.0_wp
      real(wp) :: ONE  = 1.0_wp

      n = A%n
      if (eignum > n) stop "Error in EigenmatVecs: eigen > n."

      if (job.ne.'b' .and. job.ne.'r' .and. job.ne.'l') &
           stop "Error in EigenmatVecs: Illegal operator."

      xr = ZERO
      xi = ZERO
      yr = ZERO
      yi = ZERO
      cond = ZERO

      if (A%type(eignum) .le. 1) then 

         ! Real eigenvalue.
         
         if (A%type(eignum) .eq. 1 .or. A%type(eignum) .lt. -1) then
            eig = A%eig(eignum)
         elseif (A%type(eignum) .eq. -1) then
            found = 0
            do i=eignum-1, 1, -1
               if (A%type(i) < -1) then
                  eig = A%eig(i)
                  found = 1
                  exit
               elseif (A%type(i) .ne. -1) then 
                  stop 'Error in EigenmatVecs: Illegal type.'
               end if
            end do
            if (found .eq. 0) then
               stop 'Error in EigenmatVecs: Illegal type.'
            end if
         else
            stop 'Error in EigenmatVecs: Illegal type.'
         end if
         
         if (job .eq. 'b' .or. job .eq. 'r') then 
            
            ! Compute the right eigenvector.
            
            xr(eignum) = ONE
            call HsvdProd(A%Z, 1, xr, n, 'ab')
            call HsvdProd(A%Y, 1, xr, n, 'ab')
            x(1:n) = cmplx(xr(1:n), xi(1:n), wp)
            if (A%type(eignum) .eq. 1) &
                 x(1:n) = x(1:n)/sqrt(dot_product(x(1:n),x(1:n)))
         end if
         
         if ( job .eq. 'b' .or. job .eq. 'l') then 
            
            ! Compute the left eigenvector.
            
            yr(eignum) = ONE
            call HsvdProd(A%Z, 1, yr, n, 'aitb')
            call HsvdProd(A%Y, 1, yr, n, 'aitb')
            y(1:n) = cmplx(yr(1:n), yi(1:n), wp)
            if (A%type(eignum) .eq. 1) &
                 y(1:n) = y(1:n)/sqrt(dot_product(y(1:n),y(1:n)))
         end if

         ! Compute the condition number.

         if (job .eq. 'b') then
            if (A%type(eignum) .eq. 1) then
               cond = ONE/dot_product(xr(:), yr(:))
            else
               cond = -1
            end if
         end if

      elseif (eignum .eq. n) then
         
         stop "Error in Eigematprod: Illegal type."

      elseif (A%type(eignum) .eq. 2) then
         if (A%type(eignum+1) .ne. 3) then
            stop "Error in Eigematprod: Illegal type."
         end if

         ! Complex eigenvalue.

         eig = cmplx(A%eig(eignum),A%eig(eignum+1), wp)

         if ( job .eq. 'b' .or. job .eq. 'r') then 

            ! Compute the right eigenvector.

            xr(eignum)   = 1.0
            xi(eignum+1) = 1.0
            call HsvdProd(A%Z, 1, xr, n, 'ab')
            call HsvdProd(A%Y, 1, xr, n, 'ab')
            call HsvdProd(A%Z, 1, xi, n, 'ab')
            call HsvdProd(A%Y, 1, xi, n, 'ab')
            x(1:n) = cmplx(xr(:), xi(:), wp)
            x(1:n) = x(1:n)/sqrt(dot_product(x(1:n),x(1:n)))
         end if

         if ( job .eq. 'b' .or. job .eq. 'l') then 

            ! Compute the left eigenvector.

            yr(eignum)   = 1.0
            yi(eignum+1) = 1.0
            call HsvdProd(A%Z, 1, yr, n, 'aitb')
            call HsvdProd(A%Y, 1, yr, n, 'aitb')
            call HsvdProd(A%Z, 1, yi, n, 'aitb')
            call HsvdProd(A%Y, 1, yi, n, 'aitb')
            y(1:n) = cmplx(yr(:), yi(:), wp)
            y(1:n) = y(1:n)/sqrt(dot_product(y(1:n),y(1:n)))
         end if

         ! Compute the condition number.

         if (job .eq. 'b') then 
            cond = ONE/abs(dot_product(y(1:n), x(1:n)))
         end if

      else
         stop "Error in Eigematprod: Illegal type."
      end if

   end subroutine EigenmatVecs

   
   subroutine hscal(u)

      !
      ! hscal scales the input vector u to norm sqrt(2).
      !

      real(wp), intent(inout) :: u(:)
      
      real(wp) :: nrm, mx

      ! Find the maximum magnitue of u

      mx = maxval(abs(u))
      
      u = u/mx

      ! Compute u'*u/mx^2

      nrm = dot_product(u, u)

      if (nrm .eq. 0.0_wp) then
         stop 'Error in hnorm: zero argument.'
      end if
      
      ! Scale u

      nrm = sqrt(2.0_wp)/sqrt(nrm)

      u = nrm*u

   end subroutine hscal
   
end module eigentest

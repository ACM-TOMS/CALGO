* Eigentest is a package that produces real test matrices with known
* eigensystems.  A test matrix, called an eigenmat, is generated in a
* factored form, in which the user can specify the eigenvalues
* (including complex conjugate eigenvalues and Jordan blocks)and has
* some control over the condition of the eigenvalues and eigenvectors.
* An eigenmat $A$ of order $n$ requires only $O(n)$ storage for its
* representation.  Auxiliary programs permit the computation of
*
*  (A - s*I)*b, (A - s*I)'*b, inv(A - s*I)*b, and inv(A -s*I)'*b
*
* in $O(n)$ operations.  Thus eigenmats are suitable for testing
* algorithms based on Krylov sequences, as well as others based on
* matrix vector products.  A special routine computes specified
* eigenvectors of an eigenmat and the condition of its eigenvalue.
*
* Contents
*
*    subroutine eminit
*       Initializes an eigenmat.
*
*    subroutine emprod
*       Computes the above products.
*
*    subroutine hsvdpr
*       Computes products of special matrices in the factorization.
*
*    subroutine hsprod
*       Computes products with a Householder transformation.
*
*    subroutine emvecs
*       Computes eigenvectors of an eigenmat.
*
* For more details see
*
*    1. Eigentest: A Test Maatrix Generator for Large Scale
*       Eigenproblem
*
*    2. Eigentest: The Fortran 95 User's Guide
*
* Pdf files are included with this distribution.

* Coded by Roger Lee and Pete Stewart


      subroutine eminit(nmax, n, iem, fem,  nblk, idy, idz, ia, fa, job)
      implicit none

      integer nmax, n, iem(*), nblk, idy, idz, ia(*)
      double precision fem(*), fa(*)
      character job*(*)


*     eminit is a utility routine to initialize the Fortran 77
*     eigentest arrays iem and fem.  These arrays encode an eigenmat, a
*     matrix of the form
*
*        A = (Y*Z)*(L - shift*I)/(Y*Z).
*
*     Here L is a block triangular with 1x1 blocks containing the real
*     eigenvalues of A and 2x2 blocks conaining the complex
*     eignevalues of A in the form
*
*                 |  mu  nu |
*                 |         |.
*                 | -nu  mu |
*
*     Y and Z are hsvd matrices (see below).
*
*     The eigemat is stored in arrays iem and fem as follows.
*
*        iem(1)     nmax     The maximum size of the matrix.
*        iem(2)     n        The order of the matrix.
*        iem(3)     type     The type of the entry in
*                            eig (see below).
*                            eig contain mu and nu as above.
*        iem(nmax+3)         The beginning of the integer
*                            array for the hsvdmat Y.
*        iem(nmax+8)         The beginning of the integer array
*                            for the hsvdmat Z.
*
*        fem(1)        eig   The array containing the eigenvalues
*                            of the matrix.
*        fem(nmax+1)         The beginning of the double array
*                            for the hsvdmat Y.
*        fem(4*nmax+1)       The beginning of the double array
*                            for the hsvdmat Z.
*        fem(7*nmax+1)       The beginning of a work array of
*                            length nmax.
*  The contents of eig and type are determined as follows.
*
*  Real eigenvalue
*     type(i) = 1    eig(i) = the eigenvalue.
*
*  Complex conjugate pair
*     type(i) = 2    eig(i) = the real part.
*     type(i) = 3    eig(i) = the complex part.
*
*  Jordan block of order k
*     type(i) = -k     eig(i)   = the eigenvalue.
*     type(i+j) = -1   eig(i+j) = the jth superdiagonal
*                i = 1,...,k-1.
*
*
*     A hsvd matrix has the form
*
*         X = diag(X1, X2, ..., Xnblk),
*
*     where
*
*         Xk = (I - uk*uk^T)*Sk*(I - vk*vk^T).
*
*     with Sk diagonal having positive diagonal entries.  The
*     hsvd is stored in the arrays iem, and id as follows.  This
*     storage is realative to the the origin of the hvdmat in
*     the description of the eigenmat above.  For example, iem(5)
*     in the table below for Y corresponds to iem(nmax+3+5-1) in
*     the eigenmat table.  Again, fem(1+nmax) for Z corresponds
*     to fem(4*nmax+1+1+nmax-1) in the eigenmat table.
*
*        iem(1)          nmax   maximum value of n.
*        iem(2)          n      The order of the hsvd.
*        iem(3)          nblk   The number of blocks in the hsvd.
*        iem(4)          bs     Beginning of the block start array.
*                               abs(bs(5))=1 and abs(bs(i+1))-abs(bs(i))
*                               is the size of the ith block.
*                               If bs(i+1) is negative, block i is
*                               an identity.
*        fem(1)          u      start of u.
*        fem(1+nmax)     v      start of v.
*        fem(1+2*nmax)   sig    start of the singular values.
*
*     The blocks uk of u are stacked in fem their natural order starting
*     at fem(1).  Similary for v and sig.
*
*     There will be sufficient storage to contain the structures
*     if the length of iem is of length 2*nmax+12 and fem is of
*     length 7*nmax.
*      
*     USAGE: Storage for an eigenmat may be defined as in the following
*     example.
* 
*        parameter nmax = <The order of the largest matrix>
*        integer iem(2*nmax+12), ia(nmax+1)
*        double precision fem(8*nmax), fa(nmax+1)
*
*     eminit can be used to fill the arrays iem and fem as follows.
*     You will need arrays
*
*        integer ia(nmax+1)
*        double precision fa(nmax)
*
*     Assign values to nmax, n, nblk, idy, idz,
*
*     where nblk is the number of blocks in the hsvdmat Z.  Then
*     starting with ia(1) fill in the bs array for Z and
*
*        call eminit(nmax, n, iem, fem,  nblk, idy, idz, ia, fa, 'setup')
*         
*     Next place the u vector for Y in fa and
*
*        call eminit(nmax, n, iem, fem,  nblk, idy, idz, ia, fa, 'yu')
*
*     Continue in this manner with the other parts of the structure.
*     Finally place type and eig in ia and fa and
*
*        call eminit(nmax, n, iem, fem,  nblk, idy, idz, ia, fa, 'eig')
*
*     and the eigenmat is initialized.
*
*     Actually, the calls, after the first, can be made
*     in any order, provided nmax and n retain
*     their original values.
*

      integer i

      if (job .eq. 'setup') then

         iem(1) = nmax
         iem(2) = n

*        Initialize Y.

         iem(nmax+3) = nmax
         iem(nmax+4) = n
         iem(nmax+5) = 1
         iem(nmax+6) = 1
         if (idy .ne. 0) then
            iem(nmax+7) = -(n+1)
         else
            iem(nmax+7) = n+1
         end if

*        initialize Z.

         iem(nmax+8) = nmax
         iem(nmax+9) = n

         if (idz .ne. 0) then
            iem(nmax+10) = 1
            iem(nmax+11) = 1
            iem(nmax+12) = -(n+1)
         else
            iem(nmax+10) = nblk
            do i=1,nblk+1
               iem(i+nmax+10) = ia(i)
            end do
         end if

*     Initialize Y arrays.

      else if(job .eq. 'yu') then
         do i=1,n
            fem(nmax+i) = fa(i)
         end do

      else if(job .eq. 'yv') then
         do i=1,n
            fem(2*nmax+i) = fa(i)
         end do

      else if(job .eq. 'ysig') then
         do i=1,n
            fem(3*nmax+i) = fa(i)
         end do

*     initialize Z arrays.

      else if(job .eq. 'zu') then
         do i=1,n
            fem(4*nmax+i) = fa(i)
         end do

      else if(job .eq. 'zv') then
         do i=1,n
            fem(5*nmax+i) = fa(i)
         end do

      else if(job .eq. 'zsig') then
         do i=1,n
            fem(6*nmax+i) = fa(i)
         end do

*     Initalize the eigenvalues.

      else if (job .eq. 'eig') then

         do i=1,n
            fem(i) = fa(i)
            iem(i+2) = ia(i)
         end do

      else
         stop 'Error in eminit: Illegal value of job.'
      end if

      return
      end


**********************************************************************
*            
**********************************************************************

      subroutine emprod(iem, fem, ncols, B, ldb, C, ldc, shift, op)
      implicit none

      integer iem(*), ncols, ldb, ldc
      double precision fem(*), B(ldb,*), C(ldc,*), shift
      character op*(*)

*     emprod computes the product, perhaps shifted and inverted, of an
*     eigenmat A contained in iem and fem with a matrix B, placing the
*     result in C.  An eigenmat has the form
*
*        A = (Y*Z)*(L - shift*I)/(Y*Z)
*
*     Here L is a block triangular with 1x1 blocks containing the real
*     eigenvalues of A and 2x2 blocks conaining the complex
*     eignevalues of A in the form
*
*                 |  mu  nu |
*                 |         |.
*                 | -nu  mu |
*
*     Y and Z are hsvd transformations.
*
*     The eigemat is stored in arrays iem and fem as follows.
*
*        iem(1)     nmax     The maximum size of the matrix.
*        iem(2)     n        The order of the matrix.
*        iem(3)     type     The beginning of an array
*                            containing the types of the
*                            eigenvalues in eig.  If type(i)
*                            is 1, the corresponding entry of
*                            eig is a real eigenvalue.  If
*                            type(i) is 2 and type(i+1) is 3
*                            then the corresponding entries of
*                            eig contain mu and nu as above.
*        iem(nmax+3)         The begining of the integer
*                            array for the hsvdemat Y.
*        iem(nmax+8)         The beginning of the integer array
*                            for the hsvd mat Z.
*        fem(1)        eig   The array containing the eigenvalues
*                            of the matrix.
*        fem(nmax+1)         The beginning of the double array
*                            for the hsvdmat Y.
*        fem(4*nmax+1)       The beginning of the double array
*                            for the hsvdmat Z.
*
*     For more on the storge of hsvdmats, see eminit and hsvdpr.
*
*     Arguments.
*
*     iem,fem   Arrays containing the eigenmat.
*     ncols     Number of columns in the matrix B.
*     B         Pointer to an array containing the matrix B.
*     ldb       The leading dimension of the array B.
*     C         Pointer to an array containing the matrix C.
*     ldc       The leading dimension of C.
*     shift     A shift.
*     op        A string specifying the operation to be performed.
*
*               'ab'    C = (A - shift*I)*B
*               'atb'   C = (A - shift*I)'*B
*               'aib'   C = (A - shift*I)\B
*               'aitb'  C = (A - shift*I)'\B
*     
*     Local variables

      integer eigp, i, j, k, endj, typep, yfp, yip, zfp, zip, nmax, n
      double precision lam, mu, nu, temp, mult, u11, u12, u22
*
*     Initialize constants from iem.
*
      nmax = iem(1)
      n = iem(2)
*
*     Initialize pointers into iem and fem.
*
      typep = 3
      yip = nmax + 3
      zip = nmax + 8
      eigp = 1
      yfp = nmax + 1
      zfp = 4*nmax + 1
*
*     Copy B into C.
*
      do j=1,ncols
         do i=1,n
            C(i,j) = B(i,j)
         end do
      end do

      if (op .eq. 'ab') then
*
*        Compute (A - shift*I)*C = Y*(Z*((L-shift*I)*(Z\(Y\C)))).
*
*        Compute Z\(Y\C).
*
         call hsvdpr(iem(yip), fem(yfp), ncols, C, ldc, 'aib')
         call hsvdpr(iem(zip), fem(zfp), ncols, C, ldc, 'aib')

         i = 1 
         do while (i .le. n)

            if (iem(typep+i-1) .eq. 1) then
*
*              Real eigenvalue.
*
               lam = fem(eigp+i-1) - shift
               do j=1,ncols
                  C(i,j) = lam*C(i,j)
               end do
               i = i+1

            else if (i .eq. n) then
               stop 'Error in emprod: 2x2 block starts at bs(n).'

            else if (iem(typep+i-1).eq.2 .and. iem(typep+i).eq.3) then
*
*              Complex eigenvalue.
*
               mu = fem(eigp+i-1) - shift
               nu = fem(eigp+i)
               do j=1,ncols
                  temp = mu*C(i,j) + nu*C(i+1,j)
                  C(i+1,j) = -nu*C(i,j) + mu*C(i+1,j)
                  C(i,j) = temp
               end do
               i = i+2

            else if (iem(typep+i-1).lt.-1) then
*
*              Jordan block.
*
               endj = i + abs(iem(typep+i-1))-1
               if (endj .gt. n) then
                  stop 'Error in emprod: Jordan block too large.'
               end if
               
               do j = i, endj-1
                  if (iem(typep+j) .ne. -1) then
                     stop 'Error in emprod: Illegal type.'
                  end if
               end do

               mu = fem(eigp+i-1) - shift

               do j=1,ncols
                  do k=i,endj-1
                     C(k,j) = mu*C(k,j)+ fem(eigp+k)*C(k+1,j)
                  end do
                  C(endj,j) = mu*C(endj,j)
               end do
               i = endj + 1

            else
               stop 'Error in emprod: Illegal type.'
            end if

         end do
*
*        Compute C = Y*(Z*C).
*
         call hsvdpr(iem(zip), fem(zfp), ncols, C, ldc, 'ab')
         call hsvdpr(iem(yip), fem(yfp), ncols, C, ldc, 'ab')

      else if (op .eq. 'atb') then
* 
*        C = (A-shift*I)'*C = Y'\(Z'\(((L-shift*I)'*(Z'*(Y'*C)))).
*
*        Compute C = Z'*(Y'*C).

         call hsvdpr(iem(yip), fem(yfp), ncols, C, ldc, 'atb')
         call hsvdpr(iem(zip), fem(zfp), ncols, C, ldc, 'atb')
*
*        Compute (L-shift*I)*C.  The index i points to successive blocks
*        of L.
*
         i=1
         do while (i .le. n)
*
*           Real eigenvalue.
*  
            if (iem(typep+i-1) .eq. 1) then
               lam = fem(eigp+i-1) - shift
               do j=1,ncols
                  C(i,j) = lam*C(i,j)
               end do
               i = i+1

            else if (i .eq. n) then
               stop 'Error in emprod: 2x2 block starts at bs(n).'

            else if (iem(typep+i-1).eq.2 .and. iem(typep+i).eq.3) then
*
*              Complex eigenvalue.
*
               mu = fem(eigp+i-1) - shift
               nu = fem(eigp+i)
               do j=1,ncols
                  temp = mu*C(i,j) - nu*C(i+1,j)
                  C(i+1,j) = nu*C(i,j) + mu*C(i+1,j)
                  C(i,j) = temp
               end do
               i = i+2

            else if (iem(typep+i-1).lt.-1) then
*
*              Jordan block.
*
               endj = i + abs(iem(typep+i-1))-1
               if (endj .gt. n) then
                  stop 'Error in emprod: Jordan block too large.'
               end if

               do j = i, endj-1
                  if (iem(typep+j) .ne. -1) then
                     stop 'Error in emprod: Illegal type.'
                  end if
               end do

               mu = fem(eigp+i-1) - shift

               do j=1,ncols
                  do k=endj,i+1,-1
                     C(k,j) = mu*C(k,j)+fem(eigp+k-1)*C(k-1,j)
                  end do
                  C(i,j) = mu*C(i,j)
               end do
               i = endj + 1

            else
               stop 'Error in emprod: Illegal type.'
            end if
         end do
*
*        Compute C = Y'\(Z'\C).
*
         call hsvdpr(iem(zip), fem(zfp), ncols, C, ldc, 'aitb')
         call hsvdpr(iem(yip), fem(yfp), ncols, C, ldc, 'aitb')

      else if (op .eq. 'aib') then
*
*        C = (A-shift*I)*C = Y*(Z*((L-shift*I)\(Z\(Y\C)))).
*
*        Compute C = Z\*(Y\C).
*
         call hsvdpr(iem(yip), fem(yfp), ncols, C, ldc, 'aib')
         call hsvdpr(iem(zip), fem(zfp), ncols, C, ldc, 'aib')
*
*        Compute (L-shift*I)\C.  The index i points to successive blocks
*        of L.
*
         i=1
         do while (i .le. n)
            if (iem(typep+i-1) .eq. 1) then
*
*              Real eigenvalue.
*
               lam = fem(eigp+i-1) - shift
               do j=1,ncols
                  C(i,j) = C(i,j)/lam
               end do
               i = i+1

            else if (i .eq. n) then
               stop 'Error in emprod: 2x2 block starts at bs(n).'

            else if (iem(typep+i-1).eq.2 .and. iem(typep+i).eq.3) then
*
*              Complex eigenvalue.  Use Gaussian elimination with
*              partial pivoting and back substitution to apply inv(L).
*
               mu = fem(eigp+i-1) - shift
               nu = fem(eigp+i)
               if (abs(mu) .ge. abs(nu)) then
*
*                 No pivoting.
*
                  u11 = mu
                  u12 = nu
                  mult = -nu/mu
                  u22 = mu - mult*nu
                  do j=1,ncols
                     C(i+1,j) = (C(i+1,j) - mult*C(i,j))/u22
                     C(i,j) = (C(i,j) - u12*C(i+1,j))/u11
                  end do
               else
*
*                 Pivoting.
*
                  u11 = -nu
                  u12 = mu
                  mult = -mu/nu
                  u22 = nu - mult*mu
                  do j=1,ncols
                     temp = C(i+1,j)
                     C(i+1,j) = (C(i,j) - mult*temp)/u22
                     C(i,j) = (temp - u12*C(i+1,j))/u11
                  end do
               end if
               i = i+2

            else if (iem(typep+i-1).lt.-1) then
*
*              Jordan block.
*
               endj = i + abs(iem(typep+i-1))-1
               if (endj .gt. n) then
                  stop 'Error in emprod: Jordan block too large.'
               end if

               do j = i, endj-1
                  if (iem(typep+j) .ne. -1) then
                     stop 'Error in emprod: Illegal type.'
                  end if
               end do

               mu = fem(eigp+i-1) - shift

               do j=1,ncols
                  C(endj,j) = C(endj,j)/mu
                  do k=endj-1, i, -1
                     C(k,j) = (C(k,j)-C(k+1,j)*fem(eigp+k))/mu
                  end do
               end do

               i = endj + 1

            else
               stop 'Error in emprod: Illegal type.'
            end if
         end do
*
*        Compute C = Y*(Z*C).
*
         call hsvdpr(iem(zip), fem(zfp), ncols, C, ldc, 'ab')
         call hsvdpr(iem(yip), fem(yfp), ncols, C, ldc, 'ab')

      else if (op .eq. 'aitb') then
*
*        compute C = A^-T*C = Y'\(Z'\((L-shift*I)'\(Z'*(Y'*C)))).
*
*        Compute C = Z'*(Y'*C).
*
         call hsvdpr(iem(yip), fem(yfp), ncols, C, ldc, 'atb')
         call hsvdpr(iem(zip), fem(zfp), ncols, C, ldc, 'atb')
*
*        Compute (L-shift*I)'\C.  The index i points to successive blocks
*        of L.

         i=1
         do while (i .le. n)
            if (iem(typep+i-1) .eq. 1) then
*
*              Real eigenvalue.
*
               lam = fem(eigp+i-1) - shift
               do j=1,ncols
                  C(i,j) = C(i,j)/lam
               end do
               i = i+1

            else if (i .eq. n) then
               stop 'Error in emprod: 2x2 block starts at bs(n).'

            else if (iem(typep+i-1).eq.2 .and. iem(typep+i).eq.3) then
*
*              Complex eigenvalue.  Use Gaussian elimination with
*              partial pivoting and back substitution to apply inv(L).
*
               mu = fem(eigp+i-1) - shift
               nu = fem(eigp+i)
               if (abs(mu) .ge. abs(nu)) then
*
*                 No pivoting.
*
                  u11 = mu
                  u12 = -nu
                  mult = nu/mu
                  u22 = mu + mult*nu
                  do j=1,ncols
                     C(i+1,j) = (C(i+1,j) - mult*C(i,j))/u22
                     C(i,j) = (C(i,j) - u12*C(i+1,j))/u11
                  end do
               else
*
*                 Pivoting.
*
                  u11 = nu
                  u12 = mu
                  mult = mu/nu
                  u22 = -nu - mult*mu
                  do j=1,ncols
                     temp = C(i+1,j)
                     C(i+1,j) = (C(i,j) - mult*temp)/u22
                     C(i,j) = (temp - u12*C(i+1,j))/u11
                  end do
               end if
               i = i+2

            else if (iem(typep+i-1) .lt. -1) then
*
*              Jordan block.
*
               endj = i + abs(iem(typep+i-1))-1
               if (endj .gt. n) then
                  stop 'Error in emprod: Jordan block too large.'
               end if

               do j = i, endj-1
                  if (iem(typep+j) .ne. -1) then
                     stop 'Error in emprod: Illegal type.'
                  end if
               end do

               mu = fem(eigp+i-1) - shift

               do j=1,ncols
                  C(i,j) = C(i,j)/mu
                  do k=i+1, endj
                     C(k,j) = (C(k,j)-C(k-1,j)*fem(eigp+k-1))/mu
                  end do
               end do

               i = endj + 1

            else
               stop 'Error in emprod: Illegal type.'
            end if
         end do
*
*        Compute C = Y'\(Z'\C).
*
         call hsvdpr(iem(zip), fem(zfp), ncols, C, ldc, 'aitb')
         call hsvdpr(iem(yip), fem(yfp), ncols, C, ldc, 'aitb')

      else
         stop 'Error in emprod: Illegal operation.'
      end if

      return
      end

**********************************************************************
*            
**********************************************************************

      subroutine hsprod(iem, fem, ncols, B, ldb, il, start, nb)
      implicit none
      
      integer iem(*), ncols, ldb, il, start, nb
      double precision fem(*), B(ldb,*)      

      double precision dot
      integer i, j
      
      do j = 1, ncols
         
         dot = 0.0D0
         do i = 1, nb
            dot = dot + fem(start+i-1)*B(il+i-1,j)
         end do

         do i = 1, nb
            B(il+i-1,j) = B(il+i-1,j) - dot*fem(start+i-1)
         end do

      end do
      
      end
      
**********************************************************************
*            
**********************************************************************

      subroutine hsvdpr(iem, fem, ncols, B, ldb, op)
      implicit none

      integer iem(*), ncols, ldb
      double precision fem(*), B(ldb,*)
      character op*(*)
      
*     hsvdprod computes the product of computes the product of an
*     hsvdmat X contained in the arrays iem and fem and a matrix B,
*     overwriting B with the product.  A hsvd matrix has the form
*
*         X = diag(X1, X2, ..., Xnblk),
*
*     where
*
*         Xk = (I - uk*uk')*Sk*(I - vk*vk').
*
*     with Sk diagonal having positive diagonal entries.  The
*     hsvd is stored in the arrays iem, and fsv as follows.
*
*        iem(1)          nmax   maximum value of n
*        iem(2)          n      The order of the hsvd
*        iem(3)          nblk   The number of blocks in the hsvd
*        iem(4)          bs     Beginning of the block-start array.
*                               abs(bs(5))=1 and abs(bs(i+1))-abs(bs(i))
*                               is the size of the ith block.
*                               If bs(i+1) is negative, block i is
*                               an identity.
*        fem(1)          u      start of u
*        fem(1+nmax)     v      start of v
*        fem(1+2*nmax)   sig    start of the singular values
*
*     The blocks uk of u are stacked in fem their natural order starting
*     at fem(1).  Similary for v and sig.
*
*     Parameters:
*
*     iem, fem (in)   The hsvd
*
*     ncols (in)    The number of columns in B.
*
*     B     (inout) The array B.
*
*     ldb   (in)    The leading dimension of B.
*
*     op    (in)    A string specifying the operation to be performed.
*
*                   'ab'   B <- X*B
*                   'atb'  B <- X^t*B
*                   'aib'  B <- X^-1*B
*                   'aitb' B <- X^-1T*B
*
*     Local variables:
*

      integer bsp, sigp, up, vp, i, il, iu, j, k, n, nblk, nmax, nb
*
*     Initialize constants from iem.
*
      nmax = iem(1)
      n = iem(2)
      nblk = iem(3)
*
*     Initialize pointers into iem and fem.
*
      bsp = 4
      up = 1
      vp = up + nmax
      sigp = vp + nmax

      il = abs(iem(bsp))

      if (op .eq. 'ab') then
*
*        Compute B = X*B.
*
         do k = bsp, bsp+nblk-1
            iu = abs(iem(k+1))
            nb = iu - il

            if (iem(k+1) .gt. 0) then
               call hsprod(iem, fem, ncols, B, ldb, il, vp, nb)
               
               do j=1,ncols
                  do i=il,iu-1
                     B(i,j) = fem(sigp+i-il)*B(i,j)
                  end do
               end do
               
               call hsprod(iem, fem, ncols, B, ldb, il, up, nb)

            end if
            vp = vp + nb
            up = up + nb
            sigp = sigp + nb
            il = iu
         end do

      else if (op .eq. 'atb') then
*
*        Compute B = X'*B.
*
         do k = bsp, bsp+nblk-1
            iu = abs(iem(k+1))
            nb = iu - il
            if (iem(k+1) .gt. 0) then

               call hsprod(iem, fem, ncols, B, ldb, il, up, nb)

                do j=1,ncols
                  do i=il,iu-1
                     B(i,j) = fem(sigp+i-il)*B(i,j)
                  end do
               end do

               call hsprod(iem, fem, ncols, B, ldb, il, vp, nb)

            end if
            vp = vp + nb
            up = up + nb
            sigp = sigp + nb
            il = iu
         end do

      else if (op .eq. 'aib') then
*
*        Compute X\B.
*
         do k = bsp, bsp+nblk-1
            iu = abs(iem(k+1))
            nb = iu - il
            if (iem(k+1) .gt. 0) then

               call hsprod(iem, fem, ncols, B, ldb, il, up, nb)
              
                do j=1,ncols
                  do i=il,iu-1
                     B(i,j) = B(i,j)/fem(sigp+i-il)
                  end do
               end do

               call hsprod(iem, fem, ncols, B, ldb, il, vp, nb)
              
            end if
            vp = vp + nb
            up = up + nb
            sigp = sigp + nb
            il = iu
         end do

      else if (op .eq. 'aitb') then
*
*        Compute X'\B.
*
         do k = bsp, bsp+nblk-1
            iu = abs(iem(k+1))
            nb = iu - il

            if (iem(k+1) .gt. 0) then

               call hsprod(iem, fem, ncols, B, ldb, il, vp, nb)

               do j=1,ncols
                  do i=il,iu-1
                     B(i,j) = B(i,j)/fem(sigp+i-il)
                  end do
               end do
               
               call hsprod(iem, fem, ncols, B, ldb, il, up, nb)
               
            end if
            vp = vp + nb
            up = up + nb
            sigp = sigp + nb
            il = iu
         end do
      else
         stop 'Error in hsvdpr: Illegal operation.'
      end if

      return
      end


**********************************************************************
*            
**********************************************************************

      subroutine emvecs(iem, fem, eignum, eig, x, y, cond, job)
      implicit none
      integer          iem(*), eignum
      double precision fem(*), cond
      double complex   eig, x(*), y(*)
      character        job

*     emvecs computes an eigenvalue the and the corresponding
*     left or right eigenvectors as specified job.  If both
*     are computed,  emvecs also returns the condition
*     number of the eigenvalue.  The eigenvectors are scaled
*     to have Euclidean norm 1.
*
*     iem, fem Arrays containint the eigenmat whose vectors
*              are to be computed.
*     eignum   The position in A%eig of the eigenvalue.
*     eig      The eigenvalue.
*     x(:)     The right eigenvector or principal vector.
*     y(:)     The left eigenvector or principal vector.
*     cond     The condition number of the eigenvalue.
*              (or -1, if the eigenvalue belongs to a
*              Jordan block).
*     job      A string specifying what to compute.
*
*              'r'  Compute the right eigenvector.
*              'l'  Compute the left eigenvector.
*              'b'  Compute both and the condition number.
*                   (Note: For Jordan blocks, principal vectors
*                   are computed and -1 is returned for the
*                   condition number.)

      integer i, found
      integer nmax, n, yip, yfp, zip, zfp, wkp
      double precision nrm


      nmax = iem(1)
      n    = iem(2)
      
*     Set up pointers to the inner and outter hsvds and
*     to the work array.

      yip  = nmax + 3
      zip  = nmax + 8
      yfp  = nmax + 1
      zfp  = 4*nmax + 1

      wkp  = 7*nmax + 1

*     Compute the eigenvectors.

      if (eignum .gt. n) then
         stop 'Error in emvecs: eignum > n.'
      end if
      
      if (job.ne.'b' .and. job.ne.'r' .and. job.ne.'l') then
         stop 'Error in emvecs: Illegal operator.'
      end if
*
*     Initialize variables.
*

      if (iem(2+eignum) .le. 1) then
*
*     Real eigevalue or in a Jordan block
*
         if (iem(2+eignum) .eq. 1 .or. iem(2+eignum) .lt. -1) then
            eig = dcmplx(fem(eignum))
         
         elseif (iem(2+eignum) .eq. -1) then
            found = 0
            do i = eignum-1, 1, -1
               if (iem(2+i) < -1) then
                  eig =  dcmplx(fem(i))
                  found = 1
                  exit
               elseif (iem(2+i) .ne. -1) then 
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
*         
*           Compute the right eigenvector.
*
            do i = 1, n
               fem(wkp+i-1) = 0.0D0
            end do
            fem(wkp+eignum-1) = 1.0D0
            call hsvdpr(iem(zip), fem(zfp), 1, fem(wkp), n, 'ab')
            call hsvdpr(iem(yip), fem(yfp), 1, fem(wkp), n, 'ab')            

            if (iem(2+eignum) .eq. 1) then
               nrm = 0.0D0
               do i = 1, n
                  nrm = nrm + fem(wkp+i-1)*fem(wkp+i-1)
               end do
               nrm = sqrt(nrm)
               do i = 1, n
                  x(i) = dcmplx(fem(wkp+i-1)/nrm)
               end do
            else
               do i = 1, n
                  x(i) = fem(wkp+i-1)
               end do
            end if

         end if
            
         if (job .eq. 'b' .or. job .eq. 'l') then
*
*           Compute the left eigenvector.         
*
            do i = 1, n
               fem(wkp+i-1) = 0.0D0
            end do

            fem(wkp+eignum-1) = 1.0D0
            call hsvdpr(iem(zip), fem(zfp), 1, fem(wkp), n, 'aitb')
            call hsvdpr(iem(yip), fem(yfp), 1, fem(wkp), n, 'aitb')

            if (iem(2+eignum) .eq. 1) then
               nrm = 0.0D0
               do i = 1, n
                  nrm = nrm + fem(wkp+i-1)*fem(wkp+i-1)
               end do
               nrm = sqrt(nrm)
               do i = 1, n
                  y(i) = dcmplx(fem(wkp+i-1)/nrm)
               end do
            else
               do i = 1, n
                  y(i) = fem(wkp+i-1)
               end do
            end if

         end if

         if (job .eq. 'b') then
            
            if (iem(2+eignum) .eq. 1) then
               cond = 0.0D0
               do i = 1, n
                  cond = cond + dble(y(i))*dble(x(i))
               end do
               cond = 1.0D0/abs(cond)
            else
               cond = -1.0D0
            endif
         end if

      else if (eignum .eq. n) then
         
         stop 'Error in emvecs: type error.'

      else if (iem(2+eignum) .eq. 2) then
         if (iem(2+eignum+1) .ne. 3) then
            stop 'Error in emvecs: type error.'
         end if
*
*        Complex eigenvalue.         
*
         eig = dcmplx(fem(eignum), fem(eignum+1))

         if (job .eq. 'b' .or. job .eq. 'r') then
*
*           Compute the right eigenvector.
*         

            do i = 1, n
               fem(wkp+i-1) = 0.0D0
            end do
            fem(wkp+eignum-1) = 1.0D0
            call hsvdpr(iem(zip), fem(zfp), 1, fem(wkp), n, 'ab')
            call hsvdpr(iem(yip), fem(yfp), 1, fem(wkp), n, 'ab')
            nrm = 0.0D0
            do i = 1, n
               nrm = nrm + fem(wkp+i-1)*fem(wkp+i-1)
               x(i) = dcmplx(fem(wkp+i-1))
               fem(wkp+i-1) = 0.0D0
            end do
            fem(wkp+eignum) = 1.0D0
            call hsvdpr(iem(zip), fem(zfp), 1, fem(wkp), n, 'ab')
            call hsvdpr(iem(yip), fem(yfp), 1, fem(wkp), n, 'ab')
            do i = 1, n
               nrm = nrm + fem(wkp+i-1)*fem(wkp+i-1)
            end do
            nrm = sqrt(nrm)
            do i = 1, n
               x(i) = dcmplx(dble(x(i)),fem(wkp+i-1))/nrm
            end do
         end if
            
         if (job .eq. 'b' .or. job .eq. 'l') then
*
*           Compute the left eigenvector.
*         
            do i = 1, n
               fem(wkp+i-1) = 0.0D0
            end do
            fem(wkp+eignum-1) = 1.0D0
            call hsvdpr(iem(zip), fem(zfp), 1, fem(wkp), n, 'aitb')
            call hsvdpr(iem(yip), fem(yfp), 1, fem(wkp), n, 'aitb')
            do i = 1, n
               nrm = nrm + fem(wkp+i-1)*fem(wkp+i-1)
               y(i) = dcmplx(fem(wkp+i-1))
               fem(wkp+i-1) = 0.0D0
            end do
            fem(wkp+eignum) = 1.0D0
            call hsvdpr(iem(zip), fem(zfp), 1, fem(wkp), n, 'aitb')
            call hsvdpr(iem(yip), fem(yfp), 1, fem(wkp), n, 'aitb')
            do i = 1, n
               nrm = nrm + fem(wkp+i-1)*fem(wkp+i-1)
            end do
            nrm = sqrt(nrm)
            do i = 1, n
               y(i) = dcmplx(dble(y(i)),fem(wkp+i-1))/nrm
            end do
         end if

         if (job .eq. 'b') then
*
*           Compute the condition number.
*
            cond = 0.0D0
            do i = 1, n
               cond = cond + dconjg(y(i))*x(i)
            end do
            cond = 1.0D0/cond
         end if

      else
         stop 'Error in emvecs: type error.'
      end if


      end

**********************************************************************
*            
**********************************************************************

      subroutine hscal(n, fa)
      implicit none
*
*     hscal scales the vector fa such that norm(fa)= sqrt(2). 
*     n is the length of fa.
*
      integer n
      double precision fa(*)

      integer i
      double precision nrm, mx, tmp

*
*     Find the maximum magnitude element of fa
*
      mx = abs(fa(1))
      do i = 2, n
         tmp = abs(fa(i))
         if (mx < tmp) then
            mx = tmp
         end if
      end do

*
*     Compute nrm = fa'*fa/max^2
*
      nrm = 0.0D0
      do i = 1, n
         tmp = fa(i)/mx
         nrm = nrm + tmp*tmp
      end do
         
      if (nrm .eq. 0.0D0) then
         stop 'Error in hscal: zero argument.'
      end if

*
*     Scale fa 
*
      nrm = sqrt(2.0D0/nrm)/mx
      do i = 1, n
         fa(i) = fa(i)*nrm
      end do
      end

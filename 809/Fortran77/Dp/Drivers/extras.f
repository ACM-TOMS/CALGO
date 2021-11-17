c
c***********************************************************************
c
c     This file contains auxiliary Fortran routines for the 
c     drivers. 
c
c     List of routines
c
c     rhsg    Simulated linear rhs. 
c     perrhs  perturbed rhs by 5% with random sign
c
c             WARNING.  This routine may require a modification. It 
c             makes use of the non intrinsic function  drand(),  to 
c             generate random numbers. Therefore, the user should 
c             replace drand() by the appropriate function supported
c             by the Fortran compiler.
c
c     dsymvs  matrix vector product in sparse format
c     dnorma  infinity norm of a sparse matrix
c
c     dsymvs and  dnorma come from MINPACK-2
c
c     sparse  computes the sparse structure of a matrix in dense format
c     dnrmif  computation of the infinity norm of a vector
c
c***********************************************************************
c
      subroutine rhsg ( n, g )
       
      integer  n
      double   precision  g(n)
c
c***********************************************************************
c
c     This routine generates right hand sides depending on the value
c     of n, the number of variables. For case n = 170, the rhs has
c     zero components everywhere except the positions 34*i, i=1,...,5.
c     The entries at these positions have the value -8000.0d0
c     For any other value of n, the routine produces a right hand side
c     whose components depend almost linearly on the index i.
c     The first and last components are zero. The remaining components 
c     are defined as follows
c    
c             g_i = (i-1/n-2)x100,  i=2,...,n-1
c
c     PARAMETERS
c
c     N    integer variable specifying the dimension of the rhs
c
c     G    real array of dimension N containing the rhs on output.
c
c
c***********************************************************************
c
      integer  i
      double   precision fn2, zero
      data     zero /0.0d0/
c
c     Construct a rhs for the specific case n = 170
c 
      if ( n.eq.170 ) then 
         do 90 i=1, n
            g(i) = zero
 90      continue
c
         do 100 i=1, 5
            g(34*i) = -8000.0d0  
 100     continue
c
         return
      end if
c
c     Construct right hand sides for n<>170 
c 
      g(1) = zero
      g(n) = zero
      fn2  = real(n-2)
c
      do 200 i=2, n-1
         g(i) = (real(i-1)/fn2)*1.0d2
 200  continue
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine perrhs ( n, b )
       
      integer  n
      double   precision b(n), drand
c
c***********************************************************************
c
c     This routine introduces a 5% perturbation with random sign 
c     in each component of the input vector b
c
c     PARAMETERS
c
c     N    integer variable specifying the dimension of the vector b
c
c     B    real array containing the rhs on input. On output, B 
c          holds the new rhs.
c 
c     MACHINE DEPENDENCIES
c
c     This routine may require a modification. It makes use of the non 
c     intrinsic function drand() to generate random numbers. Therefore, 
c     the user should  replace drand() by the appropriate function 
c     supported by the Fortran compiler.
c
c***********************************************************************
c
      integer  i, iexp, iexp2, nexp       
      double   precision    bi, delta
c
      do 10 i=1, n
         bi  = b(i)
         if ( bi.ne.0.0d0 ) then
            delta = bi*5.0d-2
c
c     WARNING, the function drand() may not be supported by the
c     Fortran 77 compiler. 
c
            iexp  = int(drand(i)*1.0d2)
c
            iexp2 = 2*(iexp/2)
            nexp  = 1
            if ( iexp.eq.iexp2 ) nexp = 0
            b(i) = bi + delta*(-1)**nexp
          end if
 10   continue
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine dsymvs ( n, a, adiag, jptr, indr, x, ax )
       
      integer  n, jptr(n+1),indr(*)
      double   precision a(*), adiag(n), x(n), ax(n)
c
c***********************************************************************
c
c     This subroutine performs the matrix-vector multiplication 
c     A*x = ax, where A is a sparse symmetric matrix with strict 
c     lower triangle stored in a sparse column format.
c
c     PARAMETERS
c
c     N   is an integer variable.
c         On entry n must specify the dimension of the system.
c         On exit n is unchanged.
c
c     A   is a real array of dimension *.  On entry a must contain 
c         the nonzero coefficients in the strict lower triangle of 
c         A stored in sparse column format.
c         On exit a is unchanged.
c
c     ADIAG is a real array of dimension n. On entry adiag must 
c         contain the diagonal elements of A. On exit adiag is 
c         unchanged.
c
c     JPTR  is an integer array of dimension n + 1. On entry jptr 
c         must contain pointers to the columns of A. The nonzeros 
c         in column j of A must be contained in the
c            jptr(j), ... , jptr(j+1) - 1 positions of a.
c         On exit jptr is unchanged.
c
c     INDR  is an integer array of dimension *.
c         On entry indr must contain row indices of the elements 
c         in a. On exit indr is unchanged.
c
c     X   is a real array of dimension n. On entry x must contain 
c         the vector x. On exit x is unchanged.
c
c     AX  is a real array of dimension n. On entry ax need not be 
c         specified. On exit ax contains the product A*x.
c
c
c***********************************************************************
c
      integer  i,j
      double   precision axt, xt, zero
      data     zero /0.0d0/
c
c     Multiplication by the diagonal.
c
      do 10 i = 1, n
         ax(i) = adiag(i)*x(i)
   10 continue
c
c     Multiplication by the strict lower and upper parts. 
c
      do 30 j = 1, n
         xt = x(j)
         axt = zero
         do 20 i = jptr(j), jptr(j+1) - 1
            axt = axt + a(i)*x(indr(i))
            ax(indr(i)) = ax(indr(i)) + a(i)*xt
   20    continue
         ax(j) = ax(j) + axt
   30 continue

      end
c
c-----------------------------------------------------------------------
c
      double precision function dnorma( n, a, adiag, jptr, indr, ax )
       
      integer  n, jptr(n+1),indr(*)
      double   precision a(*), adiag(n), ax(n)
c
c***********************************************************************
c
c     This function computes the infinite norm of a matrix A,
c     where A is a sparse symmetric matrix with strict lower triangle
c     stored in a sparse column format.
c
c     PARAMETERS
c
c     N   is an integer variable.
c         On entry n must specify the dimension of the system.
c         On exit n is unchanged.
c
c     A   is a real array of dimension *.  On entry a must contain 
c         the nonzero coefficients in the strict lower triangle of 
c         A stored in sparse column format.
c         On exit a is unchanged.
c
c     ADIAG is a real array of dimension n. On entry adiag must 
c         contain the diagonal elements of A. On exit adiag is 
c         unchanged.
c
c     JPTR  is an integer array of dimension n + 1. On entry jptr 
c         must contain pointers to the columns of A. The nonzeros 
c         in column j of A must be contained in the
c            jptr(j), ... , jptr(j+1) - 1 positions of a.
c         On exit jptr is unchanged.
c
c     INDR  is an integer array of dimension *.
c         On entry indr must contain row indices of the elements 
c         in a. On exit indr is unchanged.
c
c     AX    is a double precision array of dimension n.
c         On entry ax need not be specified.
c         On exit ax contains the product |A|*e.
c
c***********************************************************************
c
      integer  i, j
      double   precision axt, zero
      data     zero/0.0d0/
c
c     Multiplication by the diagonal.
c
      do 10 i = 1, n
         ax(i) = abs(adiag(i))
   10 continue
c
c     Multiplication by the strict lower and upper parts. 
c
      do 30 j = 1, n
         axt = zero
         do 20 i = jptr(j), jptr(j+1) - 1
            axt = axt + abs(a(i)) 
            ax(indr(i)) = ax(indr(i)) + abs(a(i)) 
   20    continue
         ax(j) = ax(j) + axt
   30 continue
c
c     Finding the maximum sum 
c
      axt = zero
      do 40 i=1, n
         if ( ax(i) .gt. axt ) axt = ax(i)
 40   continue
c
      dnorma = axt 

      end
c
c-----------------------------------------------------------------------
c
      subroutine sparse ( q, a, n, nz, adiag, jptra, indra )
       
      integer  n, nz,  jptra(n+1), indra(*)
      double   precision q(n,n), a(*), adiag(n)
c
c     Given a symmetric matrix A in dense format, this routine
c     computes the sparse structure of A.
c
c***********************************************************************
c
c     PARAMETERS
c
c     Q   is a real array of dimension nxn containing the coefficients
c         of the matrix in dense format
c
c     A   is a real array of dimension *.  On output a contains 
c         the nonzero coefficients in the strict lower triangle of 
c         A stored in sparse column format.
c
c     N   is an integer variable.
c         On entry n must specify the dimension of the system.
c         On exit n is unchanged.
c
c     NZ  is an integer variable. On exit it gives the number of
c         nonzero entries in the strict lower triangle of A.
c
c
c     ADIAG is a real array of dimension n. On exit  adiag 
c         contains the diagonal elements of A.
c
c     JPTRA is an integer array of dimension n + 1. On exit jptr 
c         contains pointers to the columns of A. The nonzeros 
c         in column j of A must be contained in the
c            jptr(j), ... , jptr(j+1) - 1 positions of a.
c
c     INDRA is an integer array of dimension *.
c         On exit indr contains row indices of the elements 
c         in a. 
c
c***********************************************************************
c
      integer  i, j
      double   precision qij
c
c     Checking symmetry & initializing the arrays
c     
      nz       = 1
      jptra(1) = 1
c
      do 10 i=1, n
         adiag(i) = q(i,i)
 10   continue
c
      do 30 i=1, n-1
         do 20 j=i+1, n
            qij = (q(i,j) + q(j,i))*0.5d0
            if ( qij.ne. 0.0d0 ) then
               a(nz)     = qij
               indra(nz) = j
               nz = nz + 1
            end if
 20      continue
         jptra(i+1) = nz
 30   continue
c
      jptra(n+1) = jptra(n)
      nz = nz - 1
c
      return
      end
c
c-----------------------------------------------------------------------
c
      double precision function dnrmif ( n, x )
       
      integer  n
      double   precision x(n)
c
c***********************************************************************
c
c     Computation of the infinity norm of a vector
c
c
c     PARAMETERS
c
c     N   is an integer variable specifying dimension of the problem
c  
c     X   is a double precision array containing the vector
c         whose norm will be computed
c
c***********************************************************************
c
      integer  i
      double   precision xnorm, dx, zero
      data     zero /0.0d0/
c
      xnorm = zero
      do 10 i=1, n
         dx = dabs(x(i))
         if ( dx.gt.xnorm ) xnorm = dx
 10   continue
c
      dnrmif = xnorm
c
      return
      end









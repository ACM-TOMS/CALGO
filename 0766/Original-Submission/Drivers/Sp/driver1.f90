!***********************************************************************
!     This is a test driver for the subroutines VECTOR_PADE,           *
!     invert_Sylv and solve_Sylv. The driver requires no input data.B  *
!***********************************************************************

      module working_area_main
!        Local variables for main.
         real, dimension (:,:,:), allocatable :: S_n, S_star_n
!        Assumed shape dummy arguments for VECTOR_PADE.
         integer, dimension (:),  allocatable :: n
         real, dimension (:),     allocatable :: gamma, gamma_star, kappa
         real, dimension (:,:),   allocatable :: A
         real, dimension (:,:,:), allocatable :: S, S_star
!        Assumed shape dummy arguments for solve_Sylv.
         real, dimension (:),     allocatable :: b, x
!        Assumed shape dummy arguments for invert_Sylv.
         real, dimension (:,:),   allocatable :: M, M_inv
      end module working_area_main

      module working_area_VECTOR_PADE
!        Automatic arrays for VECTOR_PADE.
!        Variables used to compute NPHS.
         real, dimension(:,:,:),  allocatable :: S_hat, New_S
         real, dimension(:,:),    allocatable :: T
!        Variables used to compute NSPS.
         real, dimension(:,:,:),  allocatable :: S_star_hat, New_S_star
         real, dimension(:,:,:),  allocatable :: T_star
      end module working_area_VECTOR_PADE

      module working_area_invert_Sylv
         real, dimension (:),     allocatable :: gamma, gamma_star,    &
                                                 gammap, kappa,        &
                                                 a0, a_norm, b, c
         real, dimension (:,:,:), allocatable :: S, S_star
      end module working_area_invert_Sylv

      module working_area_solve_Sylv
         real, dimension (:),     allocatable :: gamma, gamma_star,    &
                                                 gammap, kappa,        &
                                                 a_norm, c, d
         real, dimension (:,:,:), allocatable :: S, S_star
      end module working_area_solve_Sylv

      program main

      use working_area_main
      implicit none

      interface
         subroutine VECTOR_PADE(k, n, A, tau,                          &
               S, gamma, S_star, gamma_star, kappa, num_of_kappas, flag)
            integer,                    intent (in)    :: k
            integer, dimension (:),     intent (in)    :: n
            real,    dimension (:,:),   intent (in)    :: A
            real,                       intent (in)    :: tau
            real,    dimension (:,:,:), intent (out)   :: S, S_star
            real,    dimension (:),     intent (out)   :: gamma,   &
                                                          gamma_star, kappa
            integer,                    intent (out)   :: num_of_kappas
            integer,                    intent (inout) :: flag
         end subroutine VECTOR_PADE
         subroutine invert_Sylv(k, n, A, tau, Sylvinv, cond_number, flag)
            integer,                    intent(in)    :: k
            integer, dimension (:),     intent(in)    :: n
            real,    dimension (:,:),   intent(inout) :: A
            real,                       intent(in)    :: tau
            real,    dimension (:,:),   intent(out)   :: Sylvinv
            real,                       intent(out)   :: cond_number
            integer,                    intent(out)   :: flag
         end subroutine invert_Sylv
         subroutine solve_Sylv(k, n, A, b, tau, x, cond_number, flag)
            integer,                    intent(in)    :: k
            integer, dimension (:),     intent(in)    :: n
            real,    dimension (:,:),   intent(inout) :: A
            real,    dimension (:),     intent(in)    :: b
            real,                       intent(in)    :: tau
            real,    dimension (:),     intent(out)   :: x
            real,                       intent(out)   :: cond_number
            integer,                    intent(out)   :: flag
         end subroutine solve_Sylv
      end interface

      integer alpha, beta, i, j, k, l, offset, norm_n, num_of_kappas, &
              flag, S_deg, S_star_deg
      real tau, cond_number

!***********************************************************************
!                            Test Example                              *
!     Source: S. Cabay and A. Jones and G. Labahn, "Experiments with   *  
!     a weakly stable algorithm for computing Pade'-Hermite and        *   
!     simultaneous Pade' approximants", submitted to ACM Transactions  *
!     on mathematical software.                                        *   
!***********************************************************************

      k = 2
      allocate (n(0:k))
      n(0) = 2
      n(1) = 3
      n(2) = 1
      write(*,983) (n(beta),beta=0,k)

!     The norm of n is norm_n = n(0) + ... + n(k).
      norm_n = 6

!     S_deg = max{n(0), ..., n(k)} + 1, the largest degree of the 
!     polynomials in S.
      S_deg = 4

!     S_star_deg = norm_n - max{n(0), ..., n(k)} + 1, the largest degree 
!     of the polynomials in S_star.
      S_star_deg = 6

      allocate (A(0:norm_n, 0:k),                                     &
                S(0:S_deg,0:k,0:k), S_n(0:S_deg,0:k,0:k),             &
                S_star(0:norm_n,0:k,0:k), S_star_n(0:norm_n,0:k,0:k), &
                gamma(0:k), gamma_star(0:k), kappa(0:norm_n),         &
                b(0:norm_n), x(0:norm_n),                             &
                M(norm_n,norm_n), M_inv(norm_n,norm_n))

!     A is a vector of power series A = (A(0), ..., A(k)), where
!     A(beta) = A(0,beta) + A(1,beta)*z + ... + A(k,beta)*z**k.
      A(0,0) = 1
      A(1,0) = -1
      A(2,0) = 2
      A(3,0) = -2
      A(4,0) = 3
      A(5,0) = -3
      A(6,0) = 4

      A(0,1) = 0
      A(1,1) = 2
      A(2,1) = 0
      A(3,1) = 3
      A(4,1) = 0
      A(5,1) = 4
      A(6,1) = 0

      A(0,2) = -1
      A(1,2) = 1
      A(2,2) = 5
      A(3,2) = 3
      A(4,2) = 2
      A(5,2) = -2
      A(6,2) = -6

!     tau is the stability parameter described in subroutine VECTOR_PADE.
      tau = 100000.0

!     Print A.
      l = 0
      write(*,988) 
      write(*,986) (A(0,beta), beta=0, k), l
      do l=1, norm_n
         write(*,987) (A(l,beta), beta=0, k), l
      end do

!***********************************************************************
!     Obtain the scaled Pade'-Hermite system S and the scaled          *
!     simultaneous system S_star of types n for the vector of power    *  
!     series associated with A.                                        *   
!***********************************************************************
!     Since the sizes of the components of A affect the condition
!     numbers of the associated Sylvester matrices, and therefore the
!     choice of tau, it is usually advisable to scale the columns of A
!     (i.e., the power series associated with each column) before 
!     calling VECTOR_PADE.

      call VECTOR_PADE(k, n, A, tau,                                   &
               S, gamma, S_star, gamma_star, kappa, num_of_kappas, flag)

!     S_n is the normalized Pade'-Hermite system.
      do l=0, S_deg
         do beta=0, k
            do alpha=0, k
               if ((beta .eq. 0 .and. l .gt. n(alpha)+1) .or. &
                   (beta .gt. 0 .and. l .gt. n(alpha)  )) then
!                 Pad with zeroes to make the polynomial of degree 
!                 S_deg.
                  S_n(l,alpha,beta) = 0.0
               else
                  S_n(l,alpha,beta) = S(l,alpha,beta)/gamma(beta)
               endif
            end do
         end do
      end do
        
!     S_star_n is the normalized simultaneous Pade' system.
      do l=0, S_star_deg
         do beta=0, k
            do alpha=0, k
               if ((alpha .eq. 0 .and. l .gt. norm_n - n(beta)  ) .or. &
                   (alpha .gt. 0 .and. l .gt. norm_n - n(beta)+1)) then
!                 Pad with zeroes to make the polynomial of degree 
!                 S_star_deg.
                  S_star_n(l,alpha,beta) = 0.0
               else
                  S_star_n(l,alpha,beta) = S_star(l,alpha,beta)/ &
                                           gamma_star(alpha)
               endif
            end do
         end do
      end do

      write(*,990) 
      do l=0, S_deg
         write(*,991) l
         do alpha=0,k
!           Multiplication by 37 makes S_n integral.
            write(*,992) (37*S_n(l,alpha,beta), beta=0, k)
         end do
      end do

      write(*,984) 
      write(*,985) (kappa(i), i=0, num_of_kappas)

      write(*,993) 
      do l=0, S_star_deg
         write(*,995) l
         do alpha=0,k
!           Multiplication by 37 makes S_star_n integral.
            write(*,994) (37*S_star_n(l,alpha,beta), beta=0, k)
         end do
      end do


!***********************************************************************
!     Obtain the inverse of the striped Sylvester matrix M associated  *
!     with A.                                                          *  
!***********************************************************************

      call invert_Sylv(k, n, A, tau, M_inv, cond_number, flag)

!     Compute M (required for output purposes only).
      offset = 0
      do l=0, k
         do j=1, n(l)
            do i = 1, j-1
               M(i,offset+j) = 0
             end do
            do i=j, norm_n
               M(i,offset+j) = A(i-j,l)
             end do
          end do
         offset = offset + n(l)
        end do

!     Print M.
      write(*,996)
      do i = 1, norm_n
         write(*,997) (M(i,j),j=1,norm_n)
       end do

!     Print the inverse of M.
      write(*,998)
      do i = 1, norm_n
         write(*,980) (37*37*M_inv(i,j),j=1,norm_n)
       end do

!***********************************************************************
!     Solve the system Mx = b, using the subroutine solve_Sylv. The    *
!     vector b is the sum of the columns of M.                         *  
!***********************************************************************
      b(0) = 0
      b(1) = 3
      b(2) = 8
      b(3) = 8
      b(4) = 6
      b(5) = 5

      call solve_Sylv(k, n, A, b, tau, x, cond_number, flag)

!     Print b and x.
      write(*,999)
      write(*,981) (b(l),l=0,norm_n-1)
      write(*,982)
      write(*,981) (x(l),l=0,norm_n-1)

 980  format(2x,6f7.0)
 981  format(3x, 6f4.0)
 982  format(//6x, 'Solution x of Mx = b'/)
 983  format(//12x, 'n = (',i1,',',i1,',',i1,')'/)
 984  format(//3x,'Stability of points along the piecewise diagonal'/, &
             3x,'through n. Rough estimates of condition numbers'/,    &
             3x,'of the Sylvester matrices associated with these'/,    &
             3x,'points.')
 985  format(/5x, 8e9.2)
 986  format(2x,'   [', 3F4.0, ' ]  *  z**',I1)
 987  format(2x,'+  [', 3F4.0, ' ]  *  z**',I1)
 988  format(" Input vector of power series A")
 990  format(// "     Pade'-Hermite System S"/ &
      ' Normalized and multiplied by 37'/)
 991  format(" Coefficient matrix of z**",I1)
 992  format(2x, 3F7.0)
 993  format(//" Simultaneous Pade' System S_star"/ & 
             ' Normalized and multiplied by 37'/)
 994  format(5x, 3F7.0)
 995  format("    Coefficient matrix of z**",I1)
 996  format(//'  Striped Sylvester matrix M'/)
 997  format(3x, 6F4.0)
 998  format(//8x, 'Inverse of M multiplied by 37**2'/)
 999  format(//6x, 'Right-hand vector b'/)
         
      deallocate (n, A, S, S_n, S_star, S_star_n, gamma, gamma_star, & 
                                                kappa, b, x, M, M_inv)
      stop
      end program main

      subroutine invert_Sylv(k, n, A, tau, Sylvinv, cond_number, flag)

!***********************************************************************
!                                                                      *
!     For the vector of integers                                       *
!                                                                      *
!                      n = [n(0),...,n(k)],                            *
!                                                                      *
!     let                                                              *
!                                                                      *
!                  ||n|| = n(0)+...+n(k).                              *
!                                                                      *
!     Define the order ||n|| striped Sylvester matrix M to be          *
!                                                                      *
!                                                                      *
!                                                                      *
!  |   A(0,0)                     |   |   A(0,k)                     | *
!  |          .                   |   |          .                   | *
!  |              .               |   |              .               | *
!  |     .            .           |   |     .            .           | *
!  |     .              A(0,0)    |...|     .              A(0,k)    |.*
!  |     .                .       |   |     .                .       | *
!  |                      .       |   |                      .       | *
!  |                      .       |   |                      .       | *
!  |A(||n||-1,0)...A(||n||-n(0),0)|   |A(||n||-1,k)...A(||n||-n(k),k)| *
!                                                                      *
!     This subroutine computes the inverse, Sylvinv, of M and gives    *
!     a rough estimate, cond_number, of the algorithm.                 *
!                                                                      *
!     The inverse is obtained by using the formula derived in          *
!     S. Cabay and A. Jones and G. Labahn, "Computation of Numerical   *
!     Pade'-Hermite and Simultaneous Pade' Systems I: Near Inversion   *
!     of Generalized Sylvester Matrices", SIAM journal on matrix       *
!     analysis and applications, 17 (1996), 248-267.                   *
!                                                                      *
!     The formula expresses the inverse in terms of Pade'-Hermite and  *
!     Simultaneous Pade' Systems of type n for the vector of power     *
!     series that can be associated with A. These systems are computed *
!     by the subroutine VECTOR_PADE.                                   *
!                                                                      *
!***********************************************************************
!                                                                      *
!    on entry                                                          *
!       k            integer                                           *
!                    There are k+1 stripes in M.                       *
!                                                                      *
!       n            integer (0:k)                                     *
!                    The beta'th strip in M has n(beta) columns.       *
!                                                                      *
!       A            real (0:sum(n), 0:k)                              *
!                    Each column of this matrix gives one stripe in    *
!                    the striped Sylvester matrix M.                   *
!                                                                      *
!       tau          real                                              *
!                    Stability parameter required by the subroutine    *
!                    VECTOR_PADE. Very roughly speaking, the residual  *
!                    error, b - M*x, will look like                    *
!                               tau * unit-error * ||A||.              *
!                    For most efficient computation, tau should be     *
!                    chosen as large as the lack of accuracy will      *
!                    permit.                                           *
!                                                                      *
!    on exit                                                           *
!       Sylvinv      real (sum(n), sum(n))                             *
!                    The inverse  of M.                                *
!                                                                      *
!       cond_number  real                                              *
!                    Very rough estimate of the condition number of    *
!                    this algorithm.                                   *
!                                                                      *
!       flag       integer                                             *
!                    Error parameter.                                  *
!                    flag = 0, no errors                               *
!                    flag = 1, the Sylvester matrix at the point n     *
!                              is ill-conditioned; i.e.,               *
!                              cond_number >= tau.                     *
!                    flag = 2, the Sylvester matrix at the point n     *
!                              is numerically singular.                *
!                    flag = 3, input variables are incorrect.          *
!                                                                      *
!                                                                      *
!    functions and subroutines                                         *
!       divide_series     Divides one power series by another.         *
!       VECTOR_PADE       Computes scaled Pade'-Hermite and            *
!                         simultaneous Pade' systems for the power     *
!                         series associated with A.                    *
!                                                                      *
!***********************************************************************
 
      use working_area_invert_Sylv
      implicit none

      interface
         subroutine VECTOR_PADE(k, n, A, tau,                          &
               S, gamma, S_star, gamma_star, kappa, num_of_kappas, flag)
            integer,                    intent (in)    :: k
            integer, dimension (:),     intent (in)    :: n
            real,    dimension (:,:),   intent (in)    :: A
            real,                       intent (in)    :: tau
            real,    dimension (:,:,:), intent (out)   :: S, S_star
            real,    dimension (:),     intent (out)   :: gamma,   &
                                                          gamma_star, kappa
            integer,                    intent (out)   :: num_of_kappas
            integer,                    intent (inout) :: flag
         end subroutine VECTOR_PADE
      end interface

!     invert_Sylv subroutine parameters.
      integer,                    intent(in)    :: k
      integer, dimension (0:),    intent(in)    :: n
      real,    dimension (0:,0:), intent(inout) :: A
      real,                       intent(in)    :: tau
      real,    dimension (1:,1:), intent(out)   :: Sylvinv
      real,                       intent(out)   :: cond_number
      integer,                    intent(out)   :: flag
      
!     Local variables.
      integer alpha, beta, i, j, l,  offset, norm_n, num_of_kappas 
      real norm_a0_inv
      allocate (gamma(0:k), gamma_star(0:k), gammap(0:k), a_norm(0:k), &
                kappa(0:sum(n)), a0(0:sum(n)), b(0:sum(n)),c(0:sum(n)),&
                S(0:maxval(n)+1, 0:k, 0:k), S_star(0:sum(n), 0:k, 0:k))

      norm_n = sum(n)
        
!     Check the validity of input parameters.
      if (     k .lt. 1                    .or. &
               k .gt. size(n) - 1          .or. &
               0 .gt. minval(n)            .or. &
          norm_n .gt. size(A(:,0)) - 1     .or. &
               k .gt. size(A(0,:)) - 1     .or. &
             0.0 .eq. A(0,0)               .or. &
          norm_n .gt. size(Sylvinv(:,1))   .or. &
          norm_n .gt. size(Sylvinv(1,:)))  then
         flag = 3
         return
      else
         flag = 0
      endif

!     Compute the Pade'-Hermite system S and simultaneous Pade' system
!     S_star of type n for the vector of power series associated with A.
!     VECTOR_PADE requires values for the coefficients of z**norm_n in A, 
!     so temporarily assign some.
      do beta=0, k
         a_norm(beta) = A(norm_n,beta)
         A(norm_n,beta) = 0
      end do
      call VECTOR_PADE(k, n, A, tau,                                   &
               S, gamma, S_star, gamma_star, kappa, num_of_kappas, flag)
      do beta=0, k
         A(norm_n,beta) = a_norm(beta) 
      end do

      if (flag .eq. 0 .or. flag .eq. 1) then
!     VECTOR_PADE successfully computed systems of type n.

!        Some initializations.
         do j = 0, norm_n-1
            a0(j) = A(j,0)/A(0,0)
         end do
         do i = 1, norm_n
            do j = 1, norm_n
               Sylvinv(i,j) = 0.0
            end do
         end do

!        A rough estimate of the condition number of this algorithm is 
!        given in the last entry in kappa multiplied by the norm of the 
!        inverse of a0.
         b(0) = 1
         do l=1, norm_n-1
            b(l) = 0
         end do
         call divide_series(a0, b, c, norm_n-1, flag)
         norm_a0_inv = abs(c(0))
         do l=1, norm_n-1
            norm_a0_inv = norm_a0_inv + abs(c(l))
         end do
         cond_number = kappa(num_of_kappas) * norm_a0_inv

!        Now, evaluate the inverse formula.

         do beta=0, k
            gammap(beta) = gamma(beta) * gamma_star(beta)
         end do

         do beta = 0, k
!           If beta = 0, b corresponds with the power series v_star(z); 
!           otherwise, b corresponds with the power series 
!           z*q_star(beta).
            if (beta .eq. 0) then
               do l = 0, norm_n-1
                  if (l .le. norm_n-n(0)) then
                     b(l) = S_star(l,0,0)
                  else 
                     b(l) = 0.0
                  endif
               end do
            else
               do l = 0, norm_n-1
                  if (l .le. norm_n-n(0)) then
                     b(l) = S_star(l+1,beta,0)
                  else 
                     b(l) = 0.0
                  endif
               end do
            endif

!           c corresponds with the power series obtained by dividing b 
!           with the inverse of the power series associated with the 
!           first column of A.
            call divide_series(a0, b, c, norm_n-1, flag)
 
!           Multiply c on the left by the coefficients of the 
!           polynomials of the beta'th column of S.
            offset = 0
            do alpha = 0, k
               do i = 0, n(alpha)-1
                  if (beta .eq. 0) then 
                     j=i-1 
                  else 
                     j=i 
                  endif
                  do l = 0, norm_n - 1
                     Sylvinv(offset+n(alpha)-i,norm_n-l)              &
                      = Sylvinv(offset+n(alpha)-i, norm_n-l)          &
                      + S(n(alpha)-j,alpha,beta) * c(l) / gammap(beta)
                  end do
               end do
               offset = offset + n(alpha)
            end do
         end do

!        Sum over the coefficients of the polynomials in S. 
         offset = 0
         do alpha = 0, k
            do i = 1, n(alpha)-1
               do l = 1, norm_n-1
                  Sylvinv(offset+n(alpha)-i,norm_n-l)                  &
                               = Sylvinv(offset+n(alpha)-i,norm_n-l)   &
                               + Sylvinv(offset+n(alpha)+1-i,norm_n-l+1) 
               end do
            end do
            offset = offset + n(alpha)
         end do

      endif
      deallocate (gamma, gamma_star, gammap, a_norm, kappa, a0, b, c, &
                                                             S, S_star)
      return
      end subroutine invert_Sylv


      subroutine solve_Sylv(k, n, A, b, tau, x, cond_number, flag)
      
!***********************************************************************
!                                                                      *
!     For the vector of integers                                       *
!                                                                      *
!                      n = [n(0),...,n(k)],                            *
!                                                                      *
!     let                                                              *
!                                                                      *
!                  ||n|| = n(0)+...+n(k).                              *
!                                                                      *
!     Define the order ||n|| striped Sylvester matrix M to be          *
!                                                                      *
!                                                                      *
!                                                                      *
!  |   A(0,0)                     |   |   A(0,k)                     | *
!  |          .                   |   |          .                   | *
!  |              .               |   |              .               | *
!  |     .            .           |   |     .            .           | *
!  |     .              A(0,0)    |...|     .              A(0,k)    |.*
!  |     .                .       |   |     .                .       | *
!  |                      .       |   |                      .       | *
!  |                      .       |   |                      .       | *
!  |A(||n||-1,0)...A(||n||-n(0),0)|   |A(||n||-1,k)...A(||n||-n(k),k)| * 
!                                                                      *
!     Given the vector b, solve_Sylv determines the solution x of the  *
!     linear system of equations                                       *
!                                                                      *
!                  M * x = b.                                          *
!                                                                      *
!     The solution is obtained by using the inverse formula for M      *
!     derived in S. Cabay and A. Jones and G. Labahn, "Computation of  *
!     Numerical Pade'-Hermite and Simultaneous Pade' Systems I: Near   *
!     Inversion of Generalized Sylvester Matrices", SIAM journal on    *
!     matrix analysis and applications, 17 (1996), 248-267.            *
!                                                                      *
!     The formula expresses the inverse in terms of Pade'-Hermite and  *
!     Simultaneous Pade' Systems of type n for the vector of power     *
!     series that can be associated with A. These systems are computed *
!     by the subroutine VECTOR_PADE.                                   *
!                                                                      *
!***********************************************************************
!                                                                      *
!    on entry                                                          *
!       k            integer                                           *
!                    There are k+1 stripes in M.                       *
!                                                                      *
!       n            integer (0:k)                                     *
!                    The beta'th strip in M has n(beta) columns.       *
!                                                                      *
!       A            real (0:sum(n), 0:k)                              *
!                    Each column of this matrix gives one stripe in    *
!                    the striped Sylvester matrix M.                   *
!                                                                      *
!       b            real (0:sum(n))                                   *
!                    The right-hand vector.                            *
!                                                                      *
!       tau          real                                              *
!                    Stability parameter required by the subroutine    *
!                    VECTOR_PADE. Very roughly speaking, the residual  *
!                    error, b - M*x, will look like                    *
!                               tau * unit-error * ||A||.              *
!                    For most efficient computation, tau should be     *
!                    chosen as large as the lack of accuracy will      *
!                    permit.                                           *
!                                                                      *
!    on exit                                                           *
!       x            real (0:sum(n))                                   *
!                    The solution vector.                              *
!                                                                      *
!       cond_number  real                                              *
!                    Very rough estimate of the condition number of    *
!                    this algorithm.                                   *
!                                                                      *
!       flag       integer                                             *
!                    Error parameter.                                  *
!                    flag = 0, no errors                               *
!                    flag = 1, the Sylvester matrix at the point n     *
!                              is ill-conditioned; i.e.,               *
!                              cond_number >= tau.                     *
!                    flag = 2, the Sylvester matrix at the point n     *
!                              is numerically singular.                *
!                    flag = 3, input variables are incorrect.          *
!                                                                      *
!    functions and subroutines                                         *
!       divide_series     Divides one power series by another.         *
!       VECTOR_PADE       Computes scaled Pade'-Hermite and            *
!                         simultaneous Pade' systems for the power     *
!                         series associated with A.                    *
!                                                                      *
!***********************************************************************
 
      use working_area_solve_Sylv
      implicit none

      interface
         subroutine VECTOR_PADE(k, n, A, tau,                          &
               S, gamma, S_star, gamma_star, kappa, num_of_kappas, flag)
            integer,                    intent (in)    :: k
            integer, dimension (:),     intent (in)    :: n
            real,    dimension (:,:),   intent (in)    :: A
            real,                       intent (in)    :: tau
            real,    dimension (:,:,:), intent (out)   :: S, S_star
            real,    dimension (:),     intent (out)   :: gamma,   &
                                                          gamma_star, kappa
            integer,                    intent (out)   :: num_of_kappas
            integer,                    intent (inout) :: flag
         end subroutine VECTOR_PADE
      end interface

!     solve_Sylv subroutine parameters.
      integer,                    intent(in)    :: k
      integer, dimension (0:),    intent(in)    :: n
      real,    dimension (0:,0:), intent(inout) :: A
      real,    dimension (0:),    intent(in)    :: b
      real,                       intent(in)    :: tau
      real,    dimension (0:),    intent(out)   :: x
      real,                       intent(out)   :: cond_number
      integer,                    intent(out)   :: flag

!     Local variables.
      integer i, l, num_of_kappas, offset, alpha, beta, norm_n, maxn, ni
      real    temp, norm_c
      allocate (gamma(0:k), gamma_star(0:k), gammap(0:k), a_norm(0:k), &
                kappa(0:sum(n)), c(0:sum(n)), d(0:sum(n)),             &
                S(0:maxval(n)+1, 0:k, 0:k), S_star(0:sum(n), 0:k, 0:k))

      norm_n = sum(n)

!     Check the validity of input parameters.
      if (     k .lt. 1                    .or. &
               k .gt. size(n) - 1          .or. &
               0 .gt. minval(n)            .or. &
          norm_n .gt. size(b)              .or. &
          norm_n .gt. size(x)              .or. &
          norm_n .gt. size(A(:,0)) - 1     .or. &
               k .gt. size(A(0,:)) - 1     .or. &
             0.0 .eq. A(0,0))              then
         flag = 3
         return
      else
         flag = 0
      endif

      maxn = maxval(n)
 
!     Compute the Pade'-Hermite system S and simultaneous Pade' system
!     S_star of type n for the vector of power series associated with A.
!     VECTOR_PADE requires values for the coefficients of z**norm_n in A,
!     so temporarily assign some.
      do beta=0, k
         a_norm(beta) = A(norm_n,beta)
         A(norm_n,beta) = 0
      end do
      call VECTOR_PADE(k, n, A, tau,                                   &
               S, gamma, S_star, gamma_star, kappa, num_of_kappas, flag)
      do beta=0, k
         A(norm_n,beta) = a_norm(beta) 
      end do

      if (flag .eq. 0  .or. flag .eq. 1) then
!     VECTOR_PADE has computed systems of type n.

         do l = 0, norm_n-1
            c(l) = A(l,0)/A(0,0)
         end do
!        A rough estimate of the condition number of the algorithm is 
!        given in the last entry in kappa multiplied by the norm of 
!        the inverse of c.
         d(0) = 1
         do l=1, norm_n-1
            d(l) = 0
         end do
         call divide_series(c, d, x, norm_n-1, flag)
         norm_c = abs(x(0))
         do l=1, norm_n-1
            norm_c = norm_c + abs(x(l))
         end do
         cond_number = kappa(num_of_kappas) * norm_c

         do beta=0, k
            gammap(beta) = gamma(beta) * gamma_star(beta)
         end do
 
!        Divide the power series associated with b by the inverse of 
!        the power series associated with the first column of A 
!        to get the power series d.                                
         call divide_series(c, b, d, norm_n-1, flag)

!        Apply the inverse formula.
         do i=0, norm_n-1
            x(i) = 0.0
         end do
         do beta=0, k
            do i = 1, maxn
               temp = 0.0
               ni = max(n(0),i)
               do l = 0, norm_n - ni 
                  if (beta .eq. 0) then
                     temp = temp + S_star(l,beta,0)*d(norm_n-i-l)
                  else
                     temp = temp + S_star(l+1,beta,0)*d(norm_n-i-l)
                  endif
               end do
               c(i-1) = temp/gammap(beta)
            end do
            offset = 0
            do alpha = 0, k
               do i = 1, n(alpha)
                  temp = 0.0
                  do l = 0, n(alpha)-i
                  if (beta .eq. 0) then
                     temp = temp + S(i+l+1,alpha,beta) * c(l)
                  else
                     temp = temp + S(i+l,alpha,beta) * c(l)
                  endif
                  end do
                  x(offset+i-1) = x(offset+i-1) + temp
               end do
               offset = offset + n(alpha)
            end do
         end do

      endif
      deallocate (gamma, gamma_star, gammap, a_norm, kappa, c, d, S, &
                                                               S_star)
      return
      end subroutine solve_Sylv

      subroutine divide_series(a, b, c, N, flag)

      implicit none
      integer N, flag
      real a(0:N), b(0:N), c(0:N)

!***********************************************************************
!     For the two power represented by a and b, with a(0) nonzero,     *
!     this subroutine divides b by a (modulo N) and stores the result  *
!     in c.                                                            *
!***********************************************************************

!     Local variables.
      integer i, j
      real temp

!     Check validity of input parameters.
      if (N .lt. 0 .or. a(0) .eq. 0) then
         flag = 1
         return
      else
         flag = 0
      endif

!     Solve a triangular system of equations.
      c(0) = b(0)/a(0)
      do i=1, N
         temp = b(i)
         do j=0, i-1
            temp = temp - c(j)*a(i-j)
         end do
         c(i) = temp/a(0)
         end do
      return
      end subroutine divide_series



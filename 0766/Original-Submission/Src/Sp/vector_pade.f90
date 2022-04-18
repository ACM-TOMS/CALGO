
      subroutine VECTOR_PADE(k, n, A, tau,                             &
                      S, gamma, S_star, gamma_star, kappa, last_i, flag)

!***********************************************************************
!                                                                      *
!     For the vector of integers                                       *
!                                                                      * 
!                      n = [n(0),...,n(k)],                            * 
!                                                                      *  
!     define                                                           * 
!                                                                      * 
!                  ||n|| = n(0)+...+n(k).                              *  
!                                                                      * 
!     Let A be the vector of k+1 power series                          *      
!                                                                      *      
!                      A = [A(0),...,A(k)],                            *
!                                                                      *     
!     where, for j=0,...,k, the first ||n||+1 terms of A(j) are given  *
!     on entry by                                                      *
!                                                                      *
!          A(j) = A(0,j) + A(1,j) * z + ... + A(||n||,j) * z**||n||.   * 
!                                                                      *
!     This subroutine VECTOR_PADE computes the scaled Pade'-Hermite    * 
!     system S of type n satisfying                                    *   
!                                                                      *
!                 A * S = T * z**(||n||+1) + delta_T                   *
!                                                                      *
!     and the scaled simultaneous Pade' system S_star of type n        *
!     satisfying                                                       *    
!                                                                      *
!               |-A(1) -A(2) ... -A(k)|                                * 
!               | A(0)                |                                *
!     S_star *  |       A(0)          | = T_star * z**(||n||+1)        * 
!               |            ...      |          + delta_T_star.       *    
!               |                 A(0)|                                *
!                                                                      *
!     The system S is a (k+1) x (k+1) matrix of polynomials            * 
!                                                                      *
!                             |S(0,0) ... S(0,k)|                      * 
!                        S =  |       ...       |,                     *       
!                             |S(k,0) ... S(k,k)|                      * 
!                                                                      *
!     where, for i=0,...,k and j=1,...,k,                              * 
!                                                                      *
!       S(i,j) = S(0,i,j) + S(1,i,j) * z + ... + S(n(i),i,j) * z**n(i) * 
!                                                                      *
!     and, for i=0,...,k,                                              * 
!                                                                      *
!       S(i,0) = S(0,i,0) + S(1,i,0) * z + ...                         * 
!                    ...  + S(n(i)+1,i,0) * z**(n(i)+1)                *
!     with S(0,i,0) = S(1,i,0) = 0.                                    *   
!                                                                      *
!     The system S_star is a (k+1) x (k+1) matrix of polynomials       * 
!                                                                      *
!                             |S_star(0,0) ... S_star(0,k)|            *   
!                   S_star =  |            ...            |,           *   
!                             |S_star(k,0) ... S_star(k,k)|            *      
!                                                                      *
!     where, for j=0,...,k,                                            *
!                                                                      *
!          S_star(0,j) = S_star(0,0,j)                                 * 
!                      + S_star(1,0,j) * z                             * 
!                      +                                               *   
!                      .                                               *    
!                      .                                               *    
!                      .                                               *  
!                      +                                               *     
!                      + S_star(||n||-n(j),0,j) * z**(||n||-n(j)),     *     
!                                                                      *
!     and, for i=1,...,k and j=0,...,k,                                * 
!                                                                      *
!          S_star(i,j) = S_star(0,i,j)                                 * 
!                      + S_star(1,i,j) * z                             * 
!                      +                                               *     
!                      .                                               *    
!                      .                                               *    
!                      .                                               *  
!                      +                                               *     
!                      + S_star(||n||-n(j)+1,i,j) * z**(||n||-n(j)+1), * 
!                                                                      *
!     with S_star(0,i,j) = S_star(1,i,j) = 0.                          * 
!                                                                      *
!     On exit, the residual errors very roughly satisfy                *
!                                                                      *
!         ||delta_T||, ||delta_T_star|| = tau * unit_error * ||A||,    * 
!                                                                      *
!     whereas, the relative errors in S and S_star are very roughly    * 
!     of size                                                          *
!                                                                      *
!       kappa(last_i) * tau * unit_error * ||A|| * ||A(0)^{-1}||,      *
!                                                                      *
!     where kappa(last_i) is a crude estimate of the condition         *
!     number of the associated striped Sylvester matrix at the point n.*
!                                                                      *
!     The above is to serve as a rough guide in selecting the tolerance* 
!     "tau". But, there is a trade-off between the accuracy of the     * 
!     results and the speed of the algorithm; the smaller the value of * 
!     tau, the costlier the computations can become. For most efficient*  
!     computation, tau should be chosen to be as large as the lack of  * 
!     accuracy will permit. For some problems, the best strategy for   * 
!     selecting tau may be first to run a trial in order to examine the* 
!     estimates kappa(i), i = 1, ..., last_i, of the condition numbers *
!     of the striped Sylvester matrices associated with intermediate   *
!     points along the diagonal through n. The algorithm accepts the   *
!     i'th point as "stable" if                                        *
!                                                                      *
!                           kappa(i) < tau.                            *
!                                                                      *
!     Note that the first column of S always yields a Pade'-Hermite    *
!     approximant of type (n(0)-1, ..., n(k)-1) for A, with a "small"  *
!     residual error delta_T. This is so even when the Sylvester matrix*
!     at the point n is ill-conditioned or singular (when singular,    *
!     however, the remaining columns of S may be meaningless).         *
!     Similarly, the first row of S_star always gives a simultaneous   *
!     Pade' approximant of type n with a "small" residual error        *
!     delta_T_star.                                                    *
!                                                                      *
!***********************************************************************
!                                                                      *
!    on entry                                                          *
!       A            real (0:sum(n), 0:k)                              *
!                    Vector of power series with A(0,0) nonzero.       *
!                                                                      *
!       k            integer                                           *
!                    There are k+1 power series in A.                  *
!                                                                      *
!       n            integer (0:k)                                     *
!                    Vector of degrees defining the type of S & S_star.*
!                                                                      *
!       tau          real                                              *
!                    Stability parameter. An intermediate solution is  *
!                    accepted if at that point kappa(i) < tau.         *
!                                                                      *
!    on exit                                                           *
!       S            real (0:maxval(n)+1, 0:k, 0:k),                   *
!                    Scaled Pade'-Hermite system of type n.            *
!                                                                      *
!       gamma        real (0:k)                                        *
!                    Scaling factors. To obtain the normalised NPHS    *
!                    divide the j'th column of S by gamma(j).          *
!                                                                      *
!       S_star       real (0:sum(n), 0:k, 0:k),                        *
!                    Scaled simultaneous Pade' system of type n.       *
!                                                                      *
!       gamma_star   real (0:k)                                        *
!                    Scaling factors. To obtain the normalised NSPS    *
!                    divide the i'th row of S_star by gamma_star(i).   *
!                                                                      *
!       kappa        real (0:sum(n)),                                  *
!                    kappa(i) is the estimate of condition number of   *
!                    the associated Sylvester submatrix at the i'th    *
!                    step.                                             *
!                                                                      *
!       last_i       integer                                           *
!                    The last step, corresponding to the point n.      *
!                                                                      *
!       flag       integer                                             *
!                    Error parameter.                                  *
!                    flag = 0, no errors                               *
!                    flag = 1, the Sylvester matrix at the point n     *
!                              is ill-conditioned; i.e.,               *
!                              kappa(last_i) >= tau.                   *
!                    flag = 2, the Sylvester matrix at the point n     *
!                              is numerically singular. The first      *
!                              column of S still yields a Pade'-Hermite*
!                              approximant of type (n(0)-1,...,n(k)-1) *
!                              and the first row of S_star still yields*
!                              a simultaneous Pade' approximant of     *
!                              type n. The remaining rows and columns  *
!                              are meaningless.                        *
!                    flag = 3, input variables are incorrect.          *
!                                                                      *
!    Note that the storage allocated to the subroutine array parameters*
!    can be larger than limits designated above, but never smaller.    *
!                                                                      *
!    functions and subroutines                                         *
!       build_T           Finds the residual for the NPHS.             *
!       build_T_star      Finds the residual for the NSPS.             *
!       build_S           Builds the NPHS by solving Sylvester systems.* 
!       build_S_star      Builds the NSPS by solving Sylvester systems.*
!       mult_S            Multiples two NPHS's.                        *
!       mult_S_star       Multiples two NSPS's.                        *
!       scale_S           Scales the NPHS and determines gamma.        *
!       scale_S_star      Scales the NSPS and determines gamma_star.   *
!       gen_next_vector   Generates the next point along a diagonal.   *
!                                                                      *
!     Note also that build_S and build_S_star call                     *
!       sgefa             Linpack routine to triangulate a matrix.     *
!       sgesi             Linpack routine to solve linear equations.   *
!                                                                      *
!                                                                      *
!***********************************************************************
!     The algorithm VECTOR_PADE is shown to be weakly stable in        *
!     S. Cabay and A. Jones and G. Labahn, "Computation of Numerical   *
!     Pade'-Hermite and Simultaneous Pade' Systems II: A Weakly Stable *
!     Algorithm", SIAM journal on matrix analysis and applications,    *
!     17 (1996), 268-297.                                              *
!***********************************************************************

      use working_area_VECTOR_PADE
      implicit none

      interface
         subroutine build_T(A, S, k, m, prevnorm_nus, norm_nus, T)
            integer,                    intent (in)    :: k, norm_nus, &
                                                          prevnorm_nus 
            integer, dimension (:),     intent (in)    :: m
            real,    dimension (:,:),   intent (in)    :: A
            real,    dimension (:,:,:), intent (in)    :: S
            real,    dimension (:,:),   intent (inout) :: T
         end subroutine build_T
         subroutine build_T_star(A, S_star, k, m, prevnorm_nus, &
                                 norm_nus, T_star)
            integer,                    intent (in)    :: k, norm_nus, &
                                                          prevnorm_nus 
            integer, dimension (:),     intent (in)    :: m
            real,    dimension (:,:),   intent (in)    :: A
            real,    dimension (:,:,:), intent (in)    :: S_star
            real,    dimension (:,:,:), intent (inout) :: T_star
         end subroutine build_T_star
         subroutine build_S(A, k, n, norm_n, S, singular)
            integer,                    intent (in)  :: k, norm_n
            integer, dimension (:),     intent (in)  :: n
            real,    dimension (:,:),   intent (in)  :: A
            real,    dimension (:,:,:), intent (out) :: S
            logical,                    intent (out) :: singular
         end subroutine build_S
         subroutine build_S_star(A_star, k, n, knorm_n, S_star,singular)
            integer,                    intent (in)  :: k, knorm_n
            integer, dimension (:),     intent (in)  :: n
            real,    dimension (:,:,:), intent (in)  :: A_star
            real,    dimension (:,:,:), intent (out) :: S_star
            logical,                    intent (out) :: singular
         end subroutine build_S_star
         subroutine mult_S(S, m, S_hat, nus, New_S, new_m, k)
            integer,                    intent (in)  :: k
            integer, dimension (:),     intent (in)  :: m, nus, new_m
            real,    dimension (:,:,:), intent (in)  :: S, S_hat
            real,    dimension (:,:,:), intent (out) :: New_S
         end subroutine mult_S
         subroutine mult_S_star(S_star, m, S_star_hat, nus, &
                                        New_S_star, new_m, k)
            integer,                    intent (in)  :: k
            integer, dimension (:),     intent (in)  :: m, nus, new_m
            real,    dimension (:,:,:), intent (in)  :: S_star, S_star_hat
            real,    dimension (:,:,:), intent (out) :: New_S_star
         end subroutine mult_S_star
         subroutine scale_S(S, m, k, gamma)
            integer,                    intent (in)    :: k
            integer, dimension (:),     intent (in)    :: m
            real,    dimension (:,:,:), intent (inout) :: S
            real,    dimension (:),     intent (out)   :: gamma
         end subroutine scale_S
         subroutine scale_S_star(S_star, m, k, gamma_star)
            integer,                    intent (in)    :: k
            integer, dimension (:),     intent (in)    :: m
            real,    dimension (:,:,:), intent (inout) :: S_star
            real,    dimension (:),     intent (out)   :: gamma_star
         end subroutine scale_S_star
      end interface

!     VECTOR_PADE subroutine parameters.
      integer,                       intent (in)    :: k
      integer, dimension (0:),       intent (in)    :: n
      real,    dimension (0:,0:),    intent (in)    :: A
      real,                          intent (in)    :: tau
      real,    dimension (0:,0:,0:), intent (out)   :: S, S_star
      real,    dimension (0:),       intent (out)   :: gamma,   &
                                                       gamma_star, kappa
      integer,                       intent (out)   :: last_i
      integer,                       intent (inout) :: flag

!     Local Variables.
      logical singular
      integer alpha, beta, i, l, step, norm_n, knorm_nus,         &
              m(0:k), norm_m, new_m(0:k), norm_new_m, nu(0:k),    &
              nus(0:k), norm_nus, prevnorm_nus, S_deg, S_star_deg
!     Variables used to compute NPHS.
      allocate (S_hat(0:maxval(n)+1, 0:k, 0:k), &
                New_S(0:maxval(n)+1, 0:k, 0:k), &
                T(0:sum(n), 0:k))
!     Variables used to compute NSPS.
      allocate (S_star_hat(0:sum(n), 0:k, 0:k), &
                New_S_star(0:sum(n), 0:k, 0:k), &
                T_star(0:sum(n), 0:k, 1:k))

      norm_n = sum(n)

!     Check the validity of input parameters.
      flag = 0
         
      if (     k .lt. 1                       .or. &
               k .gt. size(n) - 1             .or. &
          norm_n .gt. size(A(:,0)) - 1        .or. &
               k .gt. size(A(0,:)) - 1        .or. &
             0.0 .eq. A(0,0)                  .or. &
       maxval(n) .gt. size(S(:,0,0)) - 2      .or. &
               k .gt. size(S(0,:,0)) - 1      .or. &
               k .gt. size(S(0,0,:)) - 1      .or. &
          norm_n .gt. size(S_star(:,0,0)) - 1 .or. &
               k .gt. size(S_star(0,:,0)) - 1 .or. &
               k .gt. size(S_star(0,0,:)) - 1) flag = 3

      if (flag .ne. 0) then
         return
      else

!        The initial stable point is m = (-1,0,...,0).

         step = 1
         m(0) = -1
         do beta=1, k
            m(beta) = 0
         end do
         kappa(0) = 1.0

!        At the initial point m, set S and S_star to be identity matrices.
         do alpha=0, k
            do beta=0, k
               if (alpha .eq. beta) then
                  S(0,alpha,beta) = 1.0
                  S_star(0,alpha,beta) = 1.0
               else
                  S(0,alpha,beta) = 0.0
                  S_star(0,alpha,beta) = 0.0
               endif
            end do
               S(1,alpha,0) = 0.0
               S_star(1,alpha,0) = 0.0
         end do

!        The index i references the i'th point along the path from the 
!        initial point (-1,0,...,0) with index i=0 to the final point 
!        n with index i=last_i. 
         i = 0
         last_i = 0
         do beta=1,k
            last_i = max(last_i, n(beta))
         end do
         last_i = min(n(0), last_i) + 1


         do while (i.lt.last_i .and. flag.eq.0)
!        Main loop 
!        ***************************************************************
!        * The iteration moves from one stable point m with index i to * 
!        * the next stable point until the last point n with index     * 
!        * last_i has been reached.                                    * 
!        ***************************************************************

            norm_m = sum(m)

!           nu is the difference between the last point n and the 
!           current stable point m.
            nu(0) = n(0) - m(0) - 1
            do beta=1, k
               nu(beta) = n(beta) - m(beta)
            end do

            step = 0
            flag = 1
            prevnorm_nus = -1

            do while (i+step.lt.last_i .and. flag.gt.0)
!           Inner loop
!           ************************************************************
!           * Given a stable point m with index i, examine successive  *
!           * points along the diagonal with indices i+step, step=1,...*
!           * until a stable one is found or until the last point n has*
!           * been examined.                                           *
!           ************************************************************

               step = step + 1
               flag = 1

!              *********************************************************
!              * Compute nus. The objective is to obtain the NPHS,     * 
!              * S_hat, of type nus for the residual T corresponding to*
!              * S and the NSPS, S_star_hat, of type nus for the       * 
!              * residual T_star corresponding to S_star. Then the     * 
!              * multiplications S*S_hat and S_star_hat*S_star yield   * 
!              * the NPHS and NSPS of types m + nus + (1,0,...,0) for A* 
!              * and A_star with index i+step. Exit from the inner loop* 
!              * takes place only if it  is determined that this point *
!              * m + nus + (1,0,...,0) is stable or if it is the last  * 
!              * point n.                                              *
!              *********************************************************


               call gen_next_vector(nu, nus, last_i, i+step, k)
               norm_nus = sum(nus)
               knorm_nus = k*norm_nus

!              Compute the residuals. 
               call build_T(A, S, k, m, prevnorm_nus, norm_nus, T)
               call build_T_star(A, S_star, k, m, prevnorm_nus, &
                                 norm_nus, T_star)
               prevnorm_nus = norm_nus

!              Determine S_hat of type nus for the residual T by solving 
!              striped Sylvester systems of equations.
               call build_S(T, k, nus, norm_nus, S_hat, singular)
               if (singular) flag = 2

!              Determine S_star_hat of type nus for the residual 
!              T_star by solving mosaic Sylvester systems of equations
               call build_S_star(T_star, k, nus, knorm_nus, &
                                        S_star_hat, singular)
               if (singular) flag = 2

!              The coordinates of the point with index i+step are 
!              stored temporarily in new_m.
               call gen_next_vector(n, new_m, last_i, i+step, k)
               norm_new_m = sum(new_m)

!              Obtain the NPHS and store in New_S
               call mult_S(S, m, S_hat, nus, New_S, new_m, k)

!              Obtain the NSPS and store in New_S_star
               call mult_S_star(S_star, m, S_star_hat, nus, &
                                        New_S_star, new_m, k)

!              Scale the NPHS and NSPS.
               call scale_S(New_S, new_m, k, gamma)
               call scale_S_star(New_S_star, new_m, k, gamma_star)

!              Compute the stability parameter.
               if (flag .eq. 2) then
                  kappa(i+step) = huge(0.0)
               else
                  kappa(i+step) = 0
                  do beta=0, k
                     kappa(i+step) = kappa(i+step)                    &
                                   + 1.0/(gamma(beta)*gamma_star(beta))
                  end do
                  if (kappa(i+step) .lt. tau) flag = 0
               endif

!           *  End of inner loop.                                      *
!           ************************************************************
            enddo


!           ************************************************************
!           * The point new_m with index i+step is either a stable     *
!           * point  or the last point n (or, both). So copy new_m to  *
!           * m, New_S to S and New_S_star to S_star.                  *
!           ************************************************************

            i = i+step
            m = new_m
            norm_m = norm_new_m

            do alpha=0, k
               do beta=0, k
                  if (beta .eq. 0) then
                     S_deg = m(alpha) + 1
                  else
                     S_deg = m(alpha)
                  endif
                  do l=0, S_deg
                     S(l,alpha,beta) = New_S(l,alpha,beta)
                  end do
               end do
            end do

            do alpha=0, k
               do beta=0, k
                  if (alpha .eq. 0) then
                     S_star_deg = norm_m -  m(beta)
                  else
                     S_star_deg = norm_m -  m(beta) + 1
                  endif
                  do l=0, S_star_deg
                     S_star(l,alpha,beta) = New_S_star(l,alpha,beta)
                  end do
              end do
           end do

!        * End of main loop.                                           *
!        ***************************************************************
         enddo

      endif
      deallocate (S_hat, New_S, T, S_star_hat, New_S_star, T_star)

      return
      end subroutine VECTOR_PADE


      subroutine build_T(A, S, k, m, prevnorm_nus, norm_nus, T)

!***********************************************************************
!                                                                      *
!     Given the vector of power series A and the Pade'-Hermite system  * 
!     S of type m, build_T returns the first norm_nus+1 terms of the   * 
!     residual power series T. It is assumed that prevnorm_nus+1 terms *
!     of T are already available and need not be computed.             *
!                                                                      *
!***********************************************************************
!                                                                      *
!     on entry                                                         *
!       A            real (0:sum(n), 0:k)                              *
!                    Vector of power series with A(0,0) nonzero.       *
!                                                                      *
!       k            integer                                           *
!                    There are k+1 power series in A.                  *
!                                                                      *
!       S            real (0:maxval(n)+1, 0:k, 0:k)                    *
!                    Normalized NPHS of type m.                        *
!                                                                      *
!       m            integer(0:k)                                      *
!                    The type specification of S.                      *
!                                                                      *
!       T            real (0:sum(n), 0:k)                              *
!                    The residual power series.                        *
!                    The first prevnorm_nus+1 terms are available on   *
!                    entry.                                            *
!                                                                      *
!       prevnorm_nus integer                                           *
!                    The first prevnorm_nus+1 terms of T are available *
!                    on entry.                                         *
!                                                                      *
!       norm_nus     integer                                           *
!                    This routine is required to compute the first     *
!                    norm_nus +1 terms of T (norm_nus > prevnorm_nus). *
!                                                                      *
!     on exit                                                          *
!       T            real (0:sum(n), 0:k)                              *
!                    The first norm_nus+1 terms of the residual power  *
!                    series.                                           *
!                                                                      *
!***********************************************************************

      implicit none

!     build_T subroutine parameters.
      integer,                       intent (in)    :: k, prevnorm_nus, &
                                                       norm_nus
      integer, dimension (0:),       intent (in)    :: m
      real,    dimension (0:,0:),    intent (in)    :: A
      real,    dimension (0:,0:,0:), intent (in)    :: S
      real,    dimension (0:,0:),    intent (inout) :: T

!     Local variables.
      integer alpha, beta, i, l, ll, start, finish, norm_m

      norm_m = sum(m)

!     The first prevnormus+1 terms of T are already available. For the
!     remaining terms of T, first initialize new residual terms to 0.
      do beta=0, k
         do i=prevnorm_nus+1, norm_nus
            T(i,beta) = 0
         end do
      end do

      do alpha=0, k
         do beta=0, k
            do l=prevnorm_nus+norm_m+2, norm_nus+norm_m+1
               if (beta .eq. 0) then
                  start = max(0, l-m(alpha)-1)
               else
                  start = max(0, l-m(alpha))
               endif
               finish = min(l, norm_nus+norm_m+1)
               do ll=start, finish
                  T(l-norm_m-1, beta) = T(l-norm_m-1, beta) &
                                      + A(ll,alpha)         &
                                      * S(l-ll,alpha,beta)
               end do
            end do
         end do
      end do

      return
      end subroutine build_T


      subroutine build_T_star(A, S_star, k, m, prevnorm_nus, norm_nus, &
                              T_star)
     
!***********************************************************************
!                                                                      *
!     Given the vector of power series A and the simultaneous Pade'    *
!     system S_star of type m, build_T_star returns the first          *
!     norm_nus+1 terms of the residual power series T_star. It is      *
!     assumed that prevnorm_nus+1 terms of T_star are already          *
!     available and need not be computed.                              *
!***********************************************************************
!                                                                      *
!     on entry                                                         *
!       A            real (0:sum(n), 0:k)                              *
!                    Vector of power series with A(0,0) nonzero.       *
!                                                                      *
!       k            integer                                           *
!                    There are k+1 power series in A.                  *
!                                                                      *
!       S_star       real (0:maxval(n)+1, 0:k, 0:k)                    *
!                    Normalized NSPS of type m.                        *
!                                                                      *
!       m            integer (0:k)                                     *
!                    The type specification of S_star.                 *
!                                                                      *
!       T_star       real (0:sum(n), 0:k, k)                           *
!                    The residual power series.                        *
!                    The first prevnorm_nus+1 terms are available on   *
!                    entry.                                            *
!                                                                      *
!       prevnorm_nus integer                                           *
!                    The first prevnorm_nus+1 terms of T_star are      *
!                    available on entry.                               *
!                                                                      *
!       norm_nus     integer                                           *
!                    This routine is required to compute the first     *
!                    norm_nus+1 terms of T_star                        *
!                    (norm_nus > prevnorm_nus).                        *
!                                                                      *
!     on exit                                                          *
!       T_star       real (0:sum(n), 0:k, k)                           *
!                    The first norm_nus+1 terms of the residual power  *
!                    series.                                           *
!                                                                      *
!***********************************************************************

      implicit none

!     build_T_star subroutine parameters.
      integer,                       intent (in)    :: k, prevnorm_nus, &
                                                       norm_nus
      integer, dimension (0:),       intent (in)    :: m
      real,    dimension (0:,0:),    intent (in)    :: A
      real,    dimension (0:,0:,0:), intent (in)    :: S_star
      real,    dimension (0:,0:,1:), intent (inout) :: T_star

!     local variables.
      integer alpha, beta, i, l, finish, norm_m

      norm_m = sum(m)
 
!     Multiply the vector A by the NSPS matrix S_star.
!     The first prevnormus+1 terms of T_star are already available, so
!     compute only the remaining terms.
      do alpha=0, k
         do beta=1, k
            do i=prevnorm_nus+1, norm_nus
               T_star(i,alpha,beta) = 0.0
               if (alpha .eq. 0) then
                  finish = min(i+norm_m+1, norm_m-m(0))
               else
                  finish = min(i+norm_m+1, norm_m-m(0) + 1)
               endif
               do l=0, finish
                  T_star(i, alpha, beta) = T_star(i, alpha, beta) &
                        - S_star(l,alpha,0) * A(norm_m+i-l+1,beta)
               end do
               if (alpha .eq. 0) then
                  finish = min(i+norm_m+1, norm_m-m(beta))
               else
                  finish = min(i+norm_m+1, norm_m-m(beta) + 1)
               endif
               do l=0, finish
                  T_star(i, alpha, beta) = T_star(i, alpha, beta) &
                        + S_star(l,alpha,beta) * A(norm_m+i-l+1,0)
               end do
            end do
         end do
      end do

      return
      end subroutine build_T_star


      subroutine build_S(A, k, n, norm_n, S, singular)

!***********************************************************************
!                                                                      *
!     For the vector of integers n = [n(0),...,n(k)], build_S computes * 
!     the Pade'-Hermite system S of type n for A by solving directly   * 
!     the striped Sylvester systems of linear equations associated     *
!     with A.                                                          *
!                                                                      *
!***********************************************************************
!                                                                      *
!     on entry                                                         *
!       A            real (0:sum(n), 0:k)                              *
!                    A vector of power series with A(0,0) nonzero.     *
!                                                                      *
!       k            integer                                           *
!                    There are k+1 power series in A.                  *
!                                                                      *
!       n            integer (0:k)                                     *
!                    The type specification of the NPHS.               *
!                                                                      *
!       norm_n       integer                                           *
!                    norm_n = n(0) + ... + n(k).                       *
!                                                                      *
!     on exit                                                          *
!       S            real (0:maxval(n)+1, 0:k, 0:k)                    *
!                    Normalized Pade'-Hermite system of type n.        *
!                                                                      *
!       singular     logical                                           *
!                    "true" if the associated striped Sylvester matrix *
!                    is singular and "false" otherwise.                *
!                                                                      *
!     functions and subroutines                                        *
!       sgefa             Linpack routine to triangulate a matrix.     *
!       sgesi             Linpack routine to solve linear equations.   *
!                                                                      *
!***********************************************************************

      implicit none

!     build_S subroutine parameters.
      integer,                       intent (in)  :: k, norm_n
      integer, dimension (0:),       intent (in)  :: n
      real,    dimension (0:,0:),    intent (in)  :: A
      real,    dimension (0:,0:,0:), intent (out) :: S
      logical,                       intent (out) :: singular

!     Local variables.
      integer alpha, beta, i, l, blockoffset, info

!     work areas - required by the subroutines sgefa and sgesi.
!     M is the striped Sylvester matrix associated with A.
      integer ipvt(norm_n)
      real M(norm_n, norm_n), B(norm_n)

!     Initialize S(z).
      do alpha=0, k
         S(0,alpha,0) = 0.0
         S(1,alpha,0) = 0.0
      end do
      do alpha=1, k
         do beta=1, k
            if (alpha .eq. beta) then
               S(0,alpha,beta) = 1.0
            else
               S(0,alpha,beta) = 0
            endif
         end do
      end do
      do beta=1, k
         S(0,0,beta) = -A(0,beta) / A(0,0)
      end do

      if (norm_n .eq. 0) then

!        Special case.
!        S is a diagonal matrix with modified first row.
         S(1,0,0) = 1 / A(0,0)
         singular = .false.

      else

!        Build the striped Sylvester matrix M.
         blockoffset = 0
         do beta=0, k
            do i=0, norm_n-1
               do l=0, n(beta)-1
                  if (l .le. i) then
                     M(i+1, blockoffset+l+1) = A(i-l, beta)
                  else
                     M(i+1, blockoffset+l+1) = 0
                  endif
               end do
             end do
             blockoffset = blockoffset + n(beta)
         end do

!        Reduce the system to triangular form.
         call sgefa(M, norm_n, norm_n, ipvt, info)

!        Compute the first column of S.
         do i=1, norm_n-1
            B(i) = 0
         end do
         B(norm_n) = 1
         call sgesl(M, norm_n, norm_n, ipvt, B, 0)
         blockoffset = 0
         do alpha=0, k
            do l=1, n(alpha)
               S(l+1,alpha,0) = B(blockoffset+l)
            end do
            blockoffset = blockoffset + n(alpha)
         end do

!        Compute the last k columns of S.
         if (info .gt. 0) then
            singular = .true.
            do beta=1, k
               do alpha=0, k
                  do l=0, n(alpha)
                     S(l,alpha,beta) = 0
                  end do
               end do
            end do
         else
            singular = .false.
            do beta=1, k
               do i=1, norm_n
                  B(i) = -A(i,beta)+(A(i,0) * A(0,beta)/A(0,0))
               end do
               call sgesl(M, norm_n, norm_n, ipvt, B, 0)
               blockoffset = 0
               do alpha=0, k
                  do l=1, n(alpha)
                     S(l,alpha,beta) = B(blockoffset+l)
                  end do
                  blockoffset = blockoffset + n(alpha)
               end do
            end do
         endif
      endif

      return
      end subroutine build_S


      subroutine build_S_star(A_star, k, n, knorm_n, S_star, singular)

!***********************************************************************
!                                                                      *
!     For the vector of integers n = [n(0),...,n(k)], build_S_star     *
!     computes the simultaneous Pade' system S_star of type n for      *
!     A_star by solving directly the mosaic Sylvester systems of       *
!     linear equations associated with A_star.                         *
!                                                                      *
!***********************************************************************
!                                                                      *
!    on entry                                                          *
!       A_star       real(0:sum(n), 0:k,1:k)                           *
!                    k+1 x k matrix of power series.                   *
!                                                                      *
!       k            integer                                           *
!                    Dimension parameter.                              *
!                                                                      *
!       n            integer(0:k)                                      *
!                    The type specification of the NSPS.               *
!                                                                      *
!       knorm_n      integer                                           *
!                    knorm_n = k * (n(0) + ... + n(k)).                *
!                                                                      *
!    on exit                                                           *
!       S_star       real(0:sum(n), 0:k, 0:k)                          *
!                    Normalized Pade'-Hermite system of type n.        *
!                                                                      *
!       singular     logical                                           *
!                    "true" if the associated mosaic Sylvester matrix  *
!                    is singular and "false" otherwise.                *
!                                                                      *
!    functions and subroutines                                         *
!       sgefa             Linpack routine to triangulate a matrix.     *
!       sgesi             Linpack routine to solve linear equations.   *
!                                                                      *
!***********************************************************************

      implicit none

!     build_S_star subroutine parameters.
      integer,                       intent (in)  :: k, knorm_n
      integer, dimension (0:),       intent (in)  :: n
      real,    dimension (0:,0:,1:), intent (in)  :: A_star
      real,    dimension (0:,0:,0:), intent (out) :: S_star
      logical,                       intent (out) :: singular

!     Local variables.
      integer alpha, beta, l, i, j, info, norm_n,                   &
               blockoffset, cblockoffset, rblockoffset

!     work areas - required by the subroutines sgefa and sgesi.
!     M_star is the mosaic Sylvester matrix associated with A_star.
      integer ipvt(knorm_n)
      real M_star(knorm_n,knorm_n), B_star(knorm_n)

      norm_n = sum(n)
      
!     Initialize S_star to zero.
      do alpha=1, k
         do beta=0, k
            do l=0, 1
               S_star(l,alpha,beta) = 0
            end do
         end do
      end do

!     Set the constant terms of the first row of S_star.
      S_star(0,0,0) = 1.0
      do beta=1, k
         S_star(0,0,beta) = - A_star(0,0,beta) / A_star(0,beta,beta)
      end do
         
      if (norm_n .eq. 0) then

!        Special case.
!        S_star is a diagonal matrix with modified first row.
         do alpha=1, k
            S_star(1, alpha, alpha) = 1 / A_star(0,alpha,alpha)
         end do

      else

!        Build the mosaic matrix M_star of order knorm_n.
         rblockoffset = 0
         do beta=1, k
            cblockoffset = 0
            do alpha=0, k
               do i=1, norm_n
                  do j=1, norm_n - n(alpha)
                     If (i .lt. j) then
                        M_star(rblockoffset+i,cblockoffset+j) = 0
                     else
                        M_star(rblockoffset+i,cblockoffset+j)         &
                                             = A_star(i-j, alpha, beta)
                     endif
                  end do
               end do
               cblockoffset = cblockoffset + norm_n - n(alpha)
            end do
            rblockoffset = rblockoffset + norm_n
         end do

!        Reduce M_star into triangular form.
         call sgefa(M_star,knorm_n,knorm_n,ipvt,info)

!        Compute the first row of S_star.
         blockoffset = 0
         do beta=1, k
            do i=1, norm_n
               B_star(blockoffset+i) = -A_star(i,0,beta)
               do j=1, k
                  B_star(blockoffset+i) = B_star(blockoffset+i) &
                             - A_star(i,j,beta)*S_star(0,0,j) 
               end do
            end do
            blockoffset = blockoffset + norm_n
         end do
         call sgesl(M_star, knorm_n, knorm_n, ipvt, B_star, 0)
         blockoffset = 0
         do beta=0, k
            do l=1, norm_n-n(beta)
               S_star(l,0,beta) = B_star(blockoffset + l)
            end do
            blockoffset = blockoffset + norm_n - n(beta)
         end do

         if (info .gt. 0) then
!           M_star is singular
            singular = .true.
            do alpha=1, k
               do beta=0, k
                  do l=0, norm_n-n(beta)+1
                     S_star(l,alpha,beta) = 0.0
                  end do
               end do
            end do
         else
            singular = .false.
            do alpha=1, k
               do i=1, knorm_n
                  B_star(i) = 0
               end do
               B_star(alpha*norm_n) = 1
               call sgesl(M_star, knorm_n, knorm_n, ipvt, B_star, 0)
               blockoffset = 0
               do beta=0, k
                  do l=1, norm_n-n(beta)
                     S_star(l+1,alpha,beta)                      &
                                         = B_star(blockoffset + l)
                  end do
                  blockoffset = blockoffset + norm_n - n(beta)
               end do
            end do

         endif
      endif

      return
      end subroutine build_S_star


      subroutine mult_S(S, m, S_hat, nus, New_S, new_m, k)

!***********************************************************************
!                                                                      *
!     mult_S multiplies the k+1 x k+1 matrix of polynomials S and      *
!     S_hat, where S is a NPHS of type m and S_hat is a NPHS of type   *
!     nus. The product New_S = S * S_hat is a NPHS of type new_m,      *
!     where new_m = m + nus + (1,0,...0).                              *
!                                                                      *
!     on entry                                                         *
!        S         real (0:maxval(n)+1,0:k,0:k)                        *
!                  NPHS of type m.                                     *
!                                                                      *
!        m         integer (0:k)                                       *
!                  The type specification of S.                        *
!                                                                      *
!        S_hat     real (0:maxval(n)+1,0:k,0:k)                        *
!                  NPHS of type nus.                                   *
!                                                                      *
!        nus       integer (0:k)                                       *
!                  The type specification of S_hat.                    *
!                                                                      *
!        new_m     integer (0:k)                                       *
!                  The type specification of New_S.                    *
!                  Must satisfy new_m = m + nus + (1,0,...0).          *
!                                                                      *
!        k         integer                                             *
!                  Dimension parameter.                                *
!                                                                      *
!     on exit                                                          *
!      New_S  real (0:maxval(n)+1,0:k,0:k)                             *
!                  The result of S * S_hat. This is the NPHS of type   * 
!                  new_m.                                              *
!                                                                      *
!***********************************************************************

      implicit none

!     mult_S subroutine parameters.
      integer,                       intent (in)  :: k
      integer, dimension (0:),       intent (in)  :: m, nus, new_m
      real,    dimension (0:,0:,0:), intent (in)  :: S, S_hat
      real,    dimension (0:,0:,0:), intent (out) :: New_S

!     Local variables.
      integer alpha, beta, i, j, l, S_deg, start, finish
 
      do alpha=0, k
         do beta=0, k
            if (beta .eq. 0) then
               S_deg = new_m(alpha) + 1
            else 
               S_deg = new_m(alpha)
            endif
            do l=0, S_deg 
               New_S(l,alpha,beta) = 0
               do j=0, k
                  if (beta .eq. 0) then 
                     start = max(0, l - nus(j) -1)
                  else
                     start = max(0, l - nus(j))
                  endif
                  if (j .eq. 0) then 
                     finish = min(l, m(alpha) + 1)
                  else
                     finish = min(l, m(alpha))
                  endif
                  do i=start, finish
                     New_S(l,alpha,beta) = New_S(l,alpha,beta)      &
                                   + S(i,alpha,j) * S_hat(l-i,j,beta)
                  end do
               end do
            end do
         end do
      end do

      return
      end subroutine mult_S


      subroutine mult_S_star(S_star, m, S_star_hat, nus, &
                                     New_S_star, new_m, k)

!***********************************************************************
!                                                                      *
!     mult_S_star multiplies the k+1 x k+1 matrix of polynomials       *
!     S_star and S_star_hat, where S_star is a NSPS of type m and      *
!     S_star_hat is a NSPS of type nus. The product                    *
!     New_S_star = S_star_hat * S_star is a NSPS of type new_m where   *
!     new_m = m + nus + (1,0,...0).                                    *
!                                                                      *
!     on entry                                                         *
!        S_star    real (0:sum(n),0:k,0:k)                             *
!                  NSPS of type m.                                     *
!                                                                      *
!        m         integer (0:k)                                       *
!                  The type specification of S_star.                   *
!                                                                      *
!      S_star_hat  real (0:sum(n),0:k,0:k)                             *
!                  NSPS of type nus.                                   *
!                                                                      *
!        nus       integer (0:k)                                       *
!                  The type specification of S_star_hat.               *
!                                                                      *
!        new_m     integer (0:k)                                       *
!                  The type specification of New_S_star.               *
!                  Must satisfy new_m = m + nus + (1,0,...0).          *
!                                                                      *
!        k         integer                                             *
!                  Dimension parameter.                                *
!                                                                      *
!     on exit                                                          *
!      New_S_star  real (0:sum(n),0:k,0:k)                             *
!                  The result of S_star_hat * S_star. This is the NSPS *
!                  of type m_new.                                      *
!                                                                      *
!***********************************************************************

      implicit none

!     mult_S_star subroutine parameters.
      integer,                       intent (in)  :: k
      integer, dimension (0:),       intent (in)  :: m, nus, new_m
      real,    dimension (0:,0:,0:), intent (in)  :: S_star, S_star_hat
      real,    dimension (0:,0:,0:), intent (out) :: New_S_star

!     Local variables.
      integer alpha, beta, i, j, l, S_star_deg, start, finish, &
              norm_m, norm_nus, norm_new_m
 
      norm_m = sum(m)
      norm_nus = sum(nus)
      norm_new_m = sum(new_m)

      do alpha=0, k
         do beta=0, k
            if (alpha .eq. 0) then
               S_star_deg = norm_new_m - new_m(beta) 
            else 
               S_star_deg = norm_new_m - new_m(beta) + 1
            endif
            do l=0, S_star_deg 
               New_S_star(l,alpha,beta) = 0.0
               do j=0, k
                  if (j .eq. 0) then 
                     start = max(0, l - norm_m + m(beta))
                  else
                     start = max(0, l - norm_m + m(beta) - 1)
                  endif
                  if (alpha .eq. 0) then 
                     finish = min(l, norm_nus - nus(j))
                  else
                     finish = min(l, norm_nus - nus(j) + 1)
                  endif
                  do i=start, finish
                     New_S_star(l,alpha,beta) = New_S_star(l,alpha,beta) &
                            + S_star_hat(i,alpha,j) * S_star(l-i,j,beta)
                  end do
               end do
            end do
         end do
      end do

      return
      end subroutine mult_S_star


      subroutine scale_S(S, m, k, gamma)

!***********************************************************************
!                                                                      *
!     Scale the NPHS so that the 1-norm of each column is equal to 1.  *
!                                                                      *
!***********************************************************************
!                                                                      *
!     on entry                                                         *
!       S            real (0:maxval(n)+1, 0:k, 0:k)                    *
!                    NPHS with r(0) = 1.                               *
!                                                                      *
!       k            integer                                           *
!                    S is a k+1 x k+1 matrix of polynomials.           *
!                                                                      *
!     on exit                                                          *
!       S            real (0:maxval(n)+1, 0:k, 0:k)                    *
!                    The scaled NPHS.                                  *
!                                                                      *
!       gamma        real(0:k)                                         *   
!                    The vector of scaling factors. The normalized NPHS*
!                    can be obtained from the scaled one by dividing   *
!                    the beta'th column by gamma(beta).                *
!                                                                      *
!***********************************************************************

      implicit none

!     scale_S subroutine parameters.
      integer,                       intent (in)    :: k
      integer, dimension (0:),       intent (in)    :: m
      real,    dimension (0:,0:,0:), intent (inout) :: S
      real,    dimension (0:),       intent (out)   :: gamma

!     Local variables.
      integer alpha, beta, S_deg, l

!     Initially, determine gamma(beta) to be the 1-norm of the beta'th 
!     column of S.
      do beta=0, k
         gamma(beta) = 0
         do alpha=0, k
            if (beta .eq. 0) then
               S_deg = m(alpha) + 1
            else
               S_deg = m(alpha) 
            endif
            do l=0, S_deg
               gamma(beta) = gamma(beta) + abs(S(l,alpha,beta))
            end do
          end do
      end do


!     Scale each column of S.
      do beta=0, k
         if (gamma(beta) .ne. 0.0) then
            do alpha=0, k
               if (beta .eq. 0) then
                  S_deg = m(alpha) + 1
               else
                  S_deg = m(alpha) 
               endif
               do l=0, S_deg
                  S(l,alpha,beta) = S(l,alpha,beta) / gamma(beta)
               end do
            end do
         endif
      end do

!     It is assumed on entry that r(0)=1, so gamma(0) is simply the 
!     inverse of the norm of the first column of S prior to scaling.
!     The remaining scaling factors are determined from the diagonal
!     of V(0).
      gamma(0) = 1/gamma(0)
      do beta=1, k
         gamma(beta) = abs(S(0,beta,beta))
      end do

      return
      end subroutine scale_S


      subroutine scale_S_star(S_star, m, k, gamma_star)

!***********************************************************************
!                                                                      *
!     Scale the NSPS so that the norm of each row is equal to 1.       *
!                                                                      *
!***********************************************************************
!                                                                      *
!     on entry                                                         *
!       S_star       real (0:sum(n), 0:k, 0:k)                         *
!                    NSPS with R_star(0) = I.                          *
!                                                                      *
!       k            integer                                           *
!                    S_star is a k+1 x k+1 matrix of polynomials.      *
!                                                                      *
!     on exit                                                          *
!       S_star       real (0:sum(n), 0:k, 0:k)                         *
!                    The scaled NPHS.                                  *
!                                                                      *
!       gamma_star   real(0:k)                                         *
!                    The vector of scaling factors. The normalized NSPS*
!                    can be obtained from the scaled one by dividing   *
!                    the alpha'th row by gamma(alpha).                 *
!                                                                      *
!***********************************************************************

      implicit none

!     scale_S_star subroutine parameters.
      integer,                       intent (in)    :: k
      integer, dimension (0:),       intent (in)    :: m
      real,    dimension (0:,0:,0:), intent (inout) :: S_star
      real,    dimension (0:),       intent (out)   :: gamma_star

!     Local variables.
      integer norm_m, alpha, beta, S_star_deg, l

      norm_m = sum(m)

!     Initially, determine gamma(alpha) to be the 1-norm of the alpha'th
!     row of S_star.
      do alpha=0, k
         gamma_star(alpha) = 0
         do beta=0, k
            if (alpha .eq. 0) then
               S_star_deg = norm_m - m(beta)
            else
               S_star_deg = norm_m - m(beta)+1
            endif
            do l=0, S_star_deg
               gamma_star(alpha) = gamma_star(alpha)       &
                                 + abs(S_star(l,alpha,beta))
            end do
         end do
      end do

!     Scale each row of S_star.
      do alpha=0, k
         if (gamma_star(alpha) .ne. 0.0) then
            do beta=0, k
               if (alpha .eq. 0) then
                  S_star_deg = norm_m - m(beta)
               else
                  S_star_deg = norm_m - m(beta)+1
               endif
               do l=0, S_star_deg
                  S_star(l,alpha,beta) = S_star(l,alpha,beta) &
                                       / gamma_star(alpha)
               end do
            end do
         endif
      end do
      
!     It is assumed on entry that R_star(0) = I, so that, 
!     for alpha=1,...,k, gamma_star(alpha) is simply the inverse 
!     of the norm of the first row of S_star prior to scaling.
      gamma_star(0) = abs(S_star(0,0,0))
      do alpha=1, k
         if (gamma_star(alpha) .ne. 0.0) then
            gamma_star(alpha) = 1.0/gamma_star(alpha) 
          endif
      end do

      return
      end subroutine scale_S_star


      subroutine gen_next_vector(n, m, last_i, i, k)

!***********************************************************************
!                                                                      *
!     gen_next_point generates the next point m beyond the point with  * 
!     index i along the diagonal of the Pade' table passing through    *
!     the point n.                                                     *
!                                                                      *
!     on entry                                                         *
!        n      integer (0:k)                                          *
!               The final point in the Pade' table.                    *
!                                                                      *
!        last_i integer                                                *
!               The index of the point n.                              *
!                                                                      *
!        i      integer                                                *
!               The current index.                                     *
!                                                                      *
!        k      integer                                                *
!               The index of the last entry in n.                      *
!                                                                      *
!     on exit                                                          *
!        m      integer (0:k)                                          *
!               The next point along the diagonal through n.           *
!                                                                      *
!***********************************************************************

      implicit none
      integer k
      integer  n(0:k), m(0:k), i, last_i

!     Local variable.
      integer beta

      do beta=0, k
         m(beta) = max(0, n(beta) - last_i + i)
      end do
      return
      end subroutine gen_next_vector



!***********************************************************************
!     This is a test driver for the subroutines VECTOR_PADE,           *
!     build_delta_T and build_delta_T_star. The driver requires input  *
!     data.                                                            *
!***********************************************************************

      module working_area_main
!        Assummed shape dummy arguments for build_delta_T.
         real, dimension (:,:),   allocatable :: delta_T
!        Assummed shape dummy arguments for build_delta_T_star.
         real, dimension (:,:,:), allocatable :: delta_T_star
!        Assumed shape dummy arguments for VECTOR_PADE.
         integer, dimension (:),  allocatable :: n
         real, dimension (:),     allocatable :: gamma, gamma_star, kappa
         real, dimension (:,:),   allocatable :: A
         real, dimension (:,:,:), allocatable :: S, S_star
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

      program main

      use working_area_main
      implicit none

      interface
         subroutine VECTOR_PADE(k, n, A, tau,                          &
                   S, gamma, S_star, gamma_star, kappa, num_steps, flag)
            integer,                    intent (in)    :: k
            integer, dimension (:),     intent (in)    :: n
            real,    dimension (:,:),   intent (in)    :: A
            real,                       intent (in)    :: tau
            real,    dimension (:,:,:), intent (out)   :: S, S_star
            real,    dimension (:),     intent (out)   :: gamma,   &
                                                          gamma_star, kappa
            integer,                    intent (out)   :: num_steps
            integer,                    intent (inout) :: flag
         end subroutine VECTOR_PADE
         subroutine build_delta_T(k, n, A, S, delta_T, delta_T_norm)
            integer,                    intent (in)  :: k
            integer, dimension (:),     intent (in)  :: n
            real,    dimension (:,:),   intent (in)  :: A
            real,    dimension (:,:,:), intent (in)  :: S
            real,                       intent (out) :: delta_T_norm
            real,    dimension (:,:),   intent (out) :: delta_T
         end subroutine build_delta_T
         subroutine build_delta_T_star(k, n, A, S_star,               &
                                       delta_T_star, delta_T_star_norm)
            integer,                    intent (in)  :: k
            integer, dimension (:),     intent (in)  :: n
            real,    dimension (:,:),   intent (in)  :: A
            real,    dimension (:,:,:), intent (in)  :: S_star
            real,                       intent (out) :: delta_T_star_norm
            real,    dimension (:,:,:), intent (out) :: delta_T_star
         end subroutine build_delta_T_star
      end interface

!     Local variables.
      integer k, flag, num_steps, S_deg, S_star_deg, alpha, beta, l
      integer norm_n
      real tau, delta_T_norm, delta_T_star_norm 

!     Input the tolerance tau and the dimension of the problem, k.
      read *, tau
      read *, k
      allocate (n(0:k))

!     Input n.
      do beta=0, k
         read *, n(beta)
      end do
      norm_n = sum(n)

      allocate (A(0:norm_n, 0:k),                                     &
                S(0:maxval(n)+1,0:k,0:k), gamma(0:k),                 &
                S_star(0:norm_n,0:k,0:k), gamma_star(0:k),            &
                kappa(0:norm_n),                                      &
                delta_T(0:norm_n, 0:k), delta_T_star(0:norm_n, 0:k, k))

!     Input the power series A.
      do beta=0, k
         do l=0, norm_n
            read *, A(l,beta)
         end do
      end do

!     Compute the Pade'-Hermite system S of type n and the
!     simultaneous Pade' system S_star of type n.
      call VECTOR_PADE(k, n, A, tau,                                   &
                   S, gamma, S_star, gamma_star, kappa, num_steps, flag)

!     Compute the errors in the order conditions.
      call build_delta_T(k, n, A, S, delta_T, delta_T_norm)
      call build_delta_T_star(k, n, A, S_star,               &
                              delta_T_star, delta_T_star_norm)

      print 190, flag
      if (flag .lt. 3) then
         print 191, tau
         print 192
         do l=0, num_steps
            print 193, l, kappa(l)
         end do
         print 194, delta_T_norm
         print 195, delta_T_star_norm
         do alpha=0, k
            do beta=0, k
               if (beta .eq. 0) then
                  S_deg = n(alpha) + 1
               else
                  S_deg = n(alpha)
               endif
               print 196, alpha, beta
               do l=0, S_deg
                  print 197, l, S(l,alpha,beta)
              end do
            end do
         end do

         do alpha=0, k
            do beta=0, k
               if (alpha.eq. 0) then
                  S_star_deg = norm_n - n(beta)
               else
                  S_star_deg = norm_n - n(beta) + 1
               endif
               print 198, alpha, beta
               do l=0, S_star_deg
                  print 199, l, S_star(l,alpha,beta)
               end do
            end do
         end do

      endif
      deallocate (n, A, S, gamma, S_star, gamma_star, kappa, delta_T, &
                                                          delta_T_star)
      stop

 190  format(' flag = ',I1)
 191  format(' Stability tolerance = ',D9.2)
 192  format(' Step   kappa')
 193  format(I4, D11.2)
 194  format(' Norm of error in NPHS order condition = ',D9.2)
 195  format(' Norm of error in NSPS order condition = ',D9.2)
 196  format(/'        S(',I3,',',I3,')')
 197  format('   z**',I3,':  ',D15.8)
 198  format(/'   S_star(',I3,',',I3,')')
 199  format('   z**',I3,':  ',D15.8)

      end program main

  
      subroutine build_delta_T(k, n, A, S, delta_T, delta_T_norm)

!***********************************************************************
!                                                                      *
!     Computes the error in the residual given the vector of           *
!     power series A and the NPHS S of type n. That is,                *
!                                                                      *
!               delta_T = A * S (mod z**(||n||+1)).                    *
!                                                                      *
!                                                                      *
!     On entry:                                                        *
!        A      real (0:sum(n),0:k)                                    *
!               Vector of power series.                                *
!                                                                      *
!        k      integer                                                *
!               There are k+1 power series in A.                       *
!                                                                      *
!        S      real (0:maxval(n)+1, 0:k, 0:k)                         *
!               Pade Hermite system for A of type n.                   *
!                                                                      *
!        n      integer (0:k)                                          *
!               The type specification of S.                           *
!                                                                      *
!     on exit                                                          *
!       delta_T real (0:sum(n), 0:k)                                   *
!               Error in the residual.                                 *
!                                                                      *
!       delta_T_norm  real                                             *
!               1-norm of the row vector delta_T.                      *
!                                                                      *
!***********************************************************************

      implicit none

!     build_delta_T subroutine parameters.
      integer,                       intent (in)  :: k
      integer, dimension (0:),       intent (in)  :: n
      real,    dimension (0:,0:),    intent (in)  :: A
      real,    dimension (0:,0:,0:), intent (in)  :: S
      real,                          intent (out) :: delta_T_norm
      real,    dimension (0:,0:),    intent (out) :: delta_T

!     Local variables.
      integer alpha, beta, l, i, start, finish, S_deg, norm_n
      real colsum 

      norm_n = sum(n)
      delta_T = 0

!     Multiply the vector A by the matrix S.
      do alpha=0, k
         do beta=0, k
            do l=0, norm_n
               if (beta .eq. 0) then
                  S_deg = n(alpha) + 1
               else
                  S_deg = n(alpha)
               endif
               start = max(0, l-S_deg)
               finish = min(l, norm_n)
               do i=start, finish
                  delta_T(l, beta) = delta_T(l, beta)              & 
                                   + A(i,alpha)  * S(l-i,alpha,beta)
               end do
            end do
         end do
      end do

!     Compute the 1-norm of delta_T.
      delta_T_norm = 0.0
      do beta=0, k
         colsum = 0.0
         do l=0, norm_n
            colsum = colsum + abs(delta_T(l,beta))
         end do
            delta_T_norm = max(colsum, delta_T_norm)
      end do

      return
      end subroutine build_delta_T


      subroutine build_delta_T_star(k, n, A, S_star,               &
                                    delta_T_star, delta_T_star_norm)

!***********************************************************************
!                                                                      *
!     Computes the error in the residual given the vector of           *
!     power series A and the NSPS S_star of type n. That is,           *
!                                                                      *
!            delta_T_star = S_star * A_star (mod z**(||n||+1)),        *
!                                                                      *
!     where                                                            *
!                                                                      *
!                            |-A(1) -A(2) ... -A(k)|                   * 
!                            | A(0)                |                   *
!                  A_star *  |       A(0)          |.                  * 
!                            |            ...      |                   *    
!                            |                 A(0)|                   *
!                                                                      *
!                                                                      *
!                                                                      *
!     On entry:                                                        *
!        A      real (0:sum(n),0:k)                                    *
!               Vector of power series.                                *
!                                                                      *
!        k      integer                                                *
!               There are k+1 power series in A.                       *
!                                                                      *
!        S_star real (0:sum(n), 0:k, 0:k)                              *
!               Simultaneous Pade' system corresponding to A of type n.*
!                                                                      *
!        n      integer (0:k)                                          *
!               The type specification of S_star.                      *
!                                                                      *
!     on exit                                                          *
!       delta_T_star     real (0:sum(n), 0:k, k)                       *
!                        Error in the residual.                        *
!                                                                      *
!       delta_T_star_norm  real                                        *
!                          1-norm of delta_T_star.                     *
!                                                                      *
!***********************************************************************

      implicit none

!     build_delta_T_star subroutine parameters.
      integer,                       intent (in)  :: k
      integer, dimension (0:),       intent (in)  :: n
      real,    dimension (0:,0:),    intent (in)  :: A
      real,    dimension (0:,0:,0:), intent (in)  :: S_star
      real,                          intent (out) :: delta_T_star_norm
      real,    dimension (0:,0:,1:), intent (out) :: delta_T_star
!     Local variables.

!     Local variables.
      integer alpha, beta, l, i, finish, S_star_deg, norm_n
      real series_norm, colsum

      norm_n = sum(n)
      delta_T_star = 0

!     Multiply the matrix S_star by the matrix A_star.
      do alpha=0, k
         do beta=1, k
            do l=0, norm_n
               if (alpha .eq. 0) then
                  S_star_deg = norm_n - n(beta)
               else
                  S_star_deg = norm_n - n(beta) + 1
               endif
               finish = min(l, S_star_deg)
               do i=0, finish
                  delta_T_star(l, alpha, beta)                   &
                                = delta_T_star(l, alpha, beta)   &
                                + S_star(i,alpha,beta) * A(l-i,0)
               end do
               if (alpha .eq. 0) then
                  S_star_deg = norm_n - n(0)
               else
                  S_star_deg = norm_n - n(0) + 1
               endif
               finish = min(l, S_star_deg)
               do i=0, finish
               delta_T_star(l, alpha, beta)                          &
                                    = delta_T_star(l, alpha, beta)   &
                                    - S_star(i,alpha,0) * A(l-i,beta)
               end do
            end do
         end do
      end do

!     Compute the 1_norm of delta_T_star.

      delta_T_star_norm = 0.0
      do beta=1, k
         colsum = 0
         do alpha=0, k
            series_norm = 0
            do l=0, norm_n
               series_norm = series_norm                   & 
                           + abs(delta_T_star(l,alpha,beta))
            end do
            colsum = colsum + series_norm 
         end do
         delta_T_star_norm = max(delta_T_star_norm, colsum)
      end do

      return
      end subroutine build_delta_T_star



program testeigentest


!  This program tests the package Eigentest.  It runs 64 test cases
!  probing various aspects of the package, as described below.  The
!  numbers in the output should be within two or so orders of
!  magnitude of the rounding unit.

   use eigentest

   implicit none
 
   ! These values may be adjusted at compile time to change the
   ! tests, which involve computing (A - shift*I) op B,
   ! where op is *, '*,  \, or '\ (' denoting transpose).

   integer,  parameter :: n = 10        ! Order of A.
   integer,  parameter :: ncols = 3     ! Number of columns in B.
   real(wp), parameter :: shift = -1.0  ! A shift.
   integer,  parameter :: spinrand = 5  ! Number of times to call the random
                                        ! number generator initially.

   type(Eigenmat) :: A
   real(wp):: B(n,n), C(n,n), D(n,n)
   complex(wp):: BB(n,n), CC(n,n), DD(n,n), ev(n)
   real(wp) :: temp, cond, norm, ev1(n)
   integer :: cn, i, j, ztype, atype, nb(4) = (/1,2,3,2/), i1, i2
   logical :: idy
   

   ! Spin the random number generator.

   do i=1,spinrand
      call Random_Number(temp)
   end do

   ! Loop on test case number.
   
   do cn = 0, 63
      print *, ""
      print *, "test case ", cn
   
      ! Set up the test case.

      ! Choose the type of the outer hsvd Y.  It can be an identity
      !  (idy == .true.) or a random hsvd (idy == .false.).

      idy = (mod(cn,2) == 0)

      ! Choose the typep of the inner hsvd Z.
      !
      !   If ztype == 0, Z is an identity.
      !   If ztype == 1, Z has two blocks: (1:6)(7:10).
      !   The first block is identity.
      !   If ztype == 2, Z has three blocks: (1:2)(3:8)(9:10),
      !     and the second block is an indentity.
      !   If ztype == 3, Z has two blocks: (1:4)(5::10),
      !     and the second block is identity.

      ztype = mod(cn/2,4)

      ! Chose the types of eigenvalues.
      !
      !   If atype == 0, all eigenvalues are real.
      !   If atype == 1, all eigenvalues are complex.
      !   If atype == 2, eigenvalues 1,2,5,6,9,10 are complex.
      !   If atype == 3, eigenvalues 3,4,7,8 are complex.
      !   If atype == 4, all eigenvalues are in a Jordan block.
      !   If atype == 5, eigenvalues (3,4),(7,8,9) are in Jordan blocks.
      !   If atype == 6, eigenvalues (1-3),(8-10) are in a Jordan block.
      !   If atype == 7, eigenvalues (2-4) are in a Jordan block.
      !                  eigenvalues 6,7,8,9 are complex.

      atype = mod(cn/8,8)

      ! Initialize A.

      call EigenmatFree(A)
      call EigenmatAlloc(A, n, nb(ztype+1), idy, ztype == 0)
      
      ! Generate eigenvalues.

      call Random_Number(A%eig(1:n))
      A%eig = A%eig + 0.5

      select case(atype)
      case(0)
         A%type = 1

      case(1)
         A%type(1:n:2) = 2
         A%type(2:n:2) = 3

      case(2)
         A%type((/3,4,7,8/)) = 1
         A%type((/1,5,9/))   = 2
         A%type((/2,6,10/))  = 3

      case(3)
         A%type((/1,2,5,6,9,10/)) = 1
         A%type((/3,7/))  = 2
         A%type((/4,8/))  = 3

      case(4)
         A%type(1)   = -n
         A%type(2:n) = -1

      case (5)
         A%type(3)   = -2
         A%type(7)   = -3
         A%type((/4,8,9/)) = -1
         A%type((/1,2,5,6,10/)) = 1

      case (6)
         A%type((/1,8/)) = -3
         A%type(2:3)    = -1
         A%type(9:10)   = -1
         A%type(4:7)    =  1
               
      case (7)
         A%type(2)      = -3
         A%type(3:4)    = -1
         A%type((/1,5,10/)) = 1
         A%type((/6,8/)) = 2
         A%type((/7,9/)) = 3

      end select      

      print *, "A%type:", A%type

      ! Generate Y.

      if (.not. idy) then

         call Random_Number(A%Y%u)
         call Random_Number(A%Y%v)
         call Random_Number(A%Y%sig)

         A%Y%u = A%Y%u - 0.5
         A%Y%v = A%Y%v - 0.5
         A%Y%sig = A%Y%sig + 1.0

         call hscal(A%Y%u)
         call hscal(A%Y%v)

         A%Y%bs(1) = 1
      end if

      print *, "Y%bs:", A%Y%bs(1), A%Y%bs(2);

      ! Generate Z.

      if (ztype > 0) then

         call Random_Number(A%Z%u)
         call Random_Number(A%Z%v)
         call Random_Number(A%Z%sig)
         
         A%Z%sig = A%Z%sig + 0.5
         A%Z%u = A%Z%u - 0.5
         A%Z%v = A%Z%v - 0.5
         
         select case(ztype)
            
         case(1)
            ! If ztype == 1, Z has two blocks: (1:6)(7:10).
            ! The first block is identity.

            A%Z%bs(1) = 1
            A%Z%bs(2) = -7
            
         case(2)
            ! If ztype == 2, Z has three blocks: (1:2)(3:8)(9:10).
            ! The second block is indentity.         

            A%Z%bs(2) = 3
            A%Z%bs(3) = -9
            
         case(3)
            ! If ztype == 3, Z has two blocks: (1:4)(5::10).
            ! The second block is identity.

            A%Z%bs(2) = 5
            A%Z%bs(3) = -n-1
            
         end select

         do j = 1, nb(ztype+1)
            i1 =  abs(A%Z%bs(j))
            i2 =  A%Z%bs(j+1)-1
            if (i2 > 0) then
               call hscal(A%Z%u(i1:i2))
               call hscal(A%Z%v(i1:i2))
            end if
         end do
      end if

      print *, "Z%bs:", A%Z%bs

      ! Generate B for A, Y, and Z to operate on.

      call Random_Number(B(1:n,1:ncols))
      B(1:n,1:ncols) = B(1:n,1:ncols) - 0.5;

      ! Test Z.

      C(1:n,1:ncols) = B(1:n,1:ncols);
      call HsvdProd(A%Z, ncols, C, n, 'ab')
      call HsvdProd(A%Z, ncols, C, n, 'aib')
      print *, "|inv(Z)*Z*x-x|  =", sum(abs((B(1:n,1:ncols)-C(1:n,1:ncols))));
      
      C(1:n,1:ncols) = B(1:n,1:ncols);
      call HsvdProd(A%Z, ncols, C, n, 'aib')
      call HsvdProd(A%Z, ncols, C, n, 'ab')
      print *, "|Z*inv(Z)*x-x|  =", sum(abs((B(1:n,1:ncols)-C(1:n,1:ncols))));
      
      C(1:n,1:ncols) = B(1:n,1:ncols);
      call HsvdProd(A%Z, ncols, C, n, 'atb')
      call HsvdProd(A%Z, ncols, C, n, 'aitb')
      print *, "|inv(Z')*Z'*x-x|=",sum(abs((B(1:n,1:ncols)-C(1:n,1:ncols))));
      
      C(1:n,1:ncols) = B(1:n,1:ncols);
      call HsvdProd(A%Z, ncols, C, n, 'aitb')
      call HsvdProd(A%Z, ncols, C, n, 'atb')
      print *, "|Z'*inv(Z')*x-x|=",sum(abs((B(1:n,1:ncols)-C(1:n,1:ncols))));
      
      ! Test Y.

      C(1:n,1:ncols) = B(1:n,1:ncols);
      call HsvdProd(A%Y, ncols, C, n, 'ab')
      call HsvdProd(A%Y, ncols, C, n, 'aib')
      print *, "|inv(Y)*Y*x-x|  =",sum(abs((B(1:n,1:ncols)-C(1:n,1:ncols))));
      
      C(1:n,1:ncols) = B(1:n,1:ncols);
      call HsvdProd(A%Y, ncols, C, n, 'aib')
      call HsvdProd(A%Y, ncols, C, n, 'ab')
      print *, "|Y*inv(Y)*x-x|  =",sum(abs((B(1:n,1:ncols)-C(1:n,1:ncols))));
      
      C(1:n,1:ncols) = B(1:n,1:ncols);
      call HsvdProd(A%Y, ncols, C, n, 'atb')
      call HsvdProd(A%Y, ncols, C, n, 'aitb')
      print *, "|inv(Y')*Y'*x-x|=",sum(abs((B(1:n,1:ncols)-C(1:n,1:ncols))));
      
      C(1:n,1:ncols) = B(1:n,1:ncols);
      call HsvdProd(A%Y, ncols, C, n, 'aitb')
      call HsvdProd(A%Y, ncols, C, n, 'atb')
      print *, "|Y'*inv(Y')*x-x|=",sum(abs((B(1:n,1:ncols)-C(1:n,1:ncols))));
      
      ! Test A.

      call EigenmatProd(A, ncols, B, n, C, n, shift, 'ab')
      call EigenmatProd(A, ncols, C, n, D, n, shift, 'aib')
      print *, "|inv(A)*A*x-x|  =",sum(abs((B(1:n,1:ncols)-D(1:n,1:ncols))))
      
      call EigenmatProd(A, ncols, B, n, C, n, shift, 'aib')
      call EigenmatProd(A, ncols, C, n, D, n, shift, 'ab')
      print *, "|A*inv(A)*x-x|  =",sum(abs((B(1:n,1:ncols)-D(1:n,1:ncols))))
      
      call EigenmatProd(A, ncols, B, n, C, n, shift, 'atb')
      call EigenmatProd(A, ncols, C, n, D, n, shift, 'aitb')
      print *, "|inv(A')*A'*x-x|=",sum(abs((B(1:n,1:ncols)-D(1:n,1:ncols))))
      
      call EigenmatProd(A, ncols, B, n, C, n, shift, 'aitb')
      call EigenmatProd(A, ncols, C, n, D, n, shift, 'atb')
      print *, "|A'*inv(A')*x-x|=",sum(abs((B(1:n,1:ncols)-D(1:n,1:ncols))))
      
      ! Test the eigenvector calculations by computing all
      ! left and right eigenvectors and their residuals.

      j = 1
      do while(j .le. n)
         call EigenmatVecs(A, j, ev(j), CC(1:n,j), DD(1:n,j), cond, 'b')
         if (A%type(j) .eq. 1 .or. A%type(j) .lt. -1 ) then
            ev1(j) = 0.0D0
            j = j + 1
         elseif (A%type(j) .eq. -1) then
            ev1(j) = A%eig(j)
            j = j + 1
         else
            CC(1:n, j+1) = conjg(CC(1:n,j))
            DD(1:n, j+1) = conjg(DD(1:n,j))
            ev(j+1) = conjg(ev(j))
            ev1(j:j+1) = 0
            j = j + 2
         end if
      end do

      ! Right eigensystem.

      C = real(CC, wp)
      D = aimag(CC)
      call EigenmatProd(A, n, C, n, B, n, 0.0_wp, 'ab')
      call EigenmatProd(A, n, D, n, C, n, 0.0_wp, 'ab')
      BB = cmplx(B, C, wp)

      norm = sum(abs(ev(1)*CC(1:n,1)-BB(1:n,1)))
      do j = 2, n
         norm = norm + sum(abs(ev1(j)*CC(1:n,j-1)+ev(j)*CC(1:n,j)-BB(1:n,j)))
      end do
      print *, "|A*REV-REV*EV|  =",norm
      
      ! Left eigensystem.
      C = real(DD, wp)
      D = aimag(DD)
      call EigenmatProd(A, n, C, n, B, n, 0.0_wp, 'atb')
      call EigenmatProd(A, n, D, n, C, n, 0.0_wp, 'atb')
      BB = cmplx(B, C, wp)
      
      norm = sum(abs(conjg(ev(n))*DD(1:n,n)-BB(1:n,n)))
      do j = 1, n-1
         norm = norm + sum(abs(ev1(j+1)*DD(1:n,j+1)+conjg(ev(j))*DD(1:n,j)-BB(1:n,j)))
      end do
      print *, "|A'*LEV-LEV*EV'|=",norm

   end do

end program testeigentest


*  This program tests the package Eigentest.  It runs 64 testcases
*  probing various aspects of the package, as described below.  The
*  numbers in the output should be within one or two orders of
*  magnitude of the rounding unit.

      program testet
      implicit none


*     These values may be adjusted at compile time to change the
*     tests, which involve computing (A - shift*I) op B,
*     where op is *, '*,  \, or '\ (' denoting transpose).

      integer nmax, n, ncols, spinrand
      double precision shift

*     Maximum order of A.
      parameter (nmax = 50)

*     Order of A.
      parameter (n = 10)

*     Number of columns in B.
      parameter (ncols = 3)

*     A shift.
      parameter (shift = -1.0D0)

*     Number of times to call the random number generator initially.
      parameter (spinrand = 0)


      integer iem(2*nmax+12)
      double precision fem(8*nmax)

      double precision B(n, n), C(n, n), D(n, n)

      integer i, j, cs, ztype, ytype, atype
      integer idy, idz, nblk, i1, i2
      double precision  norm, cond, duni, temp
      double complex eig(n), X(n, n), Y(n, n), ev1(n)
      integer ia(nmax), yip, zip, yfp, zfp
      double precision fa(nmax)


*     Spin the random number generator.

      do i = 1,spinrand
         temp = duni()
      end do

*     Initialize pointers to Y, Z.
                         
      yip = nmax + 3
      zip = nmax + 8
      yfp = nmax + 1
      zfp = 4*nmax + 1

*     Loop on test case number.

      do cs = 0, 63
         print *, "test case :", cs

*     Set up the test case.

*     Choose the type of the outer hsvd Y.  It can be an identity
*     ( ytype == 0) or a random hsvd (ytype == 1).

         ytype = mod(cs, 2)

*     Choose the typep of the inner hsvd Z.
*
*     If ztype == 0, Z is an identity.
*     If ztype == 1, Z has two blocks: (1:6)(7:10).
*     The first block is identity.
*     If ztype == 2, Z has three blocks: (1:2)(3:8)(9:10),
*       and the second block is an indentity.
*     If ztype == 3, Z has two blocks: (1:4)(5::10),
*       and the second block is identity.         

         ztype = mod(cs/2, 4)

*     Chose the types of eigenvalues.
*
*     If atype == 0, all eigenvalues are real.
*     If atype == 1, all eigenvalues are complex.
*     If atype == 2, eigenvalues 1,2,5,6,9,10 are complex.
*     If atype == 3, eigenvalues 3,4,7,8 are complex.
         atype = mod(cs/8, 8)

*     Initialize A.

         if (ytype .eq. 0) then
            idy = 1
         else
            idy = 0
         end if
         
         if (ztype .eq. 0) then
            idz = 1
            nblk = 0
         elseif (ztype .eq. 1) then

*     If ztype == 1, Z has two blocks: (1:6)(7:10).
*     The first block is identity.

            nblk = 2
            ia(1) = 1
            ia(2) = -7
            ia(3) = n+1
            idz = 0
         elseif (ztype .eq. 2) then

*      If ztype == 2, Z has three blocks: (1:2)(3:8)(9:10).
*      The second block is indentity.  

            nblk = 3
            ia(1) = 1
            ia(2) = 3
            ia(3) = -9
            ia(4) = n+1
            idz = 0
         elseif (ztype .eq. 3) then

*       If ztype == 3, Z has two blocks: (1:4)(5::10).
*       The second block is identity.

            nblk = 2
            ia(1) = 1
            ia(2) = 5
            ia(3) = -n-1         
            idz = 0
         end if

         call eminit(nmax, n, iem, fem, nblk, idy, idz, ia, fa, 'setup')

*     Generate eigenvalues.

         do i = 1, n
            fa(i) = duni() + 0.5D0
         end do
         if (atype .eq. 0) then
            do i = 1, n
               ia(i) = 1
            end do
         elseif (atype .eq. 1) then
            do i = 1, n, 2
               ia(i) = 2
               ia(i+1) = 3
            end do
         elseif (atype .eq. 2) then
            ia(1) = 2
            ia(2) = 3
            ia(3) = 1
            ia(4) = 1
            ia(5) = 2
            ia(6) = 3
            ia(7) = 1
            ia(8) = 1
            ia(9) = 2
            ia(10) = 3
         elseif (atype .eq. 3) then
            ia(1) = 1
            ia(2) = 1
            ia(3) = 2
            ia(4) = 3
            ia(5) = 1
            ia(6) = 1
            ia(7) = 2
            ia(8) = 3
            ia(9) = 1
            ia(10) = 1
         elseif (atype .eq. 4) then
            ia(1) = -n
            do i = 2, n
               ia(i) = -1
            end do
         elseif (atype .eq. 5) then
            ia(1)  =  1
            ia(2)  =  1
            ia(3)  = -2
            ia(4)  = -1
            ia(5)  =  1
            ia(6)  =  1
            ia(7)  = -3
            ia(8)  = -1
            ia(9)  = -1
            ia(10) =  1
         elseif (atype .eq. 6) then
            ia(1)  = -3
            ia(2)  = -1
            ia(3)  = -1
            ia(4)  =  1
            ia(5)  =  1
            ia(6)  =  1
            ia(7)  =  1
            ia(8)  = -3
            ia(9)  = -1
            ia(10) = -1
         elseif (atype .eq. 7) then
            ia(1)  =  1
            ia(2)  = -3
            ia(3)  = -1
            ia(4)  = -1
            ia(5)  =  1
            ia(6)  =  2
            ia(7)  =  3
            ia(8)  =  2
            ia(9)  =  3
            ia(10) =  1

         end if

         call eminit(nmax, n, iem, fem, nblk, idy, idz, ia, fa, 'eig')

         print *, "A type = ", iem(3), iem(4), iem(5), iem(6), iem(7)
     &                       , iem(8),iem(9),iem(10), iem(11), iem(12)

*     Generate Y.

         if (ytype .eq. 1) then
            do i = 1, n
               fa(i) = duni() + 1.0D0
            end do
            
            call eminit
     &           (nmax, n, iem, fem, nblk, idy, idz, ia, fa, 'ysig')

            do i = 1, n
               fa(i) = duni() - 0.5D0
            end do      
            call hscal(n, fa)
            call eminit(nmax, n, iem, fem, nblk, idy, idz, ia, fa, 'yu')
            
            do i = 1, n
               fa(i) = duni() - 0.5D0
            end do      
            call hscal(n, fa)
            call eminit(nmax, n, iem, fem, nblk, idy, idz, ia, fa, 'yv')
            
         end if
         
         print *, "Y bs = ", iem(yip+3), iem(yip+4)

*     Generate Z.

         if (ztype .gt. 0) then
            do i = 1, n
               fa(i) = duni() + 1.0D0
            end do

            call eminit
     &           (nmax, n, iem, fem, nblk, idy, idz, ia, fa, 'zsig')

            do i = 1, nblk
               i1 = abs(iem(zip+2+i))
               i2 = iem(zip+i+3) - 1
               if (i2 .gt. 0) then
                  do j = i1, i2
                     fa(i) = duni() - 0.5D0
                  end do      
                  call hscal(i2-i1+1, fa(i1))
               end if
            end do
            call eminit(nmax, n, iem, fem, nblk, idy, idz, ia, fa, 'zu')
                  
            do i = 1, nblk
               i1 = abs(iem(zip+i+2))
               i2 = iem(zip+i+3) - 1
               if (i2 .gt. 0) then
                  do j = i1, i2
                     fa(i) = duni() - 0.5D0
                  end do      
                  call hscal(i2-i1+1, fa(i1))
               end if
            end do
            call eminit(nmax, n, iem, fem, nblk, idy, idz, ia, fa, 'zv')

         end if

         if (iem(zip+2) .eq. 3) then
            print *, "Z bs = ", iem(zip+3), iem(zip+4), iem(zip+5), 
     &                          iem(zip+6)
         elseif (iem(zip+2) .eq. 2) then
            print *, "Z bs = ", iem(zip+3), iem(zip+4), iem(zip+5)
         elseif (iem(zip+2) .eq. 1) then
            print *, "Z bs = ", iem(zip+3), iem(zip+4)
         endif

*     Generate B for A, Y, and Z to operate on.

      do j = 1, ncols   
         do i = 1, n
            B(i, j) = duni() - 0.5D0
            C(i, j) = B(i, j)
         end do
      end do 

*
*     Test  Z.
*
      call hsvdpr(iem(zip), fem(zfp), ncols, C, n, 'ab')
      call hsvdpr(iem(zip), fem(zfp), ncols, C, n, 'aib')
      
      norm  = 0.0D0
      do j = 1, ncols
         do i = 1, n
            norm = norm + abs(B(i, j) - C(i, j))
            C(i, j) = B(i, j)
         end do
      end do
      
      print *, "|Z*inv(Z)*x-x|    = ", norm

      call hsvdpr(iem(zip), fem(zfp), ncols, C, n, 'aib')
      call hsvdpr(iem(zip), fem(zfp), ncols, C, n, 'ab')
      
      norm  = 0.0D0
      do j = 1, ncols
         do i = 1, n
            norm = norm + abs(B(i, j) - C(i, j))
            C(i, j) = B(i, j)
         end do
      end do

      
      print *, "|inv(Z)*Z*x-x|    = ", norm

      call hsvdpr(iem(zip), fem(zfp), ncols, C, n, 'atb')
      call hsvdpr(iem(zip), fem(zfp), ncols, C, n, 'aitb')
      
      norm  = 0.0D0
      do j = 1, ncols
         do i = 1, n
            norm = norm + abs(B(i, j) - C(i, j))
            C(i, j) = B(i, j)
         end do
      end do
      print *, "|Z'*inv(Z')*x-x|  = ", norm

      call hsvdpr(iem(zip), fem(zfp), ncols, C, n, 'aitb')
      call hsvdpr(iem(zip), fem(zfp), ncols, C, n, 'atb')
      
      norm  = 0.0D0
      do j = 1, ncols
         do i = 1, n
            norm = norm + abs(B(i, j) - C(i, j))
            C(i, j) = B(i, j)
         end do
      end do      
      print *, "|inv(Z')*Z'*x-x|  = ", norm

*
*     Test Y.
*
      call hsvdpr(iem(yip), fem(yfp), ncols, C, n, 'ab')
      call hsvdpr(iem(yip), fem(yfp), ncols, C, n, 'aib')
      
      norm  = 0.0D0
      do j = 1, ncols
         do i = 1, n
            norm = norm + abs(B(i, j) - C(i, j))
            C(i, j) = B(i, j)
         end do
      end do
      
      print *, "|Y*inv(Y)*x-x|    = ", norm

      call hsvdpr(iem(yip), fem(yfp), ncols, C, n, 'aib')
      call hsvdpr(iem(yip), fem(yfp), ncols, C, n, 'ab')
      
      norm  = 0.0D0
      do j = 1, ncols
         do i = 1, n
            norm = norm + abs(B(i, j) - C(i, j))
            C(i, j) = B(i, j)
         end do
      end do      
      print *, "|inv(Y)*Y*x-x|    = ", norm

      call hsvdpr(iem(yip), fem(yfp), ncols, C, n, 'atb')
      call hsvdpr(iem(yip), fem(yfp), ncols, C, n, 'aitb')
      
      norm  = 0.0D0
      do j = 1, ncols
         do i = 1, n
            norm = norm + abs(B(i, j) - C(i, j))
            C(i, j) = B(i, j)
         end do
      end do
      
      print *, "|Y'*inv(Y')*x-x|  = ", norm

      call hsvdpr(iem(yip), fem(yfp), ncols, C, n, 'aitb')
      call hsvdpr(iem(yip), fem(yfp), ncols, C, n, 'atb')
      
      norm  = 0.0D0
      do j = 1, ncols
         do i = 1, n
            norm = norm + abs(B(i, j) - C(i, j))
         end do
      end do
      
      print *, "|inv(Y')*Y'*x-x|  = ", norm

*
*     Test A.
*
      call emprod(iem, fem, ncols, B, n, C, n, shift, 'ab')
      call emprod(iem, fem, ncols, C, n, D, n, shift, 'aib')

      norm  = 0.0D0
      do j = 1, ncols
         do i = 1, n
            norm = norm + abs(B(i, j) - D(i, j))
         end do
      end do
      
      print *, "|inv(A)*A*x-x|    = ", norm

      call emprod(iem, fem, ncols, B, n, C, n, shift, 'aib')
      call emprod(iem, fem, ncols, C, n, D, n, shift, 'ab')

      norm  = 0.0D0
      do j = 1, ncols
         do i = 1, n
            norm = norm + abs(B(i, j) - D(i, j))
         end do
      end do      
      print *, "|A*inv(A)*x-x|    = ", norm

      call emprod(iem, fem, ncols, B, n, C, n, shift, 'atb')
      call emprod(iem, fem, ncols, C, n, D, n, shift, 'aitb')

      norm  = 0.0D0
      do j = 1, ncols
         do i = 1, n
            norm = norm + abs(B(i, j) - D(i, j))
         end do
      end do
      
      print *, "|inv(A')*A'*x-x|  = ", norm

      call emprod(iem, fem, ncols, B, n, C, n, shift, 'ab')
      call emprod(iem, fem, ncols, C, n, D, n, shift, 'aib')

      norm  = 0.0D0
      do j = 1, ncols
         do i = 1, n
            norm = norm + abs(B(i, j) - D(i, j))
         end do
      end do
      
      print *, "|A'*inv(A')*x-x|  = ", norm

*
*     Test the eigenvector calculations by computing all
*     left and right eigenvectors and their residuals.
*

      i = 1
      do while(i .le. n)
         call emvecs(iem, fem, i, eig(i), X(1,i), Y(1,i), cond, 'b')
         if (iem(2+i) .eq. 1 .or. iem(2+i) .lt. -1) then
            ev1(i) = 0.0D0
            i = i + 1
         elseif (iem(2+i) .eq. -1) then
            ev1(i) = fem(i)
            i = i + 1
         else
            do j = 1, n
               X(j, i+1) = conjg(X(j, i))
               Y(j, i+1) = conjg(Y(j, i))
            end do
            ev1(i) = 0.0D0
            ev1(i+1) = 0.0D0
            eig(i+1) = conjg(eig(i))
            i = i + 2
         end if
      end do

*     Right eigensystem.

      do j = 1, n
         do i = 1, n
            C(i,j) = dreal(X(i,j))
            D(i,j) = dimag(X(i,j))
         end do
      end do

      call emprod(iem, fem, n, C, n, B, n, 0.0D0, 'ab')
      call emprod(iem, fem, n, D, n, C, n, 0.0D0, 'ab')
      
      norm = 0.0D0
      do i = 1, n
            norm = norm + abs(eig(1)*X(i,1)-dcmplx(B(i,1),C(i,1)))
      end do
      do j = 2, n
         do i = 1, n
            norm = norm + abs(ev1(j)*X(i,j-1)+eig(j)*X(i,j) 
     &                        -dcmplx(B(i,j),C(i,j)))
         end do
      end do
      print *, '|A*REV - REV*EV|  = ', norm
      
*     Left eigensystem.

      do j = 1, n
         do i = 1, n
            C(i,j) = dreal(Y(i,j))
            D(i,j) = dimag(Y(i,j))
         end do
      end do
      
      call emprod(iem, fem, n, C, n, B, n, 0.0D0, 'atb')
      call emprod(iem, fem, n, D, n, C, n, 0.0D0, 'atb')
      
      norm = 0.0D0
      do j = 1, n-1
         eig(j) = conjg(eig(j))
         do i = 1, n
            norm = norm + abs(ev1(j+1)*Y(i,j+1) + eig(j)*Y(i,j)
     &                        -dcmplx(B(i,j),C(i,j)))
         end do
      end do
      eig(n) = conjg(eig(n))
      do i = 1, n
         norm = norm + abs(eig(n)*Y(i,n)-dcmplx(B(i,n),C(i,n)))
      end do

      print *, "|A*REV - REV*EV'| = ", norm
      end do
      stop
      end

*************************************************************************
*                            subroutines                                *
*************************************************************************
      
      DOUBLE PRECISION FUNCTION DUNI()
C***BEGIN PROLOGUE  DUNI
C***DATE WRITTEN   880714 (YYMMDD)
C***REVISION DATE  880714 (YYMMDD)
C***CATEGORY NO.  L6A21
C***KEYWORDS  RANDOM NUMBERS, UNIFORM RANDOM NUMBERS
C***AUTHOR    KAHANER, DAVID, SCIENTIFIC COMPUTING DIVISION, NBS
C             MARSAGLIA, GEORGE, SUPERCOMPUTER RES. INST., FLORIDA ST. U.
C
C***PURPOSE  THIS ROUTINE GENERATES DOUBLE PRECISION UNIFORM
C             RANDOM NUMBERS ON [0,1)
C***DESCRIPTION
C        COMPUTES DOUBLE PRECISION UNIFORM NUMBERS ON [0,1).
C           FROM THE BOOK, "NUMERICAL METHODS AND SOFTWARE" BY
C                D. KAHANER, C. MOLER, S. NASH
C                PRENTICE HALL, 1988
C
C       USAGE: 
C              TO INITIALIZE THE GENERATOR
C                   USEED = DUSTAR(ISEED)
C               WHERE: ISEED IS ANY NONZERO INTEGER
C                  WILL RETURN FLOATING POINT VALUE OF ISEED.
C
C               SUBSEQUENTLY
C                       U = DUNI()
C                  WILL RETURN A REAL UNIFORM ON [0,1)
C
C                ONE INITIALIZATION IS NECESSARY, BUT ANY NUMBER OF EVALUATIONS 
C                  OF DUNI IN ANY ORDER, ARE ALLOWED.
C
C           NOTE: DEPENDING UPON THE VALUE OF K (SEE BELOW), THE OUTPUT
C                       OF DUNI MAY DIFFER FROM ONE MACHINE TO ANOTHER.
C
C           TYPICAL USAGE: 
C
C               DOUBLE PRECISION U,DUNI,DUSTAR,USEED
C               INTEGER ISEED 
CC                 SET SEED
C               ISEED = 305
C               USEED = DUSTAR(ISEED)
C               DO 1 I = 1,1000
C                   U = DUNI()
C             1 CONTINUE
CC                 NOTE: IF K=47 (THE DEFAULT, SEE BELOW) THE OUTPUT VALUE OF
CC                           U WILL BE 0.812053811384E-01...
C               WRITE(*,*) U
C               END 
C
C          NOTE ON PORTABILITY: USERS CAN CHOOSE TO RUN DUNI IN ITS DEFAULT
C               MODE (REQUIRING NO USER ACTION) WHICH WILL GENERATE THE SAME
C               SEQUENCE OF NUMBERS ON ANY COMPUTER SUPPORTING FLOATING POINT
C               NUMBERS WITH AT LEAST 47 BIT MANTISSAS, OR IN A MODE THAT
C               WILL GENERATE NUMBERS WITH A LONGER PERIOD ON COMPUTERS WITH
C               LARGER MANTISSAS.
C          TO EXERCISE THIS OPTION:  B E F O R E  INVOKING DUSTAR INSERT
C               THE INSTRUCTION        UBITS = DUNIB(K)      K >= 47
C               WHERE K IS THE NUMBER OF BITS IN THE MANTISSA OF YOUR FLOATING
C               POINT WORD (K=96 FOR CRAY, CYBER 205). DUNIB RETURNS THE
C               FLOATING POINT VALUE OF K THAT IT ACTUALLY USED.
C                    K INPUT AS .LE. 47, THEN UBITS=47.
C                    K INPUT AS .GT. 47, THEN UBITS=FLOAT(K)
C               IF K>47 THE SEQUENCE OF NUMBERS GENERATED BY DUNI MAY DIFFER
C               FROM ONE COMPUTER TO ANOTHER.
C
C
C
C***REFERENCES  MARSAGLIA G., "COMMENTS ON THE PERFECT UNIFORM RANDOM 
C                 NUMBER GENERATOR", UNPUBLISHED NOTES, WASH S. U.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE DUNI
      DOUBLE PRECISION CSAVE,CD,CM
      PARAMETER(
     *    CSAVE=0.9162596898123D+13/0.140737488355328D+15,
     *    CD=0.76543212345678D+14/0.140737488355328D+15,
     *    CM=0.140737488355213D+15/0.140737488355328D+15)
C                            2**47=0.140737488355328D+15
      DOUBLE PRECISION U(17),S,T,DUSTAR,C,DUNIB
      INTEGER I,J,II,JJ,K,KK,I1,J1,K1,L1,M1,ISEED 
C
      SAVE U,I,J,K,C
C      LOAD DATA ARRAY IN CASE USER FORGETS TO INITIALIZE.
C      THIS ARRAY IS THE RESULT OF CALLING DUNI 100000 TIMES
C         WITH ISEED=305 AND K=96.
      DATA U/
     *0.471960981577884755837789724978D+00,
     *0.930323453205669578433639632431D+00,
     *0.110161790933730836587127944899D+00,
     *0.571501996273139518362638757010D-01,
     *0.402467554779738266237538503137D+00,
     *0.451181953427459489458279456915D+00,
     *0.296076152342721102174129954053D+00,
     *0.128202189325888116466879622359D-01,
     *0.314274693850973603980853259266D+00,
     *0.335521366752294932468163594171D-02,
     *0.488685045200439371607850367840D+00,
     *0.195470426865656758693860613516D+00,
     *0.864162706791773556901599326053D+00,
     *0.335505955815259203596381170316D+00,
     *0.377190200199058085469526470541D+00,
     *0.400780392114818314671676525916D+00,
     *0.374224214182207466262750307281D+00/
      DATA I,J,K,C/17,5,47,CSAVE/
C
C   BASIC GENERATOR IS FIBONACCI
C
      DUNI = U(I)-U(J)
      IF(DUNI.LT.0.0D0)DUNI = DUNI+1.0D0
      U(I) = DUNI
      I = I-1
      IF(I.EQ.0)I = 17
      J = J-1
      IF(J.EQ.0)J = 17
C
C   SECOND GENERATOR IS CONGRUENTIAL
C
      C = C-CD
      IF(C.LT.0.0D0) C=C+CM
C
C   COMBINATION GENERATOR
C
      DUNI = DUNI-C 
      IF(DUNI.LT.0.0D0)DUNI = DUNI+1.0D0
      RETURN
C
      ENTRY DUSTAR(ISEED)
C
C          SET UP ...
C          CONVERT ISEED TO FOUR SMALLISH POSITIVE INTEGERS.
C
        I1 = MOD(ABS(ISEED),177)+1
        J1 = MOD(ABS(ISEED),167)+1
        K1 = MOD(ABS(ISEED),157)+1
        L1 = MOD(ABS(ISEED),147)+1
C
C              GENERATE RANDOM BIT PATTERN IN ARRAY BASED ON GIVEN SEED.
C
        DO 2 II = 1,17
          S = 0.0D0 
          T = 0.5D0 
C             DO FOR EACH OF THE BITS OF MANTISSA OF WORD
C             LOOP  OVER K BITS, WHERE K IS DEFAULTED TO 47 BUT CAN
C               BE CHANGED BY USER CALL TO DUNIB(K)
          DO 3 JJ = 1,K
                  M1 = MOD(MOD(I1*J1,179)*K1,179) 
                  I1 = J1
                  J1 = K1
                  K1 = M1
                  L1 = MOD(53*L1+1,169) 
                  IF(MOD(L1*M1,64).GE.32)S=S+T
    3             T = 0.5D0*T 
    2   U(II) = S
        DUSTAR = FLOAT(ISEED) 
        RETURN
C
      ENTRY DUNIB(KK)
        IF(KK.LE.47)THEN
             K=47
        ELSE
             K=KK
        ENDIF
        DUNIB=FLOAT(K)
      END


*     ============================================================


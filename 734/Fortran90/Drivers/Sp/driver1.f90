C--**--CH2468--734--P:RW--22:9:1999
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> test_functions.f90
    MODULE test_functions

      USE precisions, ONLY : stnd, short

      IMPLICIT NONE
      PUBLIC
!------              -------------!

      INTEGER, PARAMETER :: nargs = 10

      REAL (stnd) :: arguments(nargs)

      INTEGER :: funcno

!--------- DATA FOR TESTPACK FUNCTION     ARGAUS

      REAL (stnd), DIMENSION (15), PARAMETER :: argasy = (/ 9.000D-4, &
        4.400D-3, 1.750D-2, 5.400D-2, 1.295D-1, 2.420D-1, 3.521D-1, 3.989D-1, &
        3.521D-1, 2.420D-1, 1.295D-1, 5.400D-2, 1.750D-2, 4.400D-3, 9.000D-4/)

!--------- DATA FOR TESTPACK FUNCTION     BARD70

      REAL (stnd), DIMENSION (15), PARAMETER :: bard7y = (/ .14D0, .18D0, &
        .22D0, .25D0, .29D0, .32D0, .35D0, .39D0, .37D0, .58D0, .73D0, .96D0, &
        1.34D0, 2.10D0, 4.39D0/)

!--------- DATA FOR TESTPACK FUNCTION     CHNRSN

      REAL (stnd), DIMENSION (50), PARAMETER :: chnrsn = (/ 1.25D0, 1.40D0, &
        2.40D0, 1.40D0, 1.75D0, 1.20D0, 2.25D0, 1.20D0, 1.00D0, 1.10D0, &
        1.50D0, 1.60D0, 1.25D0, 1.25D0, 1.20D0, 1.20D0, 1.40D0, 0.50D0, &
        0.50D0, 1.25D0, 1.80D0, 0.75D0, 1.25D0, 1.40D0, 1.60D0, 2.00D0, &
        1.00D0, 1.60D0, 1.25D0, 2.75D0, 1.25D0, 1.25D0, 1.25D0, 3.00D0, &
        1.50D0, 2.00D0, 1.25D0, 1.40D0, 1.80D0, 1.50D0, 2.20D0, 1.40D0, &
        1.50D0, 1.25D0, 2.00D0, 1.50D0, 1.25D0, 1.40D0, 0.60D0, 1.50D0/)

!--------- DATA FOR TESTPACK FUNCTION     HIMM32

      REAL (stnd), DIMENSION (7), PARAMETER :: him32a = (/ 0.0D0, 4.28D-4, &
        1.0D-3, 1.61D-3, 2.09D-3, 3.48D-3, 5.25D-3/)

      REAL (stnd), DIMENSION (7), PARAMETER :: him32b = (/ 7.391D0, 1.118D1, &
        1.644D1, 1.62D1, 2.22D1, 2.402D1, 3.132D1/)

!--------- DATA FOR TESTPACK FUNCTION     KOWOSB

      REAL (stnd), DIMENSION (11), PARAMETER :: kowosu = (/ 4.0D0, 2.0D0, &
        1.0D0, 0.5D0, 0.25D0, 0.167D0, 0.125D0, 0.1D0, 0.0833D0, 0.0714D0, &
        0.0625D0/)

      REAL (stnd), DIMENSION (11), PARAMETER :: kowosy = (/ 0.1957D0, &
        0.1947D0, 0.1735D0, 0.1600D0, 0.0844D0, 0.0627D0, 0.0456D0, 0.0342D0, &
        0.0323D0, 0.0235D0, 0.0246D0/)

!--------- DATA FOR TESTPACK FUNCTION        MEYER

      REAL (stnd), DIMENSION (16), PARAMETER :: mey = (/ 3.478D4, 2.861D4, &
        2.365D4, 1.963D4, 1.637D4, 1.372D4, 1.154D4, 9.744D3, 8.261D3, &
        7.030D3, 6.005D3, 5.147D3, 4.427D3, 3.820D3, 3.307D3, 2.872D3/)

!--------- DATA FOR TESTPACK FUNCTION         ORTOIT

      REAL (stnd), DIMENSION (33), PARAMETER :: orbeta = (/ 1.0D0, 1.5D0, &
        1.0D0, 0.1D0, 1.5D0, 2.0D0, 1.0D0, 1.5D0, 3.0D0, 2.0D0, 1.0D0, 3.0D0, &
        0.1D0, 1.5D0, 0.15D0, 2.0D0, 1.0D0, 0.1D0, 3.0D0, 0.1D0, 1.2D0, 1.0D0, &
        0.1D0, 2.0D0, 1.2D0, 3.0D0, 1.5D0, 3.0D0, 2.0D0, 1.0D0, 1.2D0, 2.0D0, &
        1.0D0/)

      REAL (stnd), DIMENSION (33), PARAMETER :: od = (/ 5.0D0, 5.0D0, 5.0D0, &
        2.5D0, 6.0D0, 6.0D0, 5.0D0, 6.0D0, 10.0D0, 6.0D0, 5.0D0, 9.0D0, 2.0D0, &
        7.0D0, 2.5D0, 6.0D0, 5.0D0, 2.0D0, 9.0D0, 2.0D0, 5.0D0, 5.0D0, 2.5D0, &
        5.0D0, 6.0D0, 10.0D0, 7.0D0, 10.0D0, 6.0D0, 5.0D0, 4.0D0, 4.0D0, &
        4.0D0/)

      INTEGER, DIMENSION (50), PARAMETER :: a = (/ -31, -1, -2, -4, -6, -8, &
        -10, -12, + 11, + 13, -14, -16, + 9, -18, + 5, + 20, -21, -19, -23, &
        + 7, -25, -28, -29, -32, + 3, -33, -35, -36, + 30, -37, + 38, -39, &
        -40, -41, -44, -46, + 42, + 45, + 48, -50, + 26, + 34, -43, + 15, &
        + 17, + 24, -47, -49, -22, -27 /)

      INTEGER, DIMENSION (56), PARAMETER :: b = (/ -1, + 2, -3, + 4, -5, + 6, &
        -7, + 8, -9, + 10, -11, + 12, -13, + 14, -15, + 16, -17, + 18, -19, &
        -20, 0, + 22, + 23, -24, + 25, -26, + 27, -28, + 29, -30, + 31, -32, &
        + 33, -34, -35, + 21, -36, + 37, -38, -39, -40, + 41, -42, + 43, + 44, &
        -50, + 45, + 46, -47, -48, -49, 0, 0, 0, 0, 0 /)

!--------- DATA FOR TESTPACK FUNCTION        OSBRN1

      REAL (stnd), DIMENSION (33), PARAMETER :: osb1y = (/ .844D0, .908D0, &
        .932D0, .936D0, .925D0, .908D0, .881D0, .850D0, .818D0, .784D0, &
        .751D0, .718D0, .685D0, .658D0, .628D0, .603D0, .580D0, .558D0, &
        .538D0, .522D0, .506D0, .490D0, .478D0, .467D0, .457D0, .448D0, &
        .438D0, .431D0, .424D0, .420D0, .414D0, .411D0, .406D0/)

!--------- DATA FOR TESTPACK FUNCTION        OSBRN2

      REAL (stnd), DIMENSION (65), PARAMETER :: osb2y = (/ 1.366D0, 1.191D0, &
        1.112D0, 1.013D0, .991D0, .885D0, .831D0, .847D0, .786D0, .725D0, &
        .746D0, .679D0, .608D0, .655D0, .616D0, .606D0, .602D0, .626D0, &
        .651D0, .724D0, .649D0, .649D0, .694D0, .644D0, .624D0, .661D0, &
        .612D0, .558D0, .533D0, .495D0, .50D0, .423D0, .395D0, .375D0, .372D0, &
        .391D0, .396D0, .405D0, .428D0, .429D0, .523D0, .562D0, .607D0, &
        .653D0, .672D0, .708D0, .633D0, .668D0, .645D0, .632D0, .591D0, &
        .559D0, .597D0, .625D0, .739D0, .710D0, .729D0, .720D0, .636D0, &
        .581D0, .428D0, .292D0, .162D0, .098D0, .054D0/)

    CONTAINS

      SUBROUTINE functions(x,f,g,ifg)
        USE reals, ONLY : c100, c20, c200, c36, c40, c400, c600, eight, fifth, &
          five, four, one, six, ten, tenth, three, two, zero
        USE supp_codes, ONLY : justf, justg, noforg, ok
        USE num_constants, ONLY : machhuge


        IMPLICIT NONE
!-------------!

! ARGUMENTS:

        INTEGER (short), INTENT (INOUT) :: ifg

        REAL (stnd), INTENT (IN) :: x(:)
        REAL (stnd) :: work(2*SIZE(x))
        REAL (stnd), INTENT (OUT) :: f, g(SIZE(x))

! STATUS:

!    SYSTEM  DEPENDENCE:               NONE


! DESCRIPTION:

!          THIS TEST FUNCTION  EVALUATES ONE OF THE STANDARD TEST
!     FUNCTIONS PROVIDED WITH TESTPACK.  THE ARGUMENTS IN THE CALLING
!     SEQUENCE HAVE PRECISELY THE SAME MEANING AS IN THE ROUTINE EVALUATE_F.

!         THE TEST FUNCTION TO USE IS SELECTED THRU THE MODULE TEST_FUNCTIONS.
!     THE VALUE OF THE INTEGER, FUNCNO, SPECIFIES WHICH OF THE TEST FUNCTIONS
!     IS TO BE  USED; THE FUNCTION IS CHOSEN USING A CASE CONSTRUCT.

!         SOME OF THE FUNCTIONS NEED SPECIAL ARGUMENTS (OTHER THAN THE
!     VALUE OF X); THESE ARE PROVIDED THROUGH THE MODULE TEST_FUNCTIONS. A
!     MAXIMUM OF TEN ARGUMENTS ARE PROVIDED. IF THE MAXIMUM NUMBER OF
!     ARGUMENTS IS TO BE INCREASED, THE PARAMETER NARGS SHOULD BE INCREASED.

!         ALL FUNCTION ARGUMENTS ARE REAL. INTEGER VALUES MAY BE PASSED
!     BY ASSIGNING THE INTEGER VALUE TO A REAL ARGUMENT AND THEN USING
!     NINT TO RECOVER THE INTEGER VALUE.

!         THE AMOUNT OF SPACE AVAILABLE IN THE ARRAY WORK IS ALLOCATED
!     DYNAMICALLY. THIS MEANS THAT IT DOES NOT HAVE TO
!     BE PROVIDED IN THE CALL TO FUNCTIONS OR IN THE CALL TO EVALUATE_F.

! SUBROUTINES

!     PREDEFINED FUNCTIONS : SIN, COS, TAN, ACOS, ATAN, ABS, MAX, NINT
!                            EXP, LOG, MIN, MOD, SIGN, SQRT, REAL(DBLE)

! PARAMETERS

        INTEGER, PARAMETER :: alpha = 5, beta = 14, gamma = 3

! LOCAL DECLARATIONS

        INTEGER :: i, j, k, n, i1, i2

        INTEGER (short) :: ret

        LOGICAL :: fonly, gonly

        REAL (stnd) :: pi, biggst


!--------- VARIABLES FOR THE TEST FUNCTIONS.


        REAL (stnd) :: x1, x2, x3, x4, x5, x6, g1, g2, g3, g4, g5, g6, w1, w2, &
          w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, t, ri, ti, yi, rf1, rf2

! EQUIVALENCES:         None.

! COMMON:               None.

! EXECUTION

        pi = ACOS(-one)
        biggst = LOG(machhuge)

!--------- SET LOGICAL FLAGS AND SELECT FUNCTION.

        fonly = ifg == justf
        gonly = ifg == justg
        n = SIZE(x)
        ret = ok

!     THE TEST FUNCTIONS APPEAR HERE IN ALPHABETICAL ORDER.
!        Case( 26 );    3600    BARD70
!        Case( 11 );    2100    BIGGS6
!        Case( 16 );    2600    BOX663
!        Case( 27 );    3700    CRGLVY
!        Case( 25 );    3500    ENGVL2
!        Case( 40 );    5000    PENAL1
!        Case( 41 );    5100    PENAL2
!        Case(  3 );    1300    PWSING
!        Case(  1 );    1100    ROSENB
!        Case( 24 );    3400    SCHMVT
!        Case( 81 );    9100    SARSEB ! new quadratic - not yet in Mintest
! package

!>>>>>    NOTE : IF WE SUPPOSE THAT EACH OF THESE TEST FUNCTIONS HAD
!>>>>>           BEEN CODED AS A SEPARATE ROUTINE, THEN, UNLESS
!>>>>>           OTHERWISE SPECIFIED, ALL TEST FUNCTIONS WOULD HAVE
!>>>>>           HAD AN ARGUMENT LIST AS FOLLOWS:
!>>>>>
!>>>>>                  ( X, F, G, CASE )
!>>>>>
!>>>>>           THOSE WHICH WOULD REQUIRE ADDITIONAL ARGUMENTS ARE
!>>>>>           NOTED BY GIVING A SUITABLE CALLING SEQUENCE. THIS
!>>>>>           SERVES TO DEFINE THE SPECIAL ARGUMENTS FOR THOSE TEST
!>>>>>           FUNCTIONS. SEE FOR EXAMPLE PENAL2 AT 5100.

        SELECT CASE (funcno)

!--------- TESTPACK FUNCTION     ROSENB

        CASE (1)

! 1100
          x1 = x(1)
          w1 = one - x1
          w2 = x(2) - x1*x1

          IF ( .NOT. gonly) f = c100*w2*w2 + w1*w1

          IF ( .NOT. fonly) THEN
            g(1) = -c400*w2*x1 - two*w1
            g(2) = c200*w2
          END IF



!--------- TESTPACK FUNCTION     PWSING

        CASE (3)

! 1300
          IF ( .NOT. gonly) f = zero

          IF (4*(n/4)/=n) THEN

            IF ( .NOT. fonly) g = zero

          ELSE

            DO i = 1, n/4

              j = 4*i

              w1 = x(j-3)
              w2 = x(j-2)
              w3 = x(j-1)
              w4 = x(j)

              w5 = w1 + ten*w2
              w6 = w3 - w4
              w2 = w2 - two*w3
              w3 = w2*w2*w2
              w1 = w1 - w4
              w4 = w1*w1*w1

              IF ( .NOT. gonly) f = f + w5*w5 + five*w6*w6 + w2*w3 + ten*w1*w4

              IF ( .NOT. fonly) THEN
                g(j-3) = two*w5 + c40*w4
                g(j-2) = c20*w5 + four*w3
                g(j-1) = ten*w6 - eight*w3
                g(j) = -ten*w6 - c40*w4
              END IF

            END DO

          END IF



!--------- TESTPACK FUNCTION   BIGGS6 ( X, F, G, IFG, NINT(ARGUMENTS(1)))
!--------- NINT(ARGUMENTS(1)) IS M

        CASE (11)

! 2100
          x1 = x(1)
          x2 = x(2)
          x3 = x(3)
          x4 = x(4)
          x5 = x(5)
          x6 = x(6)

          IF ( .NOT. gonly) f = zero

          IF ( .NOT. fonly) THEN
            g1 = zero
            g2 = zero
            g3 = zero
            g4 = zero
            g5 = zero
            g6 = zero
          END IF

          DO i = 1, NINT(arguments(1))
            t = REAL(i)
            ti = t/ten
            IF (MAX(-t,-ti*four,-ti*x1,-ti*x2,-ti*x5)<=biggst) THEN
              yi = EXP(-ti) - five*EXP(-t) + three*EXP(-four*ti)
              w3 = EXP(-ti*x1)
              w4 = EXP(-ti*x2)
              w5 = EXP(-ti*x5)
            ELSE
              ret = noforg
              GO TO 10
            END IF
            ri = x3*w3 - x4*w4 + x6*w5 - yi

            IF ( .NOT. gonly) f = f + ri*ri

            IF ( .NOT. fonly) THEN
              w1 = ti*ri
              g1 = g1 - w3*w1
              g2 = g2 + w4*w1
              g3 = g3 + w3*ri
              g4 = g4 - w4*ri
              g5 = g5 - w5*w1
              g6 = g6 + w5*ri
            END IF

          END DO

          IF ( .NOT. fonly) THEN
            g(1) = two*x3*g1
            g(2) = two*x4*g2
            g(3) = two*g3
            g(4) = two*g4
            g(5) = two*x6*g5
            g(6) = two*g6
          END IF



!--------- TESTPACK FUNCTION   BOX663 ( X, F, G, IFG, NINT(ARGUMENTS(1)))
!--------- NINT(ARGUMENTS(1)) IS M

        CASE (16)

! 2600
          IF ( .NOT. gonly) f = zero

          IF ( .NOT. fonly) THEN
            g1 = zero
            g2 = zero
            g3 = zero
          END IF

          DO i = 1, NINT(arguments(1))
            w2 = REAL(i)
            ti = w2/ten
            IF (MAX(-w2,-ti,-ti*x(1),-ti*x(2))<=biggst) THEN
              w3 = EXP(-ti*x(1))
              w4 = EXP(-ti*x(2))
              w5 = EXP(-ti) - EXP(-w2)
            ELSE
              ret = noforg
              GO TO 10
            END IF
            ri = w3 - w4 - w5*x(3)

            IF ( .NOT. gonly) THEN
              IF (ABS(ri)<=SQRT(machhuge-MAX(f,zero))) THEN
                f = f + ri*ri
              ELSE
                ret = noforg
                GO TO 10
              END IF
            END IF

            IF ( .NOT. fonly) THEN
              w2 = ti*ri
              g1 = g1 - w3*w2
              g2 = g2 + w4*w2
              g3 = g3 - w5*ri
            END IF

          END DO

          IF ( .NOT. fonly) THEN
            g(1) = two*g1
            g(2) = two*g2
            g(3) = two*g3
          END IF



!--------- TESTPACK FUNCTION     SCHMVT

        CASE (24)

! 3400
          x1 = x(1)
          x2 = x(2)
          x3 = x(3)

          w1 = x1 - x2
          w2 = x1 + x3

          w3 = one + w1*w1
          w4 = (pi*x2+x3)/two
          w5 = (w2/x2) - two
          IF (-w5**2<=biggst) THEN
            w6 = EXP(-w5*w5)
          ELSE
            ret = noforg
            GO TO 10
          END IF

          IF ( .NOT. gonly) f = -((one/w3)+SIN(w4)+w6)

          IF ( .NOT. fonly) THEN

            w3 = two*w1/(w3*w3)
            w4 = COS(w4)/two
            w6 = two*w5*w6/x2

            g(1) = w3 + w6
            g(2) = -w3 - pi*w4 - w6*w2/x2
            g(3) = -w4 + w6

          END IF


!--------- TESTPACK FUNCTION     ENGVL2

        CASE (25)

! 3500
          x1 = x(1)
          x2 = x(2)
          x3 = x(3)

          w1 = x1*x1
          w2 = x1*w1
          w3 = x2*x2
          w4 = x3*x3

          w5 = x3 - two
          w6 = five*x3 - x1 + one
          w7 = w1 + w3 - one

          w8 = w7 + w4
          w9 = w7 + w5*w5
          w10 = x1 + x2 + x3 - one
          w11 = x1 + x2 - x3 + one
          w12 = w2 + three*w3 + w6*w6 - c36

          IF ( .NOT. gonly) f = w8*w8 + w9*w9 + w10*w10 + w11*w11 + w12*w12

          IF ( .NOT. fonly) THEN
            w10 = w8 + w9
            g(1) = two*(two*x1*w10+two*(x1+x2)+w12*(three*w1-two*w6))
            g(2) = two*(two*x2*w10+two*(x1+x2)+six*w12*x2)
            g(3) = two*(two*(w8*x3+w5*w9)+two*x3-two+ten*w12*w6)
          END IF



!--------- TESTPACK FUNCTION     BARD70

        CASE (26)

! 3600
          x1 = x(1)
          x2 = x(2)
          x3 = x(3)

          IF ( .NOT. gonly) f = zero

          IF ( .NOT. fonly) THEN
            g1 = zero
            g2 = zero
            g3 = zero
          END IF

          DO i = 1, 15

            w1 = REAL(i)
            w2 = REAL(16-i)
            w3 = MIN(w1,w2)

            w4 = x2*w2 + x3*w3
            ri = bard7y(i) - (x1+w1/w4)
            w4 = w4*w4

            IF ( .NOT. gonly) f = f + ri*ri

            IF ( .NOT. fonly) THEN
              w4 = ri*w1/w4
              g1 = g1 - ri
              g2 = g2 + w2*w4
              g3 = g3 + w3*w4
            END IF
          END DO

          IF ( .NOT. fonly) THEN
            g(1) = g1*two
            g(2) = g2*two
            g(3) = g3*two
          END IF


!--------- TESTPACK FUNCTION     CRGLVY

        CASE (27)

! 3700
          x1 = x(1)
          x2 = x(2)
          x3 = x(3)
          x4 = x(4)

          w1 = x2 - x3
          w2 = x3 - x4
          w3 = x4 - one

          IF (x1<=biggst) THEN
            w4 = EXP(x1)
          ELSE
            ret = noforg
            GO TO 10
          END IF
          w5 = w4 - x2
          w6 = TAN(w2)

          IF ( .NOT. gonly) f = w5**4 + c100*w1**6 + w6**4 + x1**8 + w3*w3

          IF ( .NOT. fonly) THEN

            w2 = COS(w2)
            w5 = four*w5**3
            w1 = c600*w1**5
            w6 = four*w6**3/(w2*w2)

            g(1) = w4*w5 + eight*x1**7
            g(2) = -w5 + w1
            g(3) = -w1 + w6
            g(4) = -w6 + two*w3
          END IF



!--------- TESTPACK FUNCTION     PENAL1 ( X, F, G, IFG,
!                                         ARGUMENTS(1), ARGUMENTS(2)      )
!--------- ARGUMENTS(1) IS A
!--------- ARGUMENTS(2) IS B

        CASE (40)

! 5000
          rf1 = arguments(1)
          rf2 = arguments(2)

          w1 = -one/four
          w2 = zero

          DO j = 1, n
            w3 = x(j)
            w1 = w1 + w3*w3
            w3 = w3 - one
            w2 = w2 + w3*w3
          END DO

          IF ( .NOT. gonly) f = rf1*w2 + rf2*w1*w1

          IF ( .NOT. fonly) THEN
            w1 = four*rf2*w1
            w2 = two*rf1
            DO j = 1, n
              w3 = x(j)
              g(j) = w2*(w3-one) + w3*w1
            END DO
          END IF



!--------- TESTPACK FUNCTION     PENAL2 ( X, F, G, IFG )
!--------- ARGUMENTS(1) IS A
!--------- ARGUMENTS(2) IS B

        CASE (41)

! 5100
          rf1 = arguments(1)
          rf2 = arguments(2)

          IF (SIZE(work)<2*n) THEN
            f = zero
            g = zero
            GO TO 10
          END IF

          w1 = EXP(tenth)
          w2 = EXP(-tenth)
          w3 = zero

          i1 = 0
          i2 = n

          DO k = 1, n
            w4 = x(k)
            w3 = w3 + REAL(n-k+1)*w4*w4
            IF (tenth*w4<=biggst) THEN
              w5 = EXP(tenth*w4)
            ELSE
              ret = noforg
              GO TO 10
            END IF

            IF (k==1) THEN
              w6 = zero
              w7 = one

            ELSE
              w7 = w9*w1
              w10 = w5 + w8 - (w7+w9)
              w11 = w5 - w2

              IF ( .NOT. fonly) THEN
                work(i1+k) = w10
                work(i2+k) = w11
              END IF

              IF ( .NOT. gonly) w6 = w6 + w10*w10 + w11*w11

            END IF

            w8 = w5
            w9 = w7

          END DO

          w1 = x(1) - fifth
          w2 = w3 - one

          IF ( .NOT. gonly) f = rf1*w6 + rf2*(w1*w1+w2*w2)

          IF ( .NOT. fonly) THEN
            w3 = fifth*rf1
            w2 = four*rf2*w2

            DO k = 1, n

!          --NOTE THAT W8 DOES NOT NEED TO BE PRE-DEFINED WHEN K = 1.

              w4 = x(k)
              IF (tenth*w4<=biggst) THEN
                w5 = EXP(tenth*w4)
              ELSE
                ret = noforg
                GO TO 10
              END IF
              IF (k/=1) THEN
                w6 = w8
                w7 = work(i2+k)
              END IF

              IF (k<n) THEN
                w8 = work(i1+k+1)
                IF (k==1) THEN
                  g(1) = w3*w5*(w8) + w2*w4*REAL(n) + w1*two*rf2

                ELSE
                  g(k) = w3*w5*(w6+w7+w8) + w2*w4*REAL(n-k+1)

                END IF

              ELSE
                g(n) = w3*w5*(w6+w7) + w2*w4

              END IF

            END DO

          END IF


!----    SARSEB (Test problem 1 from Lootsma p 67.)

        CASE (81)

          x1 = x(1)
          x2 = x(2)
          x3 = x(3)
          x4 = x(4)

          w1 = two*x1 - x3 - one
          w2 = x2 - three
          w3 = two*x3 + x4 + one
          w4 = x4 - one

          IF ( .NOT. gonly) f = x1*w1 + x2*w2 + x3*w3 + x4*w4

          IF ( .NOT. fonly) THEN
            g(1) = w1 + two*x1
            g(2) = w2 + x2
            g(3) = w3 + two*x3 - x1
            g(4) = w4 + x4 + x3
          END IF

        CASE DEFAULT
          GO TO 10
        END SELECT


! EXIT:

10      ifg = ret

        RETURN

! FORMATS: None.

      END SUBROUTINE

    END
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> toms.f90
    PROGRAM test_minimize_f
      USE minimize, ONLY : done, long, minimizestate, minimize_f
      USE min_codes, ONLY : available, ident, normal, revcommfandg, &
        revcommrestart, sumform
      USE support, ONLY : evaluate_f
      USE supp_codes, ONLY : analytic, both
      USE supp_defs, ONLY : defaultevalstate
      USE test_functions, ONLY : arguments, funcno, functions, short, stnd
      USE reals, ONLY : c13, c1_m5, c20, half, one, ten, three, two, zero
      USE general, ONLY : cpusecs, my_date_time
      USE true_false, ONLY : false, true
      USE min_states, ONLY : tracelist

! DESCRIPTION:

!     This is a routine provided to test  the minimization routine minimize_f
!     after installation on a particular system.  It calls  minimize_f
!     to minimize a collection of 10 test functions which are provided
!     in the module Test_Functions.

!     It also serves as a model to illustrate the use of some of the
!     features of the minimization algorithm implementation.

!     For an example of the coding of a test function, see the routine
!     Functions in the file test_functions.f90.

!     Each function is minimized several times.  Tests involve analytic
!     and differenced derivatives and use both forward and reverse
!     communication. Both the conjugate gradient and quasi-Newton codes
!     are tried as well as the Nocedal updates.  These are the standard
!     tests performed with the namelist file provided.  A total of 73
!     test runs will be executed.

!     Note that test #8 runs a function with n = 500 using the quasi-Newton
!     method, if enough memory is available.  On a Sun Sparc10 (which runs at
!     about 5 megaflops), this function takes about 130 seconds to minimize at
!     double precision.  On a Sparc1 (about 0.7 megaflops), it takes over
!     660 secs.

!     If it is preferred that this long test not be run, simply set the
!     single occurrence of the argument  run8  in the file namelist.inp
!     to false.  In this case, only 72 test runs will be made instead of 73.
!     The time required on a SPARCstation 1 is then reduced to about 110 secs.

!----- Summary of Tests.

!   Standard tests:
!     1. Forward call, analytic derivatives, method = Available.
!     2. Reverse call, analytic derivatives, method = Available.
!     3. Forward call, finite differences,   method = Available.
!     4. Forward call, derivative testing,   method = Available.
!     5. Forward call, analytic derivatives, method = FixTerms.
!     6. Reverse call, analytic derivatives, method = QN.
!     7. Forward call, analytic derivatives, method = Available,Nocedal updates.
!     8. Forward call, analytic derivatives, method = Available, function #11
!        only (dimension large).
!     9. Forward call, analytic derivatives, method = Dynamic, expense = 1,
!        function #11 only (dimension large).
!    10. Forward call, analytic derivatives, method = Dynamic, expense = 300,
!        function #11 only (dimension large).

!     The auxiliary tests listed below provide extra testing, but require a more
!     extensive namelist file.  They also use a quadratic function that is not
!     used in the standard set of tests.  The auxiliary tests are *not* used
!     in testing the minimization routine after installation on a new system.

!   Auxiliary tests:
!    11. Forward call, analytic derivatives, method = ConMin.
!    12. Forward call, analytic derivatives, method = Variable, Nocedal updates.
!    13. Forward call, analytic derivatives, method = SD,
!        exact line search with quadratic function sarseb only.
!    14. Forward call, analytic derivatives, method = CG,
!        exact line search with quadratic function sarseb only.

! SUBROUTINES:

!     minimize_f
!     evaluate_f

!     My_Date_Time
!     CpuSecs

!     functions
!        bard70, biggs6, box663, crglvy, engvl2  |a subset of a full coll-
!        penal1, penal2, pwsing, rosenb, schmvt  |ection of test functions
!        sarseb

! PARAMETERS:


      IMPLICIT NONE
!-------------!

! maximal dimension
! number of tests
! number of problems
      INTEGER, PARAMETER :: mxn = 500, tests = 14, nprobs = 12, &
        ttests = tests*nprobs

! interactive input
! interactive output
! input from data file
      INTEGER (short), PARAMETER :: inpt = 7,  outpt = 8 ! output from tests
      INTEGER (short):: trminp, trmout
      INTEGER :: i1mach

! print first/last points
      INTEGER (long), PARAMETER :: freq = 1000, memory = 130000 ! available memory

! LOCAL DECLARATIONS:

!---- Places to hold some statistics and stuff from the tests.

      INTEGER :: icnts(ttests), fcnts(ttests), fnct, grct, indx(nprobs), &
        compnt(nprobs)

      REAL (stnd) :: fvals(ttests), acctim, acc

!---- Various declarations needed to run the tests.

      TYPE (tracelist) :: dotraces

      INTEGER :: n, error, i, m, contrl, from, to, pdone, tfncs, titers, &
        mprobs

      INTEGER, PARAMETER :: DIM(nprobs) = (/ 3, 6, 3, 4, 3, 10, 10, 40, 2, 3, &
        500, 4 /), ifnc(nprobs) = (/ 26, 11, 16, 27, 25, 40, 41, 3, 1, 24, 3, &
        81 /)

      REAL (stnd), PARAMETER :: rpar1(nprobs) = (/ zero, c13, ten, zero, zero, &
        c1_m5, c1_m5, zero, zero, zero, zero, zero/), rpar2(nprobs) = (/ zero, &
        zero, zero, zero, zero, one, one, zero, zero, zero, zero, zero/)

      REAL (stnd) :: x(mxn), g(mxn), dererr(nprobs), fx, time, averrs(nprobs)

!---- Declarations for remaining Minimize_f arguments.

      INTEGER (short) :: status, state, case, derv, meth, update, set_h0, &
        expens
      INTEGER (long) :: mterms, chkpoint
      LOGICAL :: dotrace, exact
      CHARACTER (30) :: chkfile

      TYPE (minimizestate) :: c

      LOGICAL :: firsttime, run8

      CHARACTER (41) :: date
!     Output identification.
      CHARACTER (60) :: title(tests) = (/ &
        ' analytic mode, forward calls; meth= Available.             ', &
        ' analytic mode, reverse calls; meth= Available.             ', &
        ' differencing,  forward calls; meth= Available.             ', &
        ' testing mode,  forward calls; meth= Available.             ', &
        ' analytic mode, forward calls; meth= FixTerms.              ', &
        ' analytic mode, reverse calls; meth= QN.                    ', &
        ' analytic mode, forward calls, meth= Available, Nocedal ups.', &
        ' analytic mode, forward calls; meth= Available, big n.      ', &
        ' analytic mode, forward calls; meth= Dynamic,   big n.      ', &
        ' analytic mode, forward calls; meth= Dynamic, big expense,n.', &
        ' analytic mode, forward calls, meth= ConMin.                ', &
        ' analytic mode, forward calls, meth= Variable, Nocedal ups. ', &
        ' analytic mode, forward calls, meth= SD, exact line search. ', &
        ' analytic mode, forward calls, meth= CG, exact line search. '/)

      NAMELIST /args/state, dotrace, chkpoint, chkfile, exact, expens, derv, &
        meth, mterms, update, set_h0, run8

! EXECUTION:

      acctim = zero !---- Initialize timing.
      trminp = i1mach(1)
      trmout = i1mach(2)

      OPEN (inpt,file='data1.7',action='READ')
      OPEN (outpt,file='res1.8')
      CALL my_date_time(date)
      WRITE (outpt,90020) ' starting test at ', date

      error = 0 ! Initialize counts.
      pdone = 0
      tfncs = 0
      titers = 0

      dererr = zero
      compnt = 0
      indx = 0
      averrs = zero

      fvals = zero
      fcnts = 0
      icnts = 0

      mprobs = 10 ! Do minimizations. Run each of the test types.
      to = mprobs
      contrl = 0

!   Initialize accuracy.

      IF (stnd==SELECTED_REAL_KIND(12)) THEN
        acc = 5.0D-04 ! accuracy for normal precision
      ELSE IF (stnd==SELECTED_REAL_KIND(6)) THEN
        acc = 4.0D-03 ! accuracy for low precision
      END IF

TESTSET: DO m = 1, tests

!       Default settings for some of the arguments of Minimize_f.
!       These arguments can be overridden by namelist input:

        state = normal
        dotrace = false
        chkpoint = 0
        chkfile = 'ChkPt'
        exact = false
        expens = 1
        derv = analytic
        meth = available
        mterms = 0
        update = sumform
        set_h0 = ident

        run8 = true

        READ (inpt,nml=args,end=30) ! read namelist data

! Write title, unless no output requested.

        WRITE (trmout,'(''  Ready for run #'',i2,'':'', a)') m, title(m)
        IF (freq/=0) WRITE (outpt,90080) m, title(m)

! Here we have a chance to choose a subset of the problems to test.

        IF ( .NOT. run8) CYCLE ! If run8 = false, do not run test #8.

10      IF (contrl/=-3) THEN
          contrl = -3
          WRITE (trmout,*) ' control: 0 quit'
          WRITE (trmout,*) '         -1 skip to next set,'
          WRITE (trmout,*) '         -2 finish this set'
          WRITE (trmout,*) '         -3 (or eof) finish full run'
          WRITE (trmout,*) '        n > 0 do problem #n'
          READ (trminp,'(bn,i2)',end=20) contrl
        END IF

20      SELECT CASE (contrl)
        CASE (1:)
          from = contrl
          to = contrl
        CASE (0)
          GO TO 30
        CASE (-1)
          CYCLE
        CASE (-2)
          from = MOD(to,mprobs) + 1
          to = mprobs
        CASE (-3)
          from = MOD(to,mprobs) + 1
          to = mprobs
        END SELECT

! For tests 8 through 10 just use function 11 (pwsing with n=500).
! For tests 13 and 14 just use function 12 (sarseb).

        SELECT CASE (m)
        CASE (8,9,10)
          from = 11
          to = 11
          mprobs = 11
        CASE (11,12)
          IF (contrl<1) THEN
            from = 1
            to = 10
          END IF
          mprobs = 10
        CASE (13,14)
          from = nprobs
          to = nprobs
          mprobs = nprobs
        END SELECT

        time = cpusecs() ! Start timing
        acctim = acctim - time

EACHFUNCTION: DO i = from, to !Repeat for each test function selected.

          pdone = pdone + 1
          IF (pdone>ttests) THEN
            WRITE (outpt,*) ' too many tests: stopping.'
            GO TO 40
          END IF

          IF (freq/=0) WRITE (outpt,90030) i

          n = DIM(i) ! Set dimension, etc.
          funcno = ifnc(i)
          arguments(1) = rpar1(i)
          arguments(2) = rpar2(i)

          SELECT CASE (i) ! Set starting point.
          CASE (1)
            x(1:n) = (/ one, one, one/)
          CASE (2)
            x(1:n) = (/ one, two, one, one, one, one/)
          CASE (3)
            x(1:n) = (/ zero, ten, c20/)
          CASE (4)
            x(1:n) = (/ one, two, two, two/)
          CASE (5)
            x(1:n) = (/ one, two, zero/)
          CASE (6)
            x(1:n) = (/ (i,i=1,n) /)
          CASE (7)
            x(1:n) = half
          CASE (8)
            x(1:n) = (/ (three,-one,zero,one,i=1,n/4) /)
          CASE (9)
            x(1:n) = (/ -1.2_stnd, one/)
          CASE (10)
            x(1:n) = half
          CASE (11)
            x(1:n) = (/ (three,-one,zero,one,i=1,n/4) /)
          CASE (12)
            x(1:n) = (/ c20, c20, c20, c20/)
          END SELECT

! ---Set up calls to minimize.

          IF (dotrace) THEN
            dotraces = tracelist(6,0,true,true,true,true,true,true,true,true, &
              true,true,true)
          ELSE
            dotraces = tracelist(6,0,false,false,false,false,false,false, &
              false,false,false,false,false)
          END IF

          SELECT CASE (m)

          CASE (1,3,4,5,7:) ! Forward calls

            status = state
            CALL minimize_f(functions,x(1:n),fx,g,acc,status,memory,c, &
              traces=dotraces,checkpoint=chkpoint,checkfile=chkfile, &
              exactls=exact,expense=expens,relativetof0=true, &
              relativetog0=true,usegrad=false,usestep=true,printunit=8_short, &
              frequency=freq,derivatives=derv,method=meth,terms=mterms, &
              updateform=update,seth0=set_h0)

            IF (m==4) THEN
              dererr(i) = c%evalers%worst
              averrs(i) = c%evalers%average
              compnt(i) = c%evalers%index
              indx(i) = c%evalers%gradcnt
            END IF

          CASE (2,6) ! Reverse communication

            status = state
            c%eval = defaultevalstate
            firsttime = true

            DO
              case = both
              CALL evaluate_f(functions,x(1:n),fx,g(1:n),case,c%eval, &
                c%evalcts,first=firsttime)
              CALL minimize_f(functions,x(1:n),fx,g,acc,status,memory,c, &
                traces=dotraces,checkpoint=chkpoint,checkfile=chkfile, &
                exactls=exact,expense=expens,relativetof0=true, &
                relativetog0=true,usegrad=false,usestep=true, &
                printunit=8_short,frequency=freq,derivatives=derv,method=meth, &
                terms=mterms,updateform=update,seth0=set_h0)

              IF (status/=revcommfandg) EXIT
              status = revcommrestart
              firsttime = false
            END DO

          END SELECT

          IF (status/=done) THEN !  Add number of errors.
            error = error + 1
          END IF

          fnct = c%evalcts%fevals
          grct = c%evalcts%gevals

          fvals(pdone) = fx
          fcnts(pdone) = c%evalcts%fevals
          icnts(pdone) = c%shared%it

          tfncs = tfncs + c%evalcts%fevals
          titers = titers + c%shared%it

          IF (freq/=0) WRITE (outpt,90060) status

        END DO EACHFUNCTION

        time = cpusecs()
        acctim = acctim + time

        IF (m==4) THEN
          IF (freq/=0) WRITE (outpt,90070) (dererr(i),compnt(i),indx(i), &
            averrs(i),i=1,mprobs)
        END IF

        IF (to/=mprobs) GO TO 10

      END DO TESTSET

30    WRITE (outpt,90000)
      WRITE (outpt,90010) (i,icnts(i),fvals(i),fcnts(i),i=1,pdone)

      WRITE (outpt,90050) pdone, error, tfncs, titers
      WRITE (outpt,90040) acctim
      WRITE (trmout,*) ' '
      WRITE (trmout,*) ' '
      WRITE (trmout,*) ' '
      WRITE (trmout,90040) acctim

      CALL my_date_time(date)
      WRITE (outpt,90020) ' test ended at ', date
      WRITE (trmout,90020) ' test ended at ', date

! EXIT:

40    STOP

! FORMATS:

90000 FORMAT (//' rn its funct. value fns |',' rn its funct. value fns |', &
        ' rn its funct. value fns'/' ------------------------|', &
        '-------------------------|','------------------------')

90010 FORMAT ((2(I3,I4,1X,E12.6,I4,' |'),(I3,I4,1X,E12.6,I4)))

90020 FORMAT (4A//)

90030 FORMAT (/' function #',I2/)

90040 FORMAT (' time taken was ',F12.3,' seconds.')

90050 FORMAT (//' *****Test Finished****   problems done ',I3, &
        '; number of errors is ',I2,'.'/ &
        '                total function calls = ',I4,' total iterations = ', &
        I4//)

90060 FORMAT (/' ************* Run Complete, status = ',I3,'.'/)

90070 FORMAT (//' Testing mode derivative estimation errors'// &
        ' max error  component  iterate  av. decimals '// &
        10(E10.2,I7,7X,I5,F9.2/))

90080 FORMAT ('1 beginning run #',I2,':',A)

!END:

    END

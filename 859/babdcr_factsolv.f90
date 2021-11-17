     SUBROUTINE BABDCR_FACTSOLV( NRWBLK, NBLOKS, MATR_A, LFTBLK, &
                                 RGTBLK, VECT_B, INFO )

     USE PRECISION
!
!  ==== BABDCR package ================
!  == subroutine BABDCR_FACTSOLV ======
!  Giuseppe Romanazzi, Pierluigi Amodio
!  Dipartimento di Matematica
!  Universita' di Bari
!  October 20, 2005
!  ====================================
!
! Purpose
! =======
!
!  BABDCR_FACTSOLV  solves a linear system Ax=b with the babd matrix 
!
!         ( Ba                                              Bb    )
!         ( S(0)   R(1)                                           )
!         (        S(1)   R(2)                                    )
!         (                                                       )
!    A =  (                 .        .                            )
!         (                        .        .                     )
!         (                               .         .             )
!         (                                                       )
!         (                                S(NBLOKS-1)  R(NBLOKS) )
!
!  of dimension NRWBLK*(NBLOKS+1) (each block is of size NRWBLK by
!  NRWBLK). The algorithm uses the cyclic reduction algorithm to 
!  factorize the coefficient matrix. The solution is obtained in three
!  phases: reduction (by using the subroutine REDUCE), solution of
!  the 2 by 2 block system (by using the subroutines DGETRF and DGETRS)
!  and back-substitution (by using the subroutine SOLVE_BLOCK).
!
!  In input, the coefficient matrix is stored in the NRWBLK by NRWBLK
!  by NBLOKS*2 array MATR_A
!
!   MATR_A  = [S(0) R(1) S(1) R(2),...., S(NBLOKS-1) R(NBLOKS)]
!
!  and in the NRWBLK by NRWBLK arrays LFTBLK = [Ba] and RGTBLK = [Bb].
!  The right hand side is a vector of dimension NRWBLK*(NBLOKS+1) and
!  it is stored in the NRWBLK by NBLOKS+1 array VECT_B in the
!  following form
!
!     VECT_B  = [f(0), f(1), ...., f(NBLOKS)].
!
!  On exit, BABDCR_FACTSOLV gives the solution of the linear system in
!  the array VECT_B. The coefficient matrix is modified but it contains
!  only a part of the factorization and cannot be used anymore.
!
!
! Parameters
! ==========
!
! Input variables:
!   NRWBLK  integer,            the dimension of each block of the
!                               babd matrix, NRWBLK>0
!
!   NBLOKS  integer,            the number of blocks S(i) and R(i)
!                               in the babd matrix, NBLOKS>0
!
!   MATR_A  double precision,   NRWBLK by NRWBLK by NBLOKS*2 array,
!                               the blocks S(i) and R(i) of the babd
!                               matrix, which are saved as previously 
!                               described
!
!   LFTBLK  double precision,   NRWBLK by NRWBLK array,
!                               the block Ba of the input babd matrix
!
!   RGTBLK  double precision,   NRWBLK by NRWBLK array,
!                               the block Bb of the input babd matrix
!
!   VECT_B  double precision,   NRWBLK by NBLOKS+1 array,
!                               the r.h.s. of the linear system
! Output variables:
!   MATR_A  double precision,   NRWBLK by NRWBLK by NBLOKS*2 array,
!                               part of the factorization of the babd 
!                               matrix
!
!   LFTBLK  double precision,   NRWBLK by NRWBLK array,
!                               part of the factorization of the last 
!                               2 by 2 reduced matrix
!
!   RGTBLK  double precision,   NRWBLK by NRWBLK array,
!                               part of the factorization of the last 
!                               2 by 2 reduced matrix
!
!   VECT_B  double precision,   NRWBLK by NBLOKS+1 array,
!                               the solution of the linear system
!   INFO    integer,
!           = 0:  successful exit
!           > 0:  if INFO = I, the I-th computed factorization is
!                 performed on a matrix of rank<NRWBLK; therefore the
!                 babd matrix is singular
!
     INTEGER,  INTENT(IN)     :: NRWBLK, NBLOKS
     REAL(DP), INTENT(IN OUT) :: MATR_A( NRWBLK, NRWBLK, NBLOKS*2 ), &
                                 LFTBLK( NRWBLK, NRWBLK ), &
                                 RGTBLK( NRWBLK, NRWBLK ), &
                                 VECT_B( NRWBLK, 0:NBLOKS )
     INTEGER,  INTENT(OUT)    :: INFO
! Local variables:
     REAL(DP) :: MATR_TEMP( NRWBLK*2, NRWBLK*2 ), TC( NRWBLK*2 )
     INTEGER  :: H, HH, Z, ZZ, JUMP, JUMP2, INDEXP, K, I, &
                 PERM( NRWBLK*2, NBLOKS ), NSOLSTEP( 20 ), NSTEPS
! Lapack routine: DGETRF
! Used subroutines: REASSEMBLE, REDUCE, SOLVE_BLOCK

!
! initialize indeces
!
      INFO= 0
      INDEXP= 0
      JUMP2= 2
      NSTEPS= 0
      K= NBLOKS
!
!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! reduction phase !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!
!
      DO WHILE ( K.GT.1 )
!
! set the number of internal reduction steps and initialize the index
! of the main block to reduce
!
        NSTEPS= NSTEPS+1
        NSOLSTEP(NSTEPS)= K/2
        HH= JUMP2
!
!! internal reduction cycle
!
        DO I= 1, K/2-1
          INDEXP= INDEXP+1
          ZZ= HH+JUMP2
!
! reshape ( Rb  Sb ) as ( Rb ) and save in ( , ) columns by columns
!                       ( Sb )
!
          CALL REASSEMBLE( NRWBLK, MATR_A(1,1,HH) )
!
! reduce 2 by 3 block system  ( Sa  Rb     ) ( Xm ) = ( Vb )
!                             (     Sb  Rc ) ( Xn )   ( Vc )
!                                            ( Xs )
! to  Sa' * Xm + Rc' * Xs = Vc'
!
          CALL REDUCE( NRWBLK, MATR_A(1,1,HH-JUMP2+1), &
                       MATR_A(1,1,HH), MATR_A(1,1,ZZ), &
                       PERM(1,INDEXP), VECT_B(1,HH/2), &
                       VECT_B(1,ZZ/2), INFO )

          IF ( INFO.NE.0 ) THEN
            INFO= INDEXP
            RETURN
          ENDIF

          HH= ZZ+JUMP2
        ENDDO
!
!! last step of internal reduction cycle
!
        INDEXP= INDEXP+1
        ZZ= HH+JUMP2
        IF ( ZZ.GT.NBLOKS*2 ) THEN
          ZZ= NBLOKS*2
        ENDIF

        CALL REASSEMBLE( NRWBLK, MATR_A(1,1,HH) )

        CALL REDUCE( NRWBLK, MATR_A(1,1,HH-JUMP2+1), MATR_A(1,1,HH), &
                     MATR_A(1,1,ZZ), PERM(1,INDEXP), VECT_B(1,HH/2), &
                     VECT_B(1,ZZ/2), INFO )
                     
        IF ( INFO.NE.0 ) THEN
          INFO= INDEXP
          RETURN
        ENDIF
!
!! end internal reduction cycle
!
        JUMP2= JUMP2*2
        K= K-K/2
      ENDDO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 2 by 2 block linear system solution !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     ( Ba    Bb    ) (x(0)     )= ( f(0)    )
!     ( S^(1) R^(1) ) (x(NBLOKS))  ( f^(1) )
!
      MATR_TEMP( 1:NRWBLK, 1:NRWBLK )= LFTBLK(:,:)
      MATR_TEMP( 1:NRWBLK, NRWBLK+1:NRWBLK*2 )= RGTBLK(:,:)
      MATR_TEMP( NRWBLK+1:NRWBLK*2, 1:NRWBLK )= MATR_A(:,:,1)
      MATR_TEMP( NRWBLK+1:NRWBLK*2, NRWBLK+1:NRWBLK*2 )= &
           MATR_A(:,:,NBLOKS*2)
!
! factorization
!
      CALL DGETRF( NRWBLK*2, NRWBLK*2, MATR_TEMP(1,1), NRWBLK*2, &
                   PERM(1,NBLOKS), INFO )

      IF ( INFO.NE.0 ) THEN
        INFO= INDEXP
        RETURN
      ENDIF
!
! solution
!
      TC= (/ VECT_B(:,0), VECT_B(:,NBLOKS) /)
      CALL DGETRS( 'N', NRWBLK*2, 1, MATR_TEMP(1,1), NRWBLK*2, &
                   PERM(1,NBLOKS), TC(1), NRWBLK*2, INFO )
      VECT_B(:,0)= TC(1:NRWBLK)
      VECT_B(:,NBLOKS)= TC(NRWBLK+1:NRWBLK*2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! back-substitution phase !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      JUMP= JUMP2/2
      INDEXP= NBLOKS-1
      DO K= NSTEPS,1,-1
!
!! first internal back-substitution step
!
        Z= JUMP*NSOLSTEP(K)
        JUMP= JUMP/2
        H= Z-JUMP
        IF ( Z.GT.NBLOKS ) THEN
          Z= NBLOKS
        ENDIF

        CALL SOLVE_BLOCK( NRWBLK, VECT_B(1,H-JUMP), VECT_B(1,Z), &
                          PERM(1,INDEXP), MATR_A(1,1,H*2), &
                          VECT_B(1,H) )
                          
        INDEXP= INDEXP-1
!
!! internal back-substitution cycle
!        
        DO I= 2, NSOLSTEP(K)
          H= H-JUMP*2
!
! determine Xn from the factorized 2 by 3 block system
!
!    ( Sa Rb    ) ( Xm )   ( Tb )
!    (    Sb Rc ) ( Xn ) = ( Tc )
!                 ( Xs )
!
! where Xm and Xs are known quantities.
!
          CALL SOLVE_BLOCK( NRWBLK, VECT_B(1,H-JUMP), &
                            VECT_B(1,H+JUMP), PERM(1,INDEXP), &
                            MATR_A(1,1,H*2), VECT_B(1,H) )
                            
          INDEXP= INDEXP-1
        ENDDO
      ENDDO

     RETURN
      
     CONTAINS
      
       SUBROUTINE REDUCE( NRWBLK, S, RS, R, P, Vb, Vc, INFO )
!
! Purpose
! =======
!
!  REDUCE  reduces a NRWBLK*2 by NRWBLK*3 block system
!
!   ( Sa  Rb     ) ( Xm )    ( Vb )
!   (     Sb  Rc ) ( Xn ) =  ( Vc )
!                  ( Xs )
!
!  to a NRWBLK by NRWBLK*2 system
!
!     ( Sa', Rc') ( Xm ) = Vc'
!                 ( Xs )
!
!  by using the following algorithm:
!
!    P * ( Rb ) = ( L ) * U;
!        ( Sb )   ( C )
!
!    P * ( Sa     ) = ( Sp  Rp ),  P * ( Vb ) = ( Vp )
!        (     Rc )   ( Sq  Rq )       ( Vc )   ( Vq )
!  
!    ( I    ) * P * ( Sa  Rb     ) = ( Sp  L*U  Rp )
!    ( G  I )       (     Sb  Rc )   ( Sa'      Rc' )
!  
!    ( I    ) * P * ( Vb ) = ( Vp )
!    ( G  I )       ( Vc )   ( V' )
!  
!  where  Sa' = Sq - G*Sp,  Rc' = Rq - G*Rp,  V' = Vq - G*Vp 
!  and  G = C*L^(-1). Each block is NRWBLK by NRWBLK.
!
!  On exit, REDUCE gives the computed blocks Sa' and Rc', the L\U
!  matrices of the LU factorization, the nonnull rows of Sp and Rp,
!  and the block G. Moreover, it gives the computed vector V' and Vp.
!
!
! Parameters
! ==========
!
! Input variables:
!   NRWBLK  integer,            the dimension of each block of the
!                               babd matrix, NRWBLK>0
!
!   S       double precision,   NRWBLK by NRWBLK array,
!                               the block Sa
!
!   RS      double precision,   NRWBLK*2 by NRWBLK array,
!                               the block ( Rb ) 
!                                         ( Sb )
!
!   R       double precision,   NRWBLK by NRWBLK array,
!                               the block Rc
!
!   Vb      double precision,   NRWBLK vector,
!                               the first block vector of the r.h.s to 
!                               reduce
!
!   Vc      double precision,   NRWBLK vector,
!                               the second block vector of the r.h.s to 
!                               reduce
! Output variables:
!   RS      double precision,   NRWBLK*2 by NRWBLK array,
!                               the L\U factors of the LU factorization
!                               of RS in the first NRWBLK rows; the
!                               transpose of the matrix
!                                    P(1:NRWBLK,:)*( Sa )
!                                                  ( Rc ),
!                               where P(1:NRWBLK,:) contains the first
!                               NRWBLK rows of the permutation matrix
!                               associated to the factorization of RS
!
!   S       double precision,   NRWBLK by NRWBLK array,
!                               the block Sa' obtained by the reduction
!
!   R       double precision,   NRWBLK by NRWBLK array,
!                               the block Rc' obtained by the reduction
!
!   P       integer,            NRWBLK*2 vector,
!                               permutation vector associated to the
!                               factorization of RS, 1<=P(I)<=NRWBLK*2
!
!   Vb      double precision,   NRWBLK vector,
!                               the vector Vp obtained by the
!                               permutation
!
!   Vc      double precision,   NRWBLK  vector,
!                               the vector V' obtained by the reduction
!
!   INFO    integer,
!           = 0:  successful exit
!           < 0:  if INFO = -I, the I-th argument of DGETRF calling has 
!                an illegal value
!           > 0:  if INFO = I, U(I,I) of factorization is exactly zero. 
!                 The factorization has been completed, but the factor 
!                 U is singular, and a division by zero occurs when the
!                 system is solved
!
       INTEGER,  INTENT (IN)    :: NRWBLK
       REAL(DP), INTENT(IN OUT) :: RS( NRWBLK*2, NRWBLK ), &
                                    S( NRWBLK, NRWBLK ), &
                                    R( NRWBLK, NRWBLK ), &
                                    Vb( NRWBLK ), Vc( NRWBLK )
       INTEGER,  INTENT(OUT)    :: P( NRWBLK*2 ), INFO
! Parameters:
       REAL(DP), PARAMETER :: ONE= 1.0D0, NONE= -1.0D0
! Local variables:
       INTEGER  :: I, J, T, P1( NRWBLK )
       REAL(DP) :: H( NRWBLK, NRWBLK ), L( NRWBLK, NRWBLK ), &
                   TEMP( NRWBLK*2 ), F( NRWBLK, NRWBLK )
! Lapack routine: DGETRF
! Blas routines: DTRSM, DGER, DGEMV

!
! factorize RS by means of the LU factorization
!
        CALL DGETRF( NRWBLK*2, NRWBLK, RS(1,1), NRWBLK*2, P1(1), INFO )
        IF ( INFO.NE.0 ) THEN
          RETURN
        ENDIF
!
! compute the block G=C*L^(-1) and save it in RS(NRWBLK+1:NRWBLK*2,:)
!
        CALL DTRSM( 'R', 'L', 'N', 'U', NRWBLK, NRWBLK, ONE, RS(1,1), &
                    NRWBLK*2, RS(NRWBLK+1,1), NRWBLK*2 )
!
! determine the permutation vector
!
        DO I= 1, NRWBLK*2
          P(I)= I
        ENDDO
        DO I= 1, NRWBLK
          J= P1(I)
          IF ( J.NE.I ) THEN
            T= P(I)
            P(I)= P(J)
            P(J)= T
          ENDIF
        ENDDO
!
! compute Sa' = Sq-G*Sp,  Rc' = Rq-G*Rp and  T'= Tq-G*Tp
! in two steps:
! 1. define H and L as S and R, respectivly; S=Sq, R=Rq
!
        H= S
        L= R
        DO I= 1, NRWBLK
          J= P(I+NRWBLK)
          IF ( J.LE.NRWBLK ) THEN
            S(I,:)= H(J,:)
            R(I,:)= 0
          ELSE
            S(I,:)= 0
            R(I,:)= L(J-NRWBLK,:)
          ENDIF
        ENDDO
!
! 2. compute  S = S-G*Sp and R = R-G*Rp, where the nonnull rows 
!    of Sp and Rp are saved in F'
!
        DO I= 1, NRWBLK
          IF ( P(I).LE.NRWBLK ) THEN
! save nonnull row of Sp in the column of F
            F(:,I)= H(P(I),:)
            CALL DGER( NRWBLK, NRWBLK, NONE, RS(NRWBLK+1,I), &
                       1, F(1,I), 1, S, NRWBLK )
          ELSE
! save nonnull row of Rp in the column of F
            F(:,I)= L(P(I)-NRWBLK,:)
            CALL DGER( NRWBLK, NRWBLK, NONE, RS(NRWBLK+1,I), &
                       1, F(1,I), 1, R, NRWBLK )
          ENDIF
        ENDDO
        
        TEMP= (/ Vb(:), Vc(:) /)
        Vb= TEMP(P(1:NRWBLK))
        Vc= TEMP(P(NRWBLK+1:NRWBLK*2))
        CALL DGEMV( 'N', NRWBLK, NRWBLK, NONE, RS(NRWBLK+1,1), &
                    NRWBLK*2, Vb, 1, ONE, Vc, 1 )
        RS(NRWBLK+1:NRWBLK*2,:)= F

       RETURN
       END SUBROUTINE REDUCE


       SUBROUTINE  REASSEMBLE( NRWBLK, MATR )
!
! Purpose
! =======
!
!  REASSEMBLE  assembles a 1 by 2 block array MATR as a 2 by 1 block
!  matrix, where each block is of dimension NRWBLK x NRWBLK
!
!  On exit, REASSEMBLE gives the permutation of the block vector MATR
!
!
! Parameters
! ==========
!
! Input variables:
!   NRWBLK  integer,            the dimension of each block of the
!                               babd matrix, NRWBLK>0
!
!   MATR    double precision,   NRWBLK by NRWBLK*2 array
!
! Output variables:
!   MATR    double precision,   NRWBLK by NRWBLK*2 array,
!                               permutation of the input matrix from
!                               ( R  S ) to ( R )
!                                           ( S )
!
       INTEGER,  INTENT(IN)     :: NRWBLK
       REAL(DP), INTENT(IN OUT) :: MATR( NRWBLK, NRWBLK*2 )
! Local variables:
       INTEGER  :: I
       REAL(DP) :: TEMP( NRWBLK, NRWBLK )

        DO I = 2, NRWBLK
          TEMP(:,I) = MATR(:,I)
        ENDDO
        DO I = 1, NRWBLK-1
          MATR(:,I*2) = MATR(:,NRWBLK+I)
        ENDDO
        DO I = 2, NRWBLK
          MATR(:,I*2-1) = TEMP(:,I)
        ENDDO

       RETURN
       END SUBROUTINE REASSEMBLE


       SUBROUTINE SOLVE_BLOCK( NRWBLK, Xm, Xs, P, RS, V )
!
! Purpose
! =======
!
!  SOLVE_BLOCK  computes the solution Xn (of length NRWBLK) of the 
!  linear system
!
!          L*U Xn = Vp - Sp*Xm - Rp*Xs
!
!  obtained from the subroutines REDUCE_BLOCK and REDUCE_RHS applied
!  to the system (see the subroutine REDUCE_BLOCK for further details)
!
!   ( Sa  Rb     )  ( Xm )   ( Vb ).
!   (     Sb  Rc )  ( Xn ) = ( Vc )
!                   ( Xs )
!
!  Each block is NRWBLK by NRWBLK. The r.h.s. Vp is obtained after a
!  row permutation of Va and Vb (see the subroutine SOLVE_BLOCK).
!
!  On exit, SOLVE_BLOCK gives the solution Xn in the array V.
!
!
! Parameters
! ==========
!
! Input variables:
!   NRWBLK  integer,            the dimension of each block of the
!                               babd matrix, NRWBLK>0
!
!   Xm      double precision,   NRWBLK vector, 
!                               the previously obtained solution Xm
!
!   Xs      double precision,   NRWBLK vector,
!                               the previously obtained solution Xs
!
!   P       integer,            NRWBLK*2 vector,
!                               permutation vector that comes out from 
!                               REDUCE_BLOCK
!
!   RS      double precision,   NRWBLK*2 by NRWBLK array,
!                               RS array that comes out from
!                               REDUCE_BLOCK
!
!   V       double precision,   NRWBLK vector, 
!                               the r.h.s. Vp
!
! Output variables:
!   V       double precision,   NRWBLK vector, 
!                               the solution Xn
!
       INTEGER,  INTENT (IN)     ::  NRWBLK, P( NRWBLK*2 )
       REAL(DP), INTENT (IN)     ::  Xm( NRWBLK ), Xs( NRWBLK ), &
                                     RS( NRWBLK*2, NRWBLK )
       REAL(DP), INTENT (IN OUT) ::  V( NRWBLK )
! Local variables:
       INTEGER :: I, INFO
       INTEGER :: P1( NRWBLK )
       REAL(DP) :: DDOT
! Blas routine: DDOT
! Lapack routine: DGETRS
! my routine: DGETRS_MOD

!
! compute V = Vp - Sp*Xm - Rp*Xs
!
        DO I= 1, NRWBLK
          IF ( P(I).LE.NRWBLK ) THEN
            V(I)= V(I) - DDOT( NRWBLK, RS(NRWBLK+1,I), 1, Xm(1), 1)
          ELSE
            V(I)= V(I) - DDOT( NRWBLK, RS(NRWBLK+1,I), 1, Xs(1), 1)
          ENDIF 
        ENDDO
!
! solve the system L*U Xn = V
!
        DO I= 1, NRWBLK
          P1( I )= I
        ENDDO
        CALL DGETRS( 'N', NRWBLK, 1, RS(1,1), NRWBLK*2, P1(1), V(1), &
                     NRWBLK, INFO )
! the solution of the above system may be improved by using an
! appropriate subroutine solving linear systems without any
! permutation
!        CALL DGETRS_MOD( 'N', NRWBLK, 1, RS(1,1), NRWBLK*2, V(1), &
!     &                   NRWBLK, INFO )

       RETURN
       END SUBROUTINE SOLVE_BLOCK     

     END SUBROUTINE BABDCR_FACTSOLV

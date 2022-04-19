     SUBROUTINE BABD_MATMULT( TRANS, NRWBLK, NBLOKS, MATR_A, LFTBLK, &
                              RGTBLK, VECT_B, SOL_Y )

     USE PRECISION
!
!  == subroutine BABD_MATMULT =========
!  Giuseppe Romanazzi, Pierluigi Amodio
!  Dipartimento di Matematica
!  Universita' di Bari
!  October 20, 2005
!  ====================================
!
! Purpose
! =======
!
!  BABD_MATMULT  computes a matrix vector product  Y = A * B  or 
!  Y = A^T * B  with the babd matrix 
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
!  NRWBLK).
!
!  In input, A is stored in the NRWBLK by NRWBLK by 2 by
!  NBLOKS array MATR_A
!
!   MATR_A  = [S(0) R(1) S(1) R(2),...., S(NBLOKS-1) R(NBLOCKS)]
!
!  and in the NRWBLK by NRWBLK arrays LFTBLK = [Ba] and RGTBLK = [Bb].
!  The array B has dimension NRWBLK*(NBLOKS+1) and is stored in the
!  NRWBLK by NBLOKS+1 array VECT_B
!
!  On exit, BABD_MATMULT the solution Y has dimension NRWBLK*(NBLOKS+1) 
!  and is stored in the NRWBLK by NBLOKS+1 array SOL_Y
!
!
! Parameters
! ==========
!
! Input variables:
!   TRANS   character,          the type of product (A or A^T)
!
!   NRWBLK  integer,            the size of blocks S, R, Ba, Bb,
!                               NRWBLK>0
!
!   NBLOKS  integer,            the number of blocks S(i) and R(i)
!                               in the babd matrix, NBLOKS>0
!
!   MATR_A  double precision,   NRWBLK by NRWBLK by 2 by NBLOKS array,
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
!                               the r.h.s. of the product
!
! Output variables:
!   SOL_Y double precision,     NRWBLK by NBLOKS+1 array,
!                               the result of the matrix-vector 
!                               multiplication
!
     CHARACTER, INTENT(IN) :: TRANS 
     INTEGER, INTENT(IN)   :: NRWBLK, NBLOKS
     REAL(DP), INTENT(IN)  :: MATR_A(NRWBLK,NRWBLK,2,NBLOKS), &
                              LFTBLK(NRWBLK,NRWBLK), &
                              RGTBLK(NRWBLK,NRWBLK), &
                              VECT_B(NRWBLK,0:NBLOKS)
     REAL(DP), INTENT(OUT) :: SOL_Y(NRWBLK,0:NBLOKS)
! Local variables:
     INTEGER :: I
     LOGICAL :: LSAME
! Lapack routine: DGEMV, LSAME
  
      SOL_Y= 0
!
! Y = A * B
!

      IF( LSAME( TRANS, 'N' ) ) THEN
        CALL DGEMV( 'N', NRWBLK, NRWBLK, 1d0, LFTBLK(1,1), NRWBLK, &
                    VECT_B(:,0), 1, 0d0, SOL_Y(:,0), 1)
        CALL DGEMV( 'N', NRWBLK, NRWBLK, 1d0, RGTBLK(1,1), NRWBLK, &
                    VECT_B(:,NBLOKS), 1, 1d0, SOL_Y(:,0), 1)
        DO I= 1,NBLOKS
          CALL DGEMV( 'N', NRWBLK, NRWBLK, 1d0, MATR_A(1,1,1,I), &
                      NRWBLK, VECT_B(:,I-1), 1, 0d0, SOL_Y(:,I), 1)
          CALL DGEMV( 'N', NRWBLK, NRWBLK, 1d0, MATR_A(1,1,2,I), &
                      NRWBLK, VECT_B(:,I), 1, 1d0, SOL_Y(:,I), 1)
        ENDDO
      ELSE
!
! Y = A^T * B
!
        CALL DGEMV( 'T', NRWBLK, NRWBLK, 1d0, LFTBLK(1,1), NRWBLK, &
                    VECT_B(:,0), 1, 0d0, SOL_Y(:,0), 1)
        DO I= 1,NBLOKS
          CALL DGEMV( 'T', NRWBLK, NRWBLK, 1d0, MATR_A(1,1,1,I), &
                      NRWBLK, VECT_B(:,I-1), 1, 1d0, SOL_Y(:,I-1), 1)
          CALL DGEMV( 'T', NRWBLK, NRWBLK, 1d0, MATR_A(1,1,2,I), &
                      NRWBLK, VECT_B(:,I), 1, 0d0, SOL_Y(:,I), 1)
        ENDDO
        CALL DGEMV( 'T', NRWBLK, NRWBLK, 1d0, RGTBLK(1,1), NRWBLK, &
                    VECT_B(:,NBLOKS), 1, 1d0, SOL_Y(:,NBLOKS), 1)
      ENDIF

     RETURN
     END SUBROUTINE BABD_MATMULT

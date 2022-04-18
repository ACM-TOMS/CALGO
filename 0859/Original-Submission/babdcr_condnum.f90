     REAL(DP) FUNCTION BABDCR_NORM1( NRWBLK, NBLOKS, MATR_A, &
                                     LFTBLK, RGTBLK )
     USE PRECISION
!
!  == subroutine BABD_NORM1 ===========
!  Giuseppe Romanazzi, Pierluigi Amodio
!  Dipartimento di Matematica
!  Universita' di Bari
!  October 20, 2005
!  ====================================
!
! Purpose
! =======
!
!  BABD_MATMULT  computes the 1-norm of the babd matrix 
!
!         ( Ba                                              Bb    )
!         ( S(0)   R(1)                                           )
!         (        S(1)   R(2)                                    )
!         (                                                       )
!         (                 .        .                            )
!         (                        .        .                     )
!         (                               .         .             )
!         (                                                       )
!         (                                S(NBLOKS-1)  R(NBLOKS) )
!
!  of dimension NRWBLK*(NBLOKS+1) (each block is of size NRWBLK by 
!  NRWBLK).
!
     INTEGER,  INTENT(IN)     :: NRWBLK, NBLOKS
     REAL(DP), INTENT(IN)     :: MATR_A( NRWBLK, NRWBLK, 2, NBLOKS ), &
                                 LFTBLK( NRWBLK, NRWBLK ), &
                                 RGTBLK( NRWBLK, NRWBLK )
!Local variables
     INTEGER  :: I, J, K
     REAL(DP) :: X( NRWBLK, 0:NBLOKS ), EST

!
! compute the 1-norm of A , ||A||_1=max_{j=1,n}(\sum_{i=1}^{n} |aij| )
!     
      X(:,0)= DABS( LFTBLK(1,:) )
      DO J= 2, NRWBLK
        X(:,0)= X(:,0)+DABS( LFTBLK(J,:) )
      ENDDO
      DO I=1,NBLOKS
        DO J= 1, NRWBLK
          DO K= 1, NRWBLK
            X(K,I-1)= X(K,I-1)+DABS( MATR_A(J,K,1,I) )
          ENDDO
        ENDDO
        DO K=1, NRWBLK
          X(K,I)= DABS( MATR_A(1,K,2,I) )
        ENDDO
        DO J= 1, NRWBLK
          DO K= 1, NRWBLK
            X(K,I)= X(K,I)+DABS( MATR_A(J,K,2,I) )
          ENDDO
        ENDDO
      ENDDO
      DO J= 1, NRWBLK
        X(:,NBLOKS)= X(:,NBLOKS)+DABS( RGTBLK(J,:) )
      ENDDO
      BABDCR_NORM1= MAXVAL( X(:,0) )
      DO I= 1, NBLOKS
        EST= MAXVAL( X(:,I) )
        BABDCR_NORM1 = MAX( BABDCR_NORM1, EST )
      ENDDO

     RETURN
     END FUNCTION BABDCR_NORM1


     REAL(DP) FUNCTION BABDCR_NORM1INV( NRWBLK, NBLOKS, MATR_A, &
                              LFTBLK, RGTBLK, PERM, FILL_IN )

     USE BABDCR
!
!  == subroutine BABD_NORM1INV ========
!  Giuseppe Romanazzi, Pierluigi Amodio
!  Dipartimento di Matematica
!  Universita' di Bari
!  October 20, 2005
!  ====================================
!
! Purpose
! =======
!
!  BABD_MATMULT  computes an estimate of the 1-norm of the inverse 
!  of a babd matrix 
!
!         ( Ba                                              Bb    )
!         ( S(0)   R(1)                                           )
!         (        S(1)   R(2)                                    )
!         (                                                       )
!         (                 .        .                            )
!         (                        .        .                     )
!         (                               .         .             )
!         (                                                       )
!         (                                S(NBLOKS-1)  R(NBLOKS) )
!
!  of dimension NRWBLK*(NBLOKS+1) (each block is of size NRWBLK by 
!  NRWBLK).
!
     INTEGER,  INTENT(IN)     :: NRWBLK, NBLOKS, &
                                 PERM( NRWBLK*2, NBLOKS )
     REAL(DP), INTENT(IN)     :: MATR_A( NRWBLK, NRWBLK, 2, NBLOKS ), &
                                 LFTBLK( NRWBLK, NRWBLK ), &
                                 RGTBLK( NRWBLK, NRWBLK ), &
                                 FILL_IN( NRWBLK, NRWBLK, NBLOKS-1 )
!Local variables
     INTEGER  :: KASE, ISIGN( NRWBLK*(NBLOKS+1) ), ISOLVE, N
     REAL(DP) :: X( NRWBLK*(NBLOKS+1) ), V( NRWBLK*(NBLOKS+1) ), EST

! initialize variables
      ISOLVE= 0
      KASE= 0
      N = NRWBLK*(NBLOKS+1)
      X= 0.0D0

      CALL DONEST( N, V, X, ISIGN, EST, KASE )
      DO WHILE (KASE .NE. 0)
        ISOLVE = ISOLVE+1
        IF (KASE .EQ. 1) THEN
!overwrite X by (A^-1)*X
          CALL BABDCR_SOLV( NRWBLK, NBLOKS, MATR_A, LFTBLK, RGTBLK, &
                            PERM, FILL_IN, X )
        ELSE
!overwrite X by (A^-t)*X
          CALL BABDCR_SOLVT( NRWBLK, NBLOKS, MATR_A, LFTBLK, RGTBLK, &
                             PERM, FILL_IN, X )
        END IF
        CALL DONEST( N, V, X, ISIGN, EST, KASE )
      END DO
      BABDCR_NORM1INV= EST

     RETURN
     END FUNCTION BABDCR_NORM1INV

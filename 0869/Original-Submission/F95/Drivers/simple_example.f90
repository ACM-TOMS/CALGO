PROGRAM ODRPACK95_EXAMPLE
   USE ODRPACK95
   USE REAL_PRECISION
   REAL (KIND=R8), ALLOCATABLE :: BETA(:),L(:),U(:),X(:,:),Y(:,:)
   INTEGER :: NP,N,M,NQ
   INTERFACE
      SUBROUTINE FCN(N,M,NP,NQ,LDN,LDM,LDNP,BETA,XPLUSD,IFIXB,IFIXX,LDIFX,&
         IDEVAL,F,FJACB,FJACD,ISTOP)
         USE REAL_PRECISION
         INTEGER :: IDEVAL,ISTOP,LDIFX,LDM,LDN,LDNP,M,N,NP,NQ
         REAL (KIND=R8) :: BETA(NP),F(LDN,NQ),FJACB(LDN,LDNP,NQ), &
            FJACD(LDN,LDM,NQ),XPLUSD(LDN,M)
         INTEGER :: IFIXB(NP),IFIXX(LDIFX,M)
      END SUBROUTINE FCN
   END INTERFACE

   NP = 2
   N  = 4
   M  = 1
   NQ = 1
   ALLOCATE(BETA(NP),L(NP),U(NP),X(N,M),Y(N,NQ))
   BETA(1:2) = (/ 2.0_R8, 0.5_R8 /)
   L(1:2)    = (/ 0.0_R8, 0.0_R8 /)
   U(1:2)    = (/ 10.0_R8, 0.9_R8 /)
   X(1:4,1)  = (/ 0.982_R8, 1.998_R8, 4.978_R8, 6.01_R8 /)
   Y(1:4,1)  = (/ 2.7_R8, 7.4_R8, 148.0_R8, 403.0_R8 /)
   CALL ODR(FCN,N,M,NP,NQ,BETA,Y,X,LOWER=L,UPPER=U)
END PROGRAM ODRPACK95_EXAMPLE


SUBROUTINE FCN(N,M,NP,NQ,LDN,LDM,LDNP,BETA,XPLUSD,IFIXB,IFIXX,LDIFX,&
   IDEVAL,F,FJACB,FJACD,ISTOP)

   USE REAL_PRECISION

   INTEGER :: IDEVAL,ISTOP,LDIFX,LDM,LDN,LDNP,M,N,NP,NQ, I
   REAL (KIND=R8) :: BETA(NP),F(LDN,NQ),FJACB(LDN,LDNP,NQ), &
      FJACD(LDN,LDM,NQ),XPLUSD(LDN,M)
   INTEGER :: IFIXB(NP),IFIXX(LDIFX,M)

   ISTOP = 0

   ! Calculate model.
   IF (MOD(IDEVAL,10).NE.0) THEN
      DO I=1,N
         F(I,1) = BETA(1)*EXP(BETA(2)*XPLUSD(I,1))
      END DO
   END IF

   ! Calculate model partials with respect to BETA.
   IF (MOD(IDEVAL/10,10).NE.0) THEN
      DO I=1,N
         FJACB(I,1,1) = EXP(BETA(2)*XPLUSD(I,1))
         FJACB(I,2,1) = BETA(1)*XPLUSD(I,1)*EXP(BETA(2)*XPLUSD(I,1))
      END DO
   END IF

   ! Calculate model partials with respect to DELTA.
   IF (MOD(IDEVAL/100,10).NE.0) THEN
      DO I=1,N
         FJACD(I,1,1) = BETA(1)*BETA(2)*EXP(BETA(2)*XPLUSD(I,1))
      END DO
   END IF
 

END SUBROUTINE FCN
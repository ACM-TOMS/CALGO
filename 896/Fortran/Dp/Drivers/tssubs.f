! SUBROUTINE TIUB14                ALL SYSTEMS                99/12/01
! PORTABILITY : ALL SYSTEMS
! 92/12/01 RA : ORIGINAL VERSION
!
! PURPOSE :
!  INITIAL VALUES OF THE VARIABLES AND STRUCTURE OF THE SPARSE HESSIAN
!  MATRIX FOR UNCONSTRAINED MINIMIZATION.
!  SPARSE VERSION WITH CHANGED TESTS 7-10, 12.
!  CHANGED FOR THE TESTS.
!
! PARAMETERS :
!  II  N  NUMBER OF VARIABLES.
!  II  NB  NUMBER OF ELEMENTS OF THE SPARSE MATRIX.
!  RO  X(N)  VECTOR OF VARIABLES.
!  IO  IH(N+1)  POINTERS OF THE DIAGONAL ELEMENTS OF THE HESSIAN MATRIX.
!  IO  JH(M)  INDICES OF THE NONZERO ELEMENTS OF THE HESSIAN MATRIX IN
!             THE PACKED ROW.
!  II  NEXT  NUMBER OF THE TEST PROBLEM.
!  IO  IERR  ERROR INDICATOR.
!
      SUBROUTINE TIUB14 (N, NB, MB, X, IH, JH, FMIN, XMAX, NEXT, IERR)
      INTEGER N,NB,MB,NEXT,IERR
      INTEGER IH(*),JH(*)
      DOUBLE PRECISION X(N),FMIN,XMAX
      DOUBLE PRECISION P,Q
      INTEGER I,J,K
      DOUBLE PRECISION ETA9
      PARAMETER  (ETA9=1.0D60)
      FMIN=0.0D0
      XMAX=1.0D3
      IERR=0
      GO TO (10,50,90,110,130,170,200,240,280,300,320,350,370,390,410,
     &430,450,470,490,510,530,550),NEXT
   10 IF (N.LT.2) GO TO 570
      N=N-MOD(N,2)
      DO 20 I=1,N
        IF (MOD(I,2).EQ.1) THEN
          X(I)=-1.2D0
        ELSE
          X(I)=1.0D0
        END IF
   20 CONTINUE
   30 NB=N-1
      DO 40 I=1,NB
        J=2*(I-1)+1
        JH(J)=I
        JH(J+1)=I+1
        IH(I)=J
   40 CONTINUE
      MB=2*NB
      IH(NB+1)=MB+1
      RETURN
   50 IF (N.LT.4) GO TO 570
      N=N-MOD(N,2)
      DO 60 I=1,N
        IF (MOD(I,2).EQ.1) THEN
          X(I)=-2.0D0
          IF (I.LE.4) X(I)=-3.0D0
        ELSE
          X(I)=0.0D0
          IF (I.LE.4) X(I)=-1.0D0
        END IF
   60 CONTINUE
   70 NB=(N-2)/2
      DO 80 I=1,NB
        J=4*(I-1)+1
        K=2*I-1
        JH(J)=K
        JH(J+1)=K+1
        JH(J+2)=K+2
        JH(J+3)=K+3
        IH(I)=J
   80 CONTINUE
      MB=4*NB
      IH(NB+1)=IH(NB)+4
      RETURN
   90 IF (N.LT.4) GO TO 570
      N=N-MOD(N,2)
      DO 100 I=1,N
        IF (MOD(I,4).EQ.1) THEN
          X(I)=3.0D0
        ELSE IF (MOD(I,4).EQ.2) THEN
          X(I)=-1.0D0
        ELSE IF (MOD(I,4).EQ.3) THEN
          X(I)=0.0D0
        ELSE
          X(I)=1.0D0
        END IF
  100 CONTINUE
      GO TO 70
  110 IF (N.LT.4) GO TO 570
      N=N-MOD(N,2)
      DO 120 I=1,N
        X(I)=2.0D0
  120 CONTINUE
      X(1)=1.0D0
!      XMAX=1.0D 1
      GO TO 70
  130 IF (N.LT.3) GO TO 570
      DO 140 I=1,N
        X(I)=-1.0D0
  140 CONTINUE
  150 MB=1
      NB=N
      DO 160 I=1,NB
        IH(I)=MB
        IF (I.GT.1) THEN
          JH(MB)=I-1
          MB=MB+1
        END IF
        JH(MB)=I
        MB=MB+1
        IF (I.LT.N) THEN
          JH(MB)=I+1
          MB=MB+1
        END IF
  160 CONTINUE
      IH(NB+1)=MB
      MB=MB-1
      RETURN
  170 IF (N.LT.7) GO TO 570
      DO 180 I=1,N
        X(I)=-1.0D0
  180 CONTINUE
      MB=1
      NB=N
      DO 190 I=1,NB
        IH(I)=MB
        IF (I.GT.5) THEN
          JH(MB)=I-5
          MB=MB+1
        END IF
        IF (I.GT.4) THEN
          JH(MB)=I-4
          MB=MB+1
        END IF
        IF (I.GT.3) THEN
          JH(MB)=I-3
          MB=MB+1
        END IF
        IF (I.GT.2) THEN
          JH(MB)=I-2
          MB=MB+1
        END IF
        IF (I.GT.1) THEN
          JH(MB)=I-1
          MB=MB+1
        END IF
        JH(MB)=I
        MB=MB+1
        IF (I.LT.N) THEN
          JH(MB)=I+1
          MB=MB+1
        END IF
  190 CONTINUE
      IH(NB+1)=MB
      MB=MB-1
      RETURN
  200 IF (N.LT.4) GO TO 570
      N=N-MOD(N,2)
      DO 210 I=1,N
        X(I)=-1.0D0
  210 CONTINUE
      MB=1
      K=N/2
      NB=N+K
      DO 220 I=1,N
        IH(I)=MB
        IF (I.GT.1) THEN
          JH(MB)=I-1
          MB=MB+1
        END IF
        JH(MB)=I
        MB=MB+1
        IF (I.LT.N) THEN
          JH(MB)=I+1
          MB=MB+1
        END IF
  220 CONTINUE
      DO 230 I=1,K
        IH(N+I)=MB
        JH(MB)=I
        MB=MB+1
        JH(MB)=I+K
        MB=MB+1
  230 CONTINUE
      IH(NB+1)=MB
      MB=MB-1
      RETURN
  240 IF (N.LT.6) GO TO 570
      DO 250 I=1,N
        X(I)=1.0D0/DBLE(N)
  250 CONTINUE
  260 MB=1
      NB=N
      K=N/2
      DO 270 I=1,NB
        IH(I)=MB
        IF (I.GT.K) THEN
          JH(MB)=I-K
          MB=MB+1
        END IF
        IF (I.GT.2) THEN
          JH(MB)=I-2
          MB=MB+1
        END IF
        IF (I.GT.1) THEN
          JH(MB)=I-1
          MB=MB+1
        END IF
        JH(MB)=I
        MB=MB+1
        IF (I.LT.N) THEN
          JH(MB)=I+1
          MB=MB+1
        END IF
        IF (I.LT.N-1) THEN
          JH(MB)=I+2
          MB=MB+1
        END IF
        IF (I.LE.K) THEN
          JH(MB)=I+K
          MB=MB+1
        END IF
  270 CONTINUE
      IH(NB+1)=MB
      MB=MB-1
      RETURN
  280 IF (N.LT.6) GO TO 570
      DO 290 I=1,N
        X(I)=1.0D0/DBLE(N)
  290 CONTINUE
      FMIN=-ETA9
      GO TO 260
  300 IF (N.LT.6) GO TO 570
      DO 310 I=1,N
        X(I)=1.0D0
  310 CONTINUE
      FMIN=-ETA9
      GO TO 260
  320 IF (N.LT.5) GO TO 570
      N=N-MOD(N,5)
      DO 330 I=0,N-5,5
        X(I+1)=-1.0D0
        X(I+2)=-1.0D0
        X(I+3)=2.0D0
        X(I+4)=-1.0D0
        X(I+5)=-1.0D0
  330 CONTINUE
      X(1)=-2.0D0
      X(2)=2.0D0
      MB=1
      NB=N/5
      DO 340 I=1,NB
        J=5*(I-1)+1
        IH(I)=MB
        MB=MB+5
        JH(MB-5)=J
        JH(MB-4)=J+1
        JH(MB-3)=J+2
        JH(MB-2)=J+3
        JH(MB-1)=J+4
  340 CONTINUE
      IH(NB+1)=MB
      MB=MB-1
      XMAX=1.0D1
      RETURN
  350 IF (N.LT.2) GO TO 570
      N=N-MOD(N,2)
      DO 360 I=2,N,2
        X(I-1)=0.0D0
        X(I)=-1.0D0
  360 CONTINUE
      XMAX=1.0D1
      GO TO 30
  370 IF (N.LT.2) GO TO 570
      N=N-MOD(N,2)
      DO 380 I=2,N,2
        X(I-1)=-1.0D0
        X(I)=1.0D0
  380 CONTINUE
!      XMAX=1.0D 0
      GO TO 30
  390 IF (N.LT.3) GO TO 570
      P=1.0D0/DBLE(N+1)
      DO 400 I=1,N
        Q=P*DBLE(I)
        X(I)=Q*(Q-1.0D0)
  400 CONTINUE
!      XMAX=1.0D 1
      GO TO 150
  410 IF (N.LT.2) GO TO 570
      P=1.0D0/DBLE(N+1)
      DO 420 I=1,N
        X(I)=DBLE(I)*DBLE(N+1-I)*P**2
  420 CONTINUE
      FMIN=-ETA9
!      XMAX=1.0D 1
      GO TO 30
  430 IF (N.LT.3) GO TO 570
      DO 440 I=1,N
        X(I)=1.0D0
  440 CONTINUE
      FMIN=-ETA9
      GO TO 150
  450 IF (N.LT.3) GO TO 570
      P=1.0D0/DBLE(N+1)
      DO 460 I=1,N
        X(I)=DBLE(I)*DBLE(N+1-I)*P**2
  460 CONTINUE
      FMIN=-ETA9
      GO TO 150
  470 IF (N.LT.3) GO TO 570
      P=1.0D0/DBLE(N+1)
      DO 480 I=1,N
        X(I)=DBLE(I)*DBLE(N+1-I)*P**2
  480 CONTINUE
      FMIN=-ETA9
      GO TO 150
  490 IF (N.LT.3) GO TO 570
      P=EXP(2.0D0)/DBLE(N+1)
      DO 500 I=1,N
        X(I)=(P*DBLE(I)+1.0D0)/3.0D0
  500 CONTINUE
      FMIN=-ETA9
      GO TO 150
  510 IF (N.LT.3) GO TO 570
      P=1.0D0/DBLE(N+1)
      DO 520 I=1,N
        X(I)=P*DBLE(I)
  520 CONTINUE
      FMIN=-ETA9
      GO TO 150
  530 IF (N.LT.3) GO TO 570
      P=1.0D0/DBLE(N+1)
      DO 540 I=1,N
        X(I)=P*DBLE(I)+1.0D0
  540 CONTINUE
      FMIN=-ETA9
      GO TO 150
  550 IF (N.LT.3) GO TO 570
      P=1.0D0/DBLE(N+1)
      DO 560 I=1,N
        X(I)=DBLE(I)*DBLE(N+1-I)*P**2
  560 CONTINUE
      FMIN=-ETA9
      GO TO 150
  570 IERR=1
      RETURN
      END
! SUBROUTINE TIUD14                ALL SYSTEMS                99/12/01
! PORTABILITY : ALL SYSTEMS
! 98/12/01 TU : ORIGINAL VERSION
!
! PURPOSE :
!  INITIAL VALUES OF THE VARIABLES AND STRUCTURE OF THE SPARSE HESSIAN
!  MATRIX FOR UNCONSTRAINED MINIMIZATION.
!
! PARAMETERS :
!  II  N  NUMBER OF VARIABLES.
!  RO  X(N)  VECTOR OF VARIABLES.
!  II  NEXT  NUMBER OF THE TEST PROBLEM.
!  IO  IERR  ERROR INDICATOR.
!
      SUBROUTINE TIUD14 (N, X, FMIN, XMAX, NEXT, IERR)
      INTEGER N,NEXT,IERR
      DOUBLE PRECISION X(*),FMIN,XMAX
      DOUBLE PRECISION P,Q
      INTEGER I
      DOUBLE PRECISION ETA9
      PARAMETER  (ETA9=1.0D60)
      FMIN=0.0D0
      XMAX=1.0D3
      IERR=0
      GO TO (10,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,
     &330,350,370,390,410,430),NEXT
   10 IF (N.LT.2) GO TO 450
      N=N-MOD(N,2)
      DO 20 I=1,N
        IF (MOD(I,2).EQ.1) THEN
          X(I)=-1.2D0
        ELSE
          X(I)=1.0D0
        END IF
   20 CONTINUE
      RETURN
   30 IF (N.LT.4) GO TO 450
      N=N-MOD(N,2)
      DO 40 I=1,N
        IF (MOD(I,2).EQ.1) THEN
          X(I)=-2.0D0
          IF (I.LE.4) X(I)=-3.0D0
        ELSE
          X(I)=0.0D0
          IF (I.LE.4) X(I)=-1.0D0
        END IF
   40 CONTINUE
      RETURN
   50 IF (N.LT.4) GO TO 450
      N=N-MOD(N,2)
      DO 60 I=1,N
        IF (MOD(I,4).EQ.1) THEN
          X(I)=3.0D0
        ELSE IF (MOD(I,4).EQ.2) THEN
          X(I)=-1.0D0
        ELSE IF (MOD(I,4).EQ.3) THEN
          X(I)=0.0D0
        ELSE
          X(I)=1.0D0
        END IF
   60 CONTINUE
      RETURN
   70 IF (N.LT.4) GO TO 450
      N=N-MOD(N,2)
      DO 80 I=1,N
        X(I)=2.0D0
   80 CONTINUE
      X(1)=1.0D0
      RETURN
   90 IF (N.LT.3) GO TO 450
      DO 100 I=1,N
        X(I)=-1.0D0
  100 CONTINUE
      RETURN
  110 IF (N.LT.7) GO TO 450
      DO 120 I=1,N
        X(I)=-1.0D0
  120 CONTINUE
      RETURN
  130 IF (N.LT.4) GO TO 450
      N=N-MOD(N,2)
      DO 140 I=1,N
        X(I)=-1.0D0
  140 CONTINUE
      RETURN
  150 IF (N.LT.6) GO TO 450
      DO 160 I=1,N
        X(I)=1.0D0/DBLE(N)
  160 CONTINUE
      RETURN
  170 IF (N.LT.6) GO TO 450
      DO 180 I=1,N
        X(I)=1.0D0/DBLE(N)
  180 CONTINUE
      FMIN=-ETA9
      RETURN
  190 IF (N.LT.6) GO TO 450
      DO 200 I=1,N
        X(I)=1.0D0
  200 CONTINUE
      FMIN=-ETA9
      RETURN
  210 IF (N.LT.5) GO TO 450
      N=N-MOD(N,5)
      DO 220 I=0,N-5,5
        X(I+1)=-1.0D0
        X(I+2)=-1.0D0
        X(I+3)=2.0D0
        X(I+4)=-1.0D0
        X(I+5)=-1.0D0
  220 CONTINUE
      X(1)=-2.0D0
      X(2)=2.0D0
      XMAX=1.0D0
      RETURN
  230 IF (N.LT.2) GO TO 450
      N=N-MOD(N,2)
      DO 240 I=2,N,2
        X(I-1)=0.0D0
        X(I)=-1.0D0
  240 CONTINUE
      XMAX=1.0D0
      RETURN
  250 IF (N.LT.2) GO TO 450
      N=N-MOD(N,2)
      DO 260 I=2,N,2
        X(I-1)=-1.0D0
        X(I)=1.0D0
  260 CONTINUE
      XMAX=1.0D0
      RETURN
  270 IF (N.LT.3) GO TO 450
      P=1.0D0/DBLE(N+1)
      DO 280 I=1,N
        Q=P*DBLE(I)
        X(I)=Q*(Q-1.0D0)
  280 CONTINUE
      RETURN
  290 IF (N.LT.3) GO TO 450
      P=1.0D0/DBLE(N+1)
      DO 300 I=1,N
        X(I)=DBLE(I)*DBLE(N+1-I)*P**2
  300 CONTINUE
      FMIN=-ETA9
      XMAX=1.0D1
      RETURN
  310 IF (N.LT.3) GO TO 450
      DO 320 I=1,N
        X(I)=1.0D0
  320 CONTINUE
      FMIN=-ETA9
      RETURN
  330 IF (N.LT.3) GO TO 450
      P=1.0D0/DBLE(N+1)
      DO 340 I=1,N
        X(I)=DBLE(I)*DBLE(N+1-I)*P**2
  340 CONTINUE
      FMIN=-ETA9
      RETURN
  350 IF (N.LT.3) GO TO 450
      P=1.0D0/DBLE(N+1)
      DO 360 I=1,N
        X(I)=DBLE(I)*DBLE(N+1-I)*P**2
  360 CONTINUE
      FMIN=-ETA9
      RETURN
  370 IF (N.LT.3) GO TO 450
      P=EXP(2.0D0)/DBLE(N+1)
      DO 380 I=1,N
        X(I)=(P*DBLE(I)+1.0D0)/3.0D0
  380 CONTINUE
      FMIN=-ETA9
      RETURN
  390 IF (N.LT.3) GO TO 450
      P=1.0D0/DBLE(N+1)
      DO 400 I=1,N
        X(I)=P*DBLE(I)
  400 CONTINUE
      FMIN=-ETA9
      RETURN
  410 IF (N.LT.3) GO TO 450
      P=1.0D0/DBLE(N+1)
      DO 420 I=1,N
        X(I)=P*DBLE(I)+1.0D0
  420 CONTINUE
      FMIN=-ETA9
      RETURN
  430 IF (N.LT.3) GO TO 450
      P=1.0D0/DBLE(N+1)
      DO 440 I=1,N
        X(I)=DBLE(I)*DBLE(N+1-I)*P**2
  440 CONTINUE
      FMIN=-ETA9
      RETURN
  450 IERR=1
      RETURN
      END
! SUBROUTINE TIUS14                ALL SYSTEMS                99/12/01
! PORTABILITY : ALL SYSTEMS
! 98/12/01 TU : ORIGINAL VERSION
!
! PURPOSE :
!  INITIAL VALUES OF THE VARIABLES AND STRUCTURE OF THE SPARSE HESSIAN
!  MATRIX FOR UNCONSTRAINED MINIMIZATION.
!  SPARSE VERSION WITH CHANGED TESTS 7-10, 12.
!  CHANGED FOR THE TESTS.
!
! PARAMETERS :
!  II  N  NUMBER OF VARIABLES.
!  II  M  NUMBER OF ELEMENTS OF THE SPARSE MATRIX.
!  RO  X(N)  VECTOR OF VARIABLES.
!  IO  IH(N+1)  POINTERS OF THE DIAGONAL ELEMENTS OF THE HESSIAN MATRIX.
!  IO  JH(M)  INDICES OF THE NONZERO ELEMENTS OF THE HESSIAN MATRIX IN
!             THE PACKED ROW.
!  II  NEXT  NUMBER OF THE TEST PROBLEM.
!  IO  IERR  ERROR INDICATOR.
!
      SUBROUTINE TIUS14 (N, M, X, IH, JH, FMIN, XMAX, NEXT, IERR)
      INTEGER N,M,NEXT,IERR
      INTEGER IH(*),JH(*)
      DOUBLE PRECISION X(N),FMIN,XMAX
      DOUBLE PRECISION P,Q
      INTEGER I,J,K,L,K1,K2,M1
      DOUBLE PRECISION ETA9
      PARAMETER  (ETA9=1.0D60)
      FMIN=0.0D0
      XMAX=1.0D3
      IERR=0
      GO TO (10,50,80,120,150,190,220,250,280,310,340,390,410,430,450,
     &470,490,510,530,550,570,590),NEXT
   10 IF (N.LT.2) GO TO 610
      N=N-MOD(N,2)
      DO 20 I=1,N
        IF (MOD(I,2).EQ.1) THEN
          X(I)=-1.2D0
        ELSE
          X(I)=1.0D0
        END IF
   20 CONTINUE
   30 DO 40 I=1,N-1
        J=2*(I-1)+1
        IH(I)=J
        JH(J)=I
        JH(J+1)=I+1
   40 CONTINUE
      J=2*(N-1)+1
      IH(N)=J
      JH(J)=N
      IH(N+1)=2*N
      M=2*N-1
      RETURN
   50 IF (N.LT.4) GO TO 610
      N=N-MOD(N,2)
      DO 60 I=1,N
        IF (MOD(I,2).EQ.1) THEN
          X(I)=-2.0D0
          IF (I.LE.4) X(I)=-3.0D0
        ELSE
          X(I)=0.0D0
          IF (I.LE.4) X(I)=-1.0D0
        END IF
   60 CONTINUE
      DO 70 I=1,N-1
        J=2*I-1
        IH(I)=J
        JH(J)=I
        JH(J+1)=I+1
        IF (MOD(I,2).EQ.0) JH(J+1)=JH(J+1)+1
   70 CONTINUE
      J=2*N-1
      IH(N)=J
      JH(J)=N
      M=J
      IH(N+1)=J+1
      RETURN
   80 IF (N.LT.4) GO TO 610
      N=N-MOD(N,2)
      DO 90 I=1,N
        IF (MOD(I,4).EQ.1) THEN
          X(I)=3.0D0
        ELSE IF (MOD(I,4).EQ.2) THEN
          X(I)=-1.0D0
        ELSE IF (MOD(I,4).EQ.3) THEN
          X(I)=0.0D0
        ELSE
          X(I)=1.0D0
        END IF
   90 CONTINUE
      IH(1)=1
      IH(2)=4
      JH(1)=1
      JH(2)=2
      JH(3)=4
      JH(4)=2
      JH(5)=3
      K=5
      DO 110 I=3,N-3,2
        IH(I)=IH(I-2)+5
        J=IH(I)
        DO 100 L=J,J+4
          JH(L)=JH(L-5)+2
          K=K+1
  100   CONTINUE
        IH(I+1)=IH(I-1)+5
  110 CONTINUE
      IH(N-1)=IH(N-3)+5
      IH(N)=IH(N-2)+4
      IH(N+1)=IH(N)+1
      JH(K+1)=JH(K)
      JH(K+2)=JH(K)+1
      JH(K+3)=JH(K+2)
      M=IH(N)
      RETURN
  120 IF (N.LT.4) GO TO 610
      N=N-MOD(N,2)
      DO 130 I=1,N
        X(I)=2.0D0
  130 CONTINUE
      X(1)=1.0D0
      DO 140 I=1,N-1
        IH(I)=2*I-1
        J=IH(I)
        JH(J)=I
        JH(J+1)=I+1
  140 CONTINUE
      IH(N)=2*N-1
      JH(J+2)=JH(J+1)
      IH(N+1)=2*N
      M=IH(N)
      RETURN
  150 IF (N.LT.3) GO TO 610
      DO 160 I=1,N
        X(I)=-1.0D0
  160 CONTINUE
  170 DO 180 I=1,N-2
        J=1+3*(I-1)
        IH(I)=J
        JH(J)=I
        JH(J+1)=I+1
        JH(J+2)=I+2
  180 CONTINUE
      J=3*(N-2)+1
      IH(N-1)=J
      JH(J)=N-1
      JH(J+1)=N
      J=IH(N-1)+2
      IH(N)=J
      JH(J)=N
      IH(N+1)=IH(N)+1
      M=3*N-3
      RETURN
  190 IF (N.LT.7) GO TO 610
      DO 200 I=1,N
        X(I)=-1.0D0
  200 CONTINUE
      DO 210 I=1,N-6
        J=7*(I-1)+1
        IH(I)=J
        JH(J)=I
        JH(J+1)=I+1
        JH(J+2)=I+2
        JH(J+3)=I+3
        JH(J+4)=I+4
        JH(J+5)=I+5
        JH(J+6)=I+6
  210 CONTINUE
      J=7*(N-6)+1
      IH(N-5)=J
      JH(J)=N-5
      JH(J+1)=N-4
      JH(J+2)=N-3
      JH(J+3)=N-2
      JH(J+4)=N-1
      JH(J+5)=N
      IH(N-4)=J+6
      JH(J+6)=N-4
      JH(J+7)=N-3
      JH(J+8)=N-2
      JH(J+9)=N-1
      JH(J+10)=N
      IH(N-3)=J+11
      JH(J+11)=N-3
      JH(J+12)=N-2
      JH(J+13)=N-1
      JH(J+14)=N
      IH(N-2)=J+15
      JH(J+15)=N-2
      JH(J+16)=N-1
      JH(J+17)=N
      IH(N-1)=J+18
      JH(J+18)=N-1
      JH(J+19)=N
      IH(N)=J+20
      JH(J+20)=N
      IH(N+1)=J+21
      M=J+20
      RETURN
  220 IF (N.LT.4) GO TO 610
      N=N-MOD(N,2)
      DO 230 I=1,N
        X(I)=-1.0D0
  230 CONTINUE
      M=0
      K=N/2
      IH(1)=1
      DO 240 I=1,N
        M=M+1
        JH(M)=I
        IF (I.LT.N) THEN
          M=M+1
          JH(M)=I+1
        END IF
        IF (I.LE.K) THEN
          M=M+1
          JH(M)=I+K
        END IF
        IH(I+1)=M+1
  240 CONTINUE
      RETURN
  250 IF (N.LT.10) GO TO 610
      DO 260 I=1,N
        X(I)=1.0D0/DBLE(N)
  260 CONTINUE
      M=0
      K=N/2
      IH(1)=1
      DO 270 I=1,N
        M=M+1
        JH(M)=I
        IF (I+1.LE.N) THEN
          M=M+1
          JH(M)=I+1
        END IF
        IF (I+2.LE.N) THEN
          M=M+1
          JH(M)=I+2
        END IF
        IF (I+3.LE.N) THEN
          M=M+1
          JH(M)=I+3
        END IF
        IF (I+4.LE.N) THEN
          M=M+1
          JH(M)=I+4
        END IF
        IF (I+K-2.LE.N) THEN
          M=M+1
          JH(M)=I+K-2
        END IF
        IF (I+K-1.LE.N) THEN
          M=M+1
          JH(M)=I+K-1
        END IF
        IF (I+K.LE.N) THEN
          M=M+1
          JH(M)=I+K
        END IF
        IF (I+K+1.LE.N) THEN
          M=M+1
          JH(M)=I+K+1
        END IF
        IF (I+K+2.LE.N) THEN
          M=M+1
          JH(M)=I+K+2
        END IF
        IH(I+1)=M+1
  270 CONTINUE
      RETURN
  280 IF (N.LT.4) GO TO 610
      DO 290 I=1,N
        X(I)=1.0D0/DBLE(N)
  290 CONTINUE
      FMIN=-ETA9
      M=0
      K=N/2
      IH(1)=1
      DO 300 I=1,N
        M=M+1
        JH(M)=I
        IH(I+1)=M+1
  300 CONTINUE
      RETURN
  310 IF (N.LT.10) GO TO 610
      DO 320 I=1,N
        X(I)=1.0D0
  320 CONTINUE
      FMIN=-ETA9
      M=0
      K=N/2
      IH(1)=1
      DO 330 I=1,N
        M=M+1
        JH(M)=I
        IF (I.LT.N) THEN
          M=M+1
          JH(M)=I+1
        END IF
        IF (I.LT.N-1) THEN
          M=M+1
          JH(M)=I+2
        END IF
        IF (I.LE.K) THEN
          M=M+1
          JH(M)=I+K
        END IF
        IH(I+1)=M+1
  330 CONTINUE
      RETURN
  340 IF (N.LT.5) GO TO 610
      N=N-MOD(N,5)
      DO 350 I=0,N-5,5
        X(I+1)=-1.0D0
        X(I+2)=-1.0D0
        X(I+3)=2.0D0
        X(I+4)=-1.0D0
        X(I+5)=-1.0D0
  350 CONTINUE
      X(1)=-2.0D0
      X(2)=2.0D0
      IH(1)=1
      K=1
      DO 380 I=1,N,5
        M1=I
        K1=I
        K2=I+4
        DO 370 M=1,5
          DO 360 J=K1,K2
            JH(K)=J
            K=K+1
  360     CONTINUE
          K1=K1+1
          M1=M1+1
          IH(M1)=6-M+IH(M1-1)
  370   CONTINUE
  380 CONTINUE
      M=IH(N)
      XMAX=1.0D0
      RETURN
  390 IF (N.LT.2) GO TO 610
      N=N-MOD(N,2)
      DO 400 I=2,N,2
        X(I-1)=0.0D0
        X(I)=-1.0D0
  400 CONTINUE
      XMAX=1.0D0
      GO TO 30
  410 IF (N.LT.2) GO TO 610
      N=N-MOD(N,2)
      DO 420 I=2,N,2
        X(I-1)=-1.0D0
        X(I)=1.0D0
  420 CONTINUE
      GO TO 30
  430 IF (N.LT.3) GO TO 610
      P=1.0D0/DBLE(N+1)
      DO 440 I=1,N
        Q=P*DBLE(I)
        X(I)=Q*(Q-1.0D0)
  440 CONTINUE
      GO TO 170
  450 IF (N.LT.3) GO TO 610
      P=1.0D0/DBLE(N+1)
      DO 460 I=1,N
        X(I)=DBLE(I)*DBLE(N+1-I)*P**2
  460 CONTINUE
      FMIN=-ETA9
      XMAX=1.0D1
      GO TO 30
  470 IF (N.LT.3) GO TO 610
      DO 480 I=1,N
        X(I)=1.0D0
  480 CONTINUE
      FMIN=-ETA9
      GO TO 170
  490 IF (N.LT.3) GO TO 610
      P=1.0D0/DBLE(N+1)
      DO 500 I=1,N
        X(I)=DBLE(I)*DBLE(N+1-I)*P**2
  500 CONTINUE
      FMIN=-ETA9
      GO TO 170
  510 IF (N.LT.3) GO TO 610
      P=1.0D0/DBLE(N+1)
      DO 520 I=1,N
        X(I)=DBLE(I)*DBLE(N+1-I)*P**2
  520 CONTINUE
      FMIN=-ETA9
      GO TO 170
  530 IF (N.LT.3) GO TO 610
      P=EXP(2.0D0)/DBLE(N+1)
      DO 540 I=1,N
        X(I)=(P*DBLE(I)+1.0D0)/3.0D0
  540 CONTINUE
      FMIN=-ETA9
      GO TO 170
  550 IF (N.LT.3) GO TO 610
      P=1.0D0/DBLE(N+1)
      DO 560 I=1,N
        X(I)=P*DBLE(I)
  560 CONTINUE
      FMIN=-ETA9
      GO TO 170
  570 IF (N.LT.3) GO TO 610
      P=1.0D0/DBLE(N+1)
      DO 580 I=1,N
        X(I)=P*DBLE(I)+1.0D0
  580 CONTINUE
      FMIN=-ETA9
      GO TO 170
  590 IF (N.LT.3) GO TO 610
      P=1.0D0/DBLE(N+1)
      DO 600 I=1,N
        X(I)=DBLE(I)*DBLE(N+1-I)*P**2
  600 CONTINUE
      FMIN=-ETA9
      GO TO 170
  610 IERR=1
      RETURN
      END
! SUBROUTINE TAFU14             ALL SYSTEMS                92/12/01
! PORTABILITY : ALL SYSTEMS
! 92/12/01 RA : ORIGINAL VERSION
!
! PURPOSE :
!  VALUES OF MODEL FUNCTIONS FOR UNCONSTRAINED MINIMIZATION.
!  SPARSE VERSION - LU VERSION WITH MODIFIED TESTS NO 7-10,12.
!
! PARAMETERS :
!  II  N  NUMBER OF VARIABLES.
!  RI  X(N)  VECTOR OF VARIABLES.
!  RO  F  VALUE OF THE MODEL FUNCTION.
!  II  NEXT  NUMBER OF THE TEST PROBLEM.
!
      SUBROUTINE TAFU14 (N, KA, X, FA, NEXT)
      INTEGER N,NEXT
      DOUBLE PRECISION X(*),FA
      DOUBLE PRECISION A,B,C,D,P,Q,R,U,V
      INTEGER I,J,K,KA
      GO TO (10,20,30,40,50,60,80,90,110,130,150,170,180,190,200,210,
     &220,230,240,250,260,270),NEXT
   10 A=X(KA)**2-X(KA+1)
      B=X(KA)-1.0D0
      FA=1.0D2*A**2+B**2
      RETURN
   20 I=2*KA-1
      A=X(I)**2-X(I+1)
      B=X(I)-1.0D0
      C=X(I+2)**2-X(I+3)
      D=X(I+2)-1.0D0
      U=X(I+1)+X(I+3)-2.0D0
      V=X(I+1)-X(I+3)
      FA=1.0D2*A**2+B**2+9.0D1*C**2+D**2+1.0D1*U**2+0.1D0*V**2
      RETURN
   30 I=2*KA-1
      A=X(I)+1.0D1*X(I+1)
      B=X(I+2)-X(I+3)
      C=X(I+1)-2.0D0*X(I+2)
      D=X(I)-X(I+3)
      FA=A**2+5.0D0*B**2+C**4+1.0D1*D**4
      RETURN
   40 I=2*KA-1
      A=EXP(X(I))
      B=A-X(I+1)
      D=X(I+1)-X(I+2)
      P=X(I+2)-X(I+3)
      Q=SIN(P)/COS(P)
      U=X(I)
      V=X(I+3)-1.0D0
      FA=B**4+1.0D2*D**6+Q**4+U**8+V**2
      RETURN
   50 P=7.0D0/3.0D0
      A=(3.0D0-2.0D0*X(KA))*X(KA)+1.0D0
      IF (KA.GT.1) A=A-X(KA-1)
      IF (KA.LT.N) A=A-X(KA+1)
      FA=ABS(A)**P
      RETURN
   60 P=7.0D0/3.0D0
      A=(2.0D0+5.0D0*X(KA)**2)*X(KA)+1.0D0
      DO 70 I=MAX(1,KA-5),MIN(N,KA+1)
        A=A+X(I)*(1.0D0+X(I))
   70 CONTINUE
      FA=ABS(A)**P
      RETURN
   80 P=7.0D0/3.0D0
      IF (KA.LE.N) THEN
        A=(3.0D0-2.0D0*X(KA))*X(KA)+1.0D0
        IF (KA.GT.1) A=A-X(KA-1)
        IF (KA.LT.N) A=A-X(KA+1)
      ELSE
        I=KA-N
        A=X(I)+X(I+N/2)
      END IF
      FA=ABS(A)**P
      RETURN
   90 K=N/2
      P=0.0D0
      DO 100 I=KA-2,KA+2
        IF (I.LT.1.OR.I.GT.N) GO TO 100
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=DBLE(I+KA)/1.0D1
        P=P+A*SIN(X(I))+B*COS(X(I))
  100 CONTINUE
      IF (KA.GT.K) THEN
        I=KA-K
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=DBLE(I+KA)/1.0D1
        P=P+A*SIN(X(I))+B*COS(X(I))
      ELSE
        I=KA+K
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=DBLE(I+KA)/1.0D1
        P=P+A*SIN(X(I))+B*COS(X(I))
      END IF
      FA=(DBLE(N+KA)-P)**2/DBLE(N)
      RETURN
  110 K=N/2
      FA=DBLE(KA)*(1.0D0-COS(X(KA)))
      DO 120 I=KA-2,KA+2
        IF (I.LT.1.OR.I.GT.N) GO TO 120
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=DBLE(I+KA)/1.0D1
        FA=FA+A*SIN(X(I))+B*COS(X(I))
  120 CONTINUE
      IF (KA.GT.K) THEN
        I=KA-K
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=DBLE(I+KA)/1.0D1
        FA=FA+A*SIN(X(I))+B*COS(X(I))
      ELSE
        I=KA+K
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=DBLE(I+KA)/1.0D1
        FA=FA+A*SIN(X(I))+B*COS(X(I))
      END IF
      FA=FA/DBLE(N)
      RETURN
  130 K=N/2
      FA=0.0D0
      Q=1.0D0+DBLE(KA)/1.0D1
      DO 140 I=KA-2,KA+2
        IF (I.LT.1.OR.I.GT.N) GO TO 140
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=1.0D0+DBLE(I)/1.0D1
        C=DBLE(I+KA)/1.0D1
        FA=FA+A*SIN(Q*X(KA)+B*X(I)+C)
  140 CONTINUE
      IF (KA.GT.K) THEN
        I=KA-K
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=1.0D0+DBLE(I)/1.0D1
        C=DBLE(I+KA)/1.0D1
        FA=FA+A*SIN(Q*X(KA)+B*X(I)+C)
      ELSE
        I=KA+K
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=1.0D0+DBLE(I)/1.0D1
        C=DBLE(I+KA)/1.0D1
        FA=FA+A*SIN(Q*X(KA)+B*X(I)+C)
      END IF
      FA=FA/DBLE(N)
      RETURN
  150 P=-0.2008D-2
      Q=-0.1900D-2
      R=-0.0261D-2
      I=5*(KA-1)
      A=1.0D0
      B=0.0D0
      DO 160 J=1,5
        A=A*X(I+J)
        B=B+X(I+J)**2
  160 CONTINUE
      A=EXP(A)
      B=B-1.0D1-P
      C=X(I+2)*X(I+3)-5.0D0*X(I+4)*X(I+5)-Q
      D=X(I+1)**3+X(I+2)**3+1.0D0-R
      FA=A+1.0D1*(B**2+C**2+D**2)
      RETURN
  170 A=X(KA)-3.0D0
      B=X(KA)-X(KA+1)
      FA=A**2+B**2+EXP(2.0D1*B)
      RETURN
  180 A=X(KA+1)**2
      B=X(KA)**2
      C=A+1.0D0
      D=B+1.0D0
      FA=B**C+A**D
      RETURN
  190 P=1.0D0/DBLE(N+1)
      Q=0.5D0*P**2
      A=2.0D0*X(KA)+Q*(X(KA)+DBLE(KA)*P+1.0D0)**3
      IF (KA.GT.1) A=A-X(KA-1)
      IF (KA.LT.N) A=A-X(KA+1)
      FA=A**2
      RETURN
  200 P=1.0D0/DBLE(N+1)
      Q=2.0D0/P
      R=2.0D0*P
      A=X(KA)-X(KA+1)
      FA=Q*X(KA)*A
      IF (ABS(A).LE.1.0D-6) THEN
        FA=FA+R*EXP(X(KA+1))*(1.0D0+A/2.0D0*(1.0D0+A/3.0D0*(1.0D0+A/
     &   4.0D0)))
      ELSE
        B=EXP(X(KA))-EXP(X(KA+1))
        FA=FA+R*B/A
      END IF
      IF (KA.EQ.1) THEN
        FA=FA+R*(EXP(X(1))-1.0D0)/X(1)
      ELSE IF (KA.EQ.N-1) THEN
        FA=FA+Q*X(N)**2+R*(EXP(X(N))-1.0D0)/X(N)
      END IF
      RETURN
  210 A=DBLE(KA)*(1.0D0-COS(X(KA)))
      IF (KA.GT.1) A=A+DBLE(KA)*SIN(X(KA-1))
      IF (KA.LT.N) A=A-DBLE(KA)*SIN(X(KA+1))
      FA=A
      RETURN
  220 P=1.0D0/DBLE(N+1)
      IF (KA.EQ.1) THEN
        FA=0.25D0*X(KA)**2/P+1.25D-1*X(KA+1)**2/P+P*(EXP(X(KA))-1.0D0)
      ELSE IF (KA.EQ.N) THEN
        FA=0.25D0*X(KA)**2/P+1.25D-1*X(KA-1)**2/P+P*(EXP(X(KA))-1.0D0)
      ELSE
        FA=1.25D-1*(X(KA+1)-X(KA-1))**2/P+P*(EXP(X(KA))-1.0D0)
      END IF
      RETURN
  230 P=1.0D0/DBLE(N+1)
      Q=DBLE(KA)*P
      IF (KA.EQ.1) THEN
        FA=0.5D0*X(KA)**2/P+0.25D0*X(KA+1)**2/P-P*(X(KA)**2+2.0D0*X(KA)*
     &   Q)
      ELSE IF (KA.EQ.N) THEN
        FA=0.5D0*X(KA)**2/P+0.25D0*X(KA-1)**2/P-P*(X(KA)**2+2.0D0*X(KA)*
     &   Q)
      ELSE
        FA=2.5D-1*(X(KA+1)-X(KA-1))**2/P-P*(X(KA)**2+2.0D0*X(KA)*Q)
      END IF
      RETURN
  240 P=1.0D0/DBLE(N+1)
      Q=EXP(2.0D0*DBLE(KA)*P)
      IF (KA.EQ.1) THEN
        R=1.0D0/3.0D0
        FA=0.5D0*(X(KA)-R)**2/P+7.0D0*R**2+2.5D-1*(X(KA+1)-R)**2/P+P*
     &   (X(KA)**2+2.0D0*X(KA)*Q)
      ELSE IF (KA.EQ.N) THEN
        R=EXP(2.0D0)/3.0D0
        FA=0.5D0*(X(KA)-R)**2/P+7.0D0*R**2+2.5D-1*(X(KA-1)-R)**2/P+P*
     &   (X(KA)**2+2.0D0*X(KA)*Q)
      ELSE
        FA=2.5D-1*(X(KA+1)-X(KA-1))**2/P+P*(X(KA)**2+2.0D0*X(KA)*Q)
      END IF
      RETURN
  250 P=1.0D0/DBLE(N+1)
      IF (KA.EQ.1) THEN
        FA=(0.5D0*X(KA)**2/P-P)+(2.5D-1*X(KA+1)**2/P-P)*EXP(-2.0D0*X(KA)
     &   **2)
      ELSE IF (KA.EQ.N) THEN
        FA=(0.5D0*X(KA)**2/P-P)*EXP(-2.0D0)+(2.5D-1*X(KA-1)**2/P-P)*
     &   EXP(-2.0D0*X(KA)**2)
      ELSE
        FA=(2.5D-1*(X(KA+1)-X(KA-1))**2/P-P)*EXP(-2.0D0*X(KA)**2)
      END IF
      RETURN
  260 P=1.0D0/DBLE(N+1)
      IF (KA.EQ.1) THEN
        A=0.5D0*(X(KA+1)-1.0D0)/P
        B=(X(KA)-1.0D0)/P
        FA=P*(X(KA)**2+A*ATAN(A)-LOG(SQRT(1.0D0+A**2)))+0.5D0*P*(1.0D0+
     &   B*ATAN(B)-LOG(SQRT(1.0D0+B**2)))
      ELSE IF (KA.EQ.N) THEN
        A=0.5D0*(2.0D0-X(KA-1))/P
        B=(2.0D0-X(KA))/P
        FA=P*(X(KA)**2+A*ATAN(A)-LOG(SQRT(1.0D0+A**2)))+0.5D0*P*(4.0D0+
     &   B*ATAN(B)-LOG(SQRT(1.0D0+B**2)))
      ELSE
        A=0.5D0*(X(KA+1)-X(KA-1))/P
        FA=P*(X(KA)**2+A*ATAN(A)-LOG(SQRT(1.0D0+A**2)))
      END IF
      RETURN
  270 P=1.0D0/DBLE(N+1)
      IF (KA.EQ.1) THEN
        A=0.5D0*X(KA+1)/P
        B=X(KA)/P
        FA=P*(1.0D2*(X(KA)-A**2)**2+(1.0D0-A)**2)+0.5D0*P*(1.0D2*B**4+
     &   (1.0D0-B)**2)
      ELSE IF (KA.EQ.N) THEN
        A=-0.5D0*X(KA-1)/P
        B=-X(KA)/P
        FA=P*(1.0D2*(X(KA)-A**2)**2+(1.0D0-A)**2)+0.5D0*P*(1.0D2*B**4+
     &   (1.0D0-B)**2)
      ELSE
        A=0.5D0*(X(KA+1)-X(KA-1))/P
        FA=P*(1.0D2*(X(KA)-A**2)**2+(1.0D0-A)**2)
      END IF
      RETURN
      END
! SUBROUTINE TAGU14                ALL SYSTEMS                92/12/01
! PORTABILITY : ALL SYSTEMS
! 92/12/01 RA : ORIGINAL VERSION
!
! PURPOSE :
!  GRADIENTS OF MODEL FUNCTIONS FOR UNCONSTRAINED MINIMIZATION.
!  SPARSE VERSION - LU VERSION WITH MODIFIED TESTS NO 7-10,12.
!
! PARAMETERS :
!  II  N  NUMBER OF VARIABLES.
!  RI  X(N)  VECTOR OF VARIABLES.
!  RO  G(N)  GRADIENT OF THE MODEL FUNCTION.
!  II  NEXT  NUMBER OF THE TEST PROBLEM.
!
      SUBROUTINE TAGU14 (N, KA, X, GA, NEXT)
      INTEGER N,NEXT
      DOUBLE PRECISION X(*),GA(*)
      DOUBLE PRECISION A,B,C,D,P,Q,R,U,V
      INTEGER I,J,K,KA
      GO TO (10,20,30,40,50,60,90,100,130,150,170,190,200,210,220,230,
     &240,250,260,270,280,290),NEXT
   10 A=X(KA)**2-X(KA+1)
      B=X(KA)-1.0D0
      GA(KA)=4.0D2*X(KA)*A+2.0D0*B
      GA(KA+1)=-2.0D2*A
      RETURN
   20 I=2*KA-1
      A=X(I)**2-X(I+1)
      B=X(I)-1.0D0
      C=X(I+2)**2-X(I+3)
      D=X(I+2)-1.0D0
      U=X(I+1)+X(I+3)-2.0D0
      V=X(I+1)-X(I+3)
      GA(I)=4.0D2*X(I)*A+2.0D0*B
      GA(I+1)=-2.0D2*A+2.0D1*U+0.2D0*V
      GA(I+2)=3.6D2*X(I+2)*C+2.0D0*D
      GA(I+3)=-1.8D2*C+2.0D1*U-0.2D0*V
      RETURN
   30 I=2*KA-1
      A=X(I)+1.0D1*X(I+1)
      B=X(I+2)-X(I+3)
      C=X(I+1)-2.0D0*X(I+2)
      D=X(I)-X(I+3)
      GA(I)=2.0D0*A+4.0D1*D**3
      GA(I+1)=2.0D1*A+4.0D0*C**3
      GA(I+2)=-8.0D0*C**3+1.0D1*B
      GA(I+3)=-4.0D1*D**3-1.0D1*B
      RETURN
   40 I=2*KA-1
      A=EXP(X(I))
      B=A-X(I+1)
      B=4.0D0*B**3
      D=X(I+1)-X(I+2)
      D=6.0D2*D**5
      P=X(I+2)-X(I+3)
      C=COS(P)
      Q=SIN(P)/COS(P)
      Q=4.0D0*Q**3/C**2
      U=X(I)
      V=X(I+3)-1.0D0
      GA(I)=A*B+8.0D0*U**7
      GA(I+1)=D-B
      GA(I+2)=Q-D
      GA(I+3)=2.0D0*V-Q
      RETURN
   50 P=7.0D0/3.0D0
      A=(3.0D0-2.0D0*X(KA))*X(KA)+1.0D0
      IF (KA.GT.1) A=A-X(KA-1)
      IF (KA.LT.N) A=A-X(KA+1)
      B=P*ABS(A)**(P-1.0D0)*SIGN(1.0D0,A)
      GA(KA)=B*(3.0D0-4.0D0*X(KA))
      IF (KA.GT.1) GA(KA-1)=-B
      IF (KA.LT.N) GA(KA+1)=-B
      RETURN
   60 P=7.0D0/3.0D0
      A=(2.0D0+5.0D0*X(KA)**2)*X(KA)+1.0D0
      DO 70 I=MAX(1,KA-5),MIN(N,KA+1)
        A=A+X(I)*(1.0D0+X(I))
   70 CONTINUE
      B=P*ABS(A)**(P-1.0D0)*SIGN(1.0D0,A)
      DO 80 I=MAX(1,KA-5),MIN(N,KA+1)
        GA(I)=B*(1.0D0+2.0D0*X(I))
   80 CONTINUE
      GA(KA)=GA(KA)+B*(2.0D0+1.5D1*X(KA)**2)
      RETURN
   90 P=7.0D0/3.0D0
      IF (KA.LE.N) THEN
        A=(3.0D0-2.0D0*X(KA))*X(KA)+1.0D0
        IF (KA.GT.1) A=A-X(KA-1)
        IF (KA.LT.N) A=A-X(KA+1)
        B=P*ABS(A)**(P-1.0D0)*SIGN(1.0D0,A)
        GA(KA)=B*(3.0D0-4.0D0*X(KA))
        IF (KA.GT.1) GA(KA-1)=-B
        IF (KA.LT.N) GA(KA+1)=-B
      ELSE
        I=KA-N
        A=X(I)+X(I+N/2)
        B=P*ABS(A)**(P-1.0D0)*SIGN(1.0D0,A)
        GA(I)=B
        GA(I+N/2)=B
      END IF
      RETURN
  100 K=N/2
      P=0.0D0
      DO 110 I=KA-2,KA+2
        IF (I.LT.1.OR.I.GT.N) GO TO 110
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=DBLE(I+KA)/1.0D1
        P=P+A*SIN(X(I))+B*COS(X(I))
  110 CONTINUE
      IF (KA.GT.K) THEN
        I=KA-K
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=DBLE(I+KA)/1.0D1
        P=P+A*SIN(X(I))+B*COS(X(I))
      ELSE
        I=KA+K
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=DBLE(I+KA)/1.0D1
        P=P+A*SIN(X(I))+B*COS(X(I))
      END IF
      P=2.0D0*(DBLE(N+KA)-P)/DBLE(N)
      DO 120 I=KA-2,KA+2
        IF (I.LT.1.OR.I.GT.N) GO TO 120
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=DBLE(I+KA)/1.0D1
        GA(I)=-P*(A*COS(X(I))-B*SIN(X(I)))
  120 CONTINUE
      IF (KA.GT.K) THEN
        I=KA-K
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=DBLE(I+KA)/1.0D1
        GA(I)=-P*(A*COS(X(I))-B*SIN(X(I)))
      ELSE
        I=KA+K
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=DBLE(I+KA)/1.0D1
        GA(I)=-P*(A*COS(X(I))-B*SIN(X(I)))
      END IF
      RETURN
  130 K=N/2
      P=1.0D0/DBLE(N)
      DO 140 I=KA-2,KA+2
        IF (I.LT.1.OR.I.GT.N) GO TO 140
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=DBLE(I+KA)/1.0D1
        GA(I)=P*(A*COS(X(I))-B*SIN(X(I)))
  140 CONTINUE
      IF (KA.GT.K) THEN
        I=KA-K
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=DBLE(I+KA)/1.0D1
        GA(I)=P*(A*COS(X(I))-B*SIN(X(I)))
      ELSE
        I=KA+K
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=DBLE(I+KA)/1.0D1
        GA(I)=P*(A*COS(X(I))-B*SIN(X(I)))
      END IF
      GA(KA)=GA(KA)+P*DBLE(KA)*SIN(X(KA))
      RETURN
  150 K=N/2
      GA(KA)=0.0D0
      Q=1.0D0+DBLE(KA)/1.0D1
      DO 160 I=KA-2,KA+2
        IF (I.LT.1.OR.I.GT.N) GO TO 160
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=1.0D0+DBLE(I)/1.0D1
        C=DBLE(I+KA)/1.0D1
        P=A*COS(Q*X(KA)+B*X(I)+C)/DBLE(N)
        GA(KA)=GA(KA)+P*Q
        IF (I.EQ.KA) THEN
          GA(I)=GA(I)+P*B
        ELSE
          GA(I)=P*B
        END IF
  160 CONTINUE
      IF (KA.GT.K) THEN
        I=KA-K
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=1.0D0+DBLE(I)/1.0D1
        C=DBLE(I+KA)/1.0D1
        P=A*COS(Q*X(KA)+B*X(I)+C)/DBLE(N)
        GA(KA)=GA(KA)+P*Q
        GA(I)=P*B
      ELSE
        I=KA+K
        A=5.0D0*(1.0D0+MOD(I,5)+MOD(KA,5))
        B=1.0D0+DBLE(I)/1.0D1
        C=DBLE(I+KA)/1.0D1
        P=A*COS(Q*X(KA)+B*X(I)+C)/DBLE(N)
        GA(KA)=GA(KA)+P*Q
        GA(I)=P*B
      END IF
      RETURN
  170 P=-0.2008D-2
      Q=-0.1900D-2
      R=-0.0261D-2
      I=5*(KA-1)
      A=1.0D0
      B=0.0D0
      DO 180 J=1,5
        A=A*X(I+J)
        B=B+X(I+J)**2
  180 CONTINUE
      A=A*EXP(A)
      B=B-1.0D1-P
      C=X(I+2)*X(I+3)-5.0D0*X(I+4)*X(I+5)-Q
      D=X(I+1)**3+X(I+2)**3+1.0D0-R
      GA(I+1)=A/X(I+1)+2.0D1*(2.0D0*B*X(I+1)+3.0D0*D*X(I+1)**2)
      GA(I+2)=A/X(I+2)+2.0D1*(2.0D0*B*X(I+2)+C*X(I+3)+3.0D0*D*X(I+2)**2)
      GA(I+3)=A/X(I+3)+2.0D1*(2.0D0*B*X(I+3)+C*X(I+2))
      GA(I+4)=A/X(I+4)+2.0D1*(2.0D0*B*X(I+4)-5.0D0*C*X(I+5))
      GA(I+5)=A/X(I+5)+2.0D1*(2.0D0*B*X(I+5)-5.0D0*C*X(I+4))
      RETURN
  190 A=X(KA)-3.0D0
      B=X(KA)-X(KA+1)
      GA(KA)=2.0D0*A+2.0D0*B+2.0D1*EXP(2.0D1*B)
      GA(KA+1)=-2.0D0*B-2.0D1*EXP(2.0D1*B)
      RETURN
  200 A=X(KA+1)**2
      B=X(KA)**2
      C=A+1.0D0
      D=B+1.0D0
      P=0.0D0
      IF (A.GT.P) P=LOG(A)
      Q=0.0D0
      IF (B.GT.Q) Q=LOG(B)
      IF (X(KA).EQ.0.0D0) THEN
        GA(KA)=0.0D0
      ELSE
        GA(KA)=2.0D0*X(KA)*(C*B**A+P*A**D)
      ENDIF
      IF (X(KA+1).EQ.0.0D0) THEN
        GA(KA+1)=0.0D0
      ELSE
        GA(KA+1)=2.0D0*X(KA+1)*(D*A**B+Q*B**C)
      ENDIF
      RETURN
  210 P=1.0D0/DBLE(N+1)
      Q=0.5D0*P**2
      A=2.0D0*X(KA)+Q*(X(KA)+DBLE(KA)*P+1.0D0)**3
      IF (KA.GT.1) A=A-X(KA-1)
      IF (KA.LT.N) A=A-X(KA+1)
      GA(KA)=A*(4.0D0+6.0D0*Q*(X(KA)+DBLE(KA)*P+1.0D0)**2.0D0)
      IF (KA.GT.1) GA(KA-1)=-2.0D0*A
      IF (KA.LT.N) GA(KA+1)=-2.0D0*A
      RETURN
  220 P=1.0D0/DBLE(N+1)
      Q=2.0D0/P
      R=2.0D0*P
      A=X(KA)-X(KA+1)
      GA(KA)=Q*(2.0D0*X(KA)-X(KA+1))
      GA(KA+1)=-Q*X(KA)
      IF (ABS(A).LE.1.0D-6) THEN
        GA(KA)=GA(KA)+R*EXP(X(KA+1))*(1.0D0/2.0D0+A*(1.0D0/3.0D0+A/
     &   8.0D0))
        GA(KA+1)=GA(KA+1)+R*EXP(X(KA+1))*(1.0D0/2.0D0+A*(1.0D0/6.0D0+A/
     &   24.0D0))
      ELSE
        B=EXP(X(KA))-EXP(X(KA+1))
        GA(KA)=GA(KA)+R*(EXP(X(KA))*A-B)/A**2
        GA(KA+1)=GA(KA+1)-R*(EXP(X(KA+1))*A-B)/A**2
      END IF
      IF (KA.EQ.1) THEN
        GA(1)=GA(1)+R*(EXP(X(1))*(X(1)-1.0D0)+1.0D0)/X(1)**2
      ELSE IF (KA.EQ.N-1) THEN
        GA(N)=GA(N)+2.0D0*Q*X(N)+R*(EXP(X(N))*(X(N)-1.0D0)+1.0D0)/X(N)**
     &   2
      END IF
      RETURN
  230 A=DBLE(KA)*SIN(X(KA))
      GA(KA)=A
      IF (KA.GT.1) GA(KA-1)=+DBLE(KA)*COS(X(KA-1))
      IF (KA.LT.N) GA(KA+1)=-DBLE(KA)*COS(X(KA+1))
      RETURN
  240 P=1.0D0/DBLE(N+1)
      IF (KA.EQ.1) THEN
        GA(KA)=0.5D0*X(KA)/P+P*EXP(X(KA))
        GA(KA+1)=0.25D0*X(KA+1)/P
      ELSE IF (KA.EQ.N) THEN
        GA(KA)=0.5D0*X(KA)/P+P*EXP(X(KA))
        GA(KA-1)=0.25D0*X(KA-1)/P
      ELSE
        A=0.25D0*(X(KA+1)-X(KA-1))/P
        GA(KA)=P*EXP(X(KA))
        GA(KA-1)=-A
        GA(KA+1)=A
      END IF
      RETURN
  250 P=1.0D0/DBLE(N+1)
      Q=DBLE(KA)*P
      IF (KA.EQ.1) THEN
        GA(KA)=X(KA)/P-2.0D0*P*(X(KA)+Q)
        GA(KA+1)=0.5D0*X(KA+1)/P
      ELSE IF (KA.EQ.N) THEN
        GA(KA)=X(KA)/P-2.0D0*P*(X(KA)+Q)
        GA(KA-1)=0.5D0*X(KA-1)/P
      ELSE
        A=0.5D0*(X(KA+1)-X(KA-1))/P
        GA(KA)=-2.0D0*P*(X(KA)+Q)
        GA(KA-1)=-A
        GA(KA+1)=A
      END IF
      RETURN
  260 P=1.0D0/DBLE(N+1)
      Q=EXP(2.0D0*DBLE(KA)*P)
      IF (KA.EQ.1) THEN
        R=1.0D0/3.0D0
        A=0.5D0*(X(KA+1)-R)/P
        GA(KA)=2.0D0*P*(X(KA)+Q)+(X(KA)-R)/P
        GA(KA+1)=A
      ELSE IF (KA.EQ.N) THEN
        R=EXP(2.0D0)/3.0D0
        A=0.5D0*(X(KA-1)-R)/P
        GA(KA)=2.0D0*P*(X(KA)+Q)+(X(KA)-R)/P
        GA(KA-1)=A
      ELSE
        A=0.5D0*(X(KA+1)-X(KA-1))/P
        GA(KA)=2.0D0*P*(X(KA)+Q)
        GA(KA-1)=-A
        GA(KA+1)=A
      END IF
      RETURN
  270 P=1.0D0/DBLE(N+1)
      A=EXP(-2.0D0*X(KA)**2)
      IF (KA.EQ.1) THEN
        B=0.5D0*X(KA+1)/P
        GA(KA)=X(KA)/P-4.0D0*X(KA)*A*P*(B**2-1.0D0)
        GA(KA+1)=A*B
      ELSE IF (KA.EQ.N) THEN
        B=0.5D0*X(KA-1)/P
        GA(KA)=X(KA)/P*EXP(-2.0D0)-4.0D0*X(KA)*A*P*(B**2-1.0D0)
        GA(KA-1)=A*B
      ELSE
        B=0.5D0*(X(KA+1)-X(KA-1))/P
        GA(KA)=-4.0D0*X(KA)*A*P*(B**2-1.0D0)
        GA(KA-1)=-A*B
        GA(KA+1)=A*B
      END IF
      RETURN
  280 P=1.0D0/DBLE(N+1)
      IF (KA.EQ.1) THEN
        A=0.5D0*(X(KA+1)-1.0D0)/P
        B=(X(KA)-1.0D0)/P
        U=0.5D0*ATAN(A)
        V=0.5D0*ATAN(B)
        GA(KA)=2.0D0*P*X(KA)+V
        GA(KA+1)=U
      ELSE IF (KA.EQ.N) THEN
        A=0.5D0*(2.0D0-X(KA-1))/P
        B=(2.0D0-X(KA))/P
        U=0.5D0*ATAN(A)
        V=0.5D0*ATAN(B)
        GA(KA)=2.0D0*P*X(KA)-V
        GA(KA-1)=-U
      ELSE
        A=0.5D0*(X(KA+1)-X(KA-1))/P
        U=0.5D0*ATAN(A)
        GA(KA)=2.0D0*P*X(KA)
        GA(KA-1)=-U
        GA(KA+1)=U
      END IF
      RETURN
  290 P=1.0D0/DBLE(N+1)
      IF (KA.EQ.1) THEN
        A=0.5D0*X(KA+1)/P
        B=X(KA)/P
        GA(KA)=2.0D2*P*(X(KA)-A**2)+2.0D2*B**3-(1.0D0-B)
        GA(KA+1)=-2.0D2*(X(KA)-A**2)*A-(1.0D0-A)
      ELSE IF (KA.EQ.N) THEN
        A=-0.5D0*X(KA-1)/P
        B=-X(KA)/P
        GA(KA)=2.0D2*P*(X(KA)-A**2)-2.0D2*B**3+(1.0D0-B)
        GA(KA-1)=2.0D2*(X(KA)-A**2)*A+(1.0D0-A)
      ELSE
        A=0.5D0*(X(KA+1)-X(KA-1))/P
        GA(KA)=2.0D2*P*(X(KA)-A**2)
        GA(KA-1)=2.0D2*(X(KA)-A**2)*A+(1.0D0-A)
        GA(KA+1)=-2.0D2*(X(KA)-A**2)*A-(1.0D0-A)
      END IF
      RETURN
      END
! SUBROUTINE TFFU14                ALL SYSTEMS                98/12/01
! PORTABILITY : ALL SYSTEMS
! 98/12/01 TU : ORIGINAL VERSION
!
! PURPOSE :
!  VALUES OF MODEL FUNCTIONS FOR UNCONSTRAINED MINIMIZATION.
!  SPARSE VERSION.
!
! PARAMETERS :
!  II  N  NUMBER OF VARIABLES.
!  RI  X(N)  VECTOR OF VARIABLES.
!  RO  F  VALUE OF THE MODEL FUNCTION.
!  II  NEXT  NUMBER OF THE TEST PROBLEM.
!
      SUBROUTINE TFFU14 (N, X, F, NEXT)
      INTEGER N,NEXT
      DOUBLE PRECISION X(*),F
      DOUBLE PRECISION A,B,C,D,P,Q,R,U,V
      INTEGER I,J,K
      F=0.0D0
      GO TO (10,30,50,70,90,110,140,160,190,220,250,280,300,320,340,360,
     &380,400,420,440,460,480),NEXT
   10 DO 20 J=2,N
        A=X(J-1)**2-X(J)
        B=X(J-1)-1.0D0
        F=F+1.0D2*A**2+B**2
   20 CONTINUE
      RETURN
   30 DO 40 J=2,N-2,2
        A=X(J-1)**2-X(J)
        B=X(J-1)-1.0D0
        C=X(J+1)**2-X(J+2)
        D=X(J+1)-1.0D0
        U=X(J)+X(J+2)-2.0D0
        V=X(J)-X(J+2)
        F=F+1.0D2*A**2+B**2+9.0D1*C**2+D**2+1.0D1*U**2+0.1D0*V**2
   40 CONTINUE
      RETURN
   50 DO 60 J=2,N-2,2
        A=X(J-1)+1.0D1*X(J)
        B=X(J+1)-X(J+2)
        C=X(J)-2.0D0*X(J+1)
        D=X(J-1)-X(J+2)
        F=F+A**2+5.0D0*B**2+C**4+1.0D1*D**4
   60 CONTINUE
      RETURN
   70 DO 80 J=2,N-2,2
        A=EXP(X(J-1))
        B=A-X(J)
        D=X(J)-X(J+1)
        P=X(J+1)-X(J+2)
        Q=SIN(P)/COS(P)
        U=X(J-1)
        V=X(J+2)-1.0D0
        F=F+B**4+1.0D2*D**6+Q**4+U**8+V**2
   80 CONTINUE
      RETURN
   90 P=7.0D0/3.0D0
      DO 100 J=1,N
        A=(3.0D0-2.0D0*X(J))*X(J)+1.0D0
        IF (J.GT.1) A=A-X(J-1)
        IF (J.LT.N) A=A-X(J+1)
        F=F+ABS(A)**P
  100 CONTINUE
      RETURN
  110 P=7.0D0/3.0D0
      DO 130 J=1,N
        A=(2.0D0+5.0D0*X(J)**2)*X(J)+1.0D0
        DO 120 I=MAX(1,J-5),MIN(N,J+1)
          A=A+X(I)*(1.0D0+X(I))
  120   CONTINUE
        F=F+ABS(A)**P
  130 CONTINUE
      RETURN
  140 P=7.0D0/3.0D0
      K=N/2
      DO 150 J=1,N
        A=(3.0D0-2.0D0*X(J))*X(J)+1.0D0
        IF (J.GT.1) A=A-X(J-1)
        IF (J.LT.N) A=A-X(J+1)
        F=F+ABS(A)**P
        IF (J.LE.K) THEN
          F=F+ABS(X(J)+X(J+K))**P
        END IF
  150 CONTINUE
      RETURN
  160 K=N/2
      DO 180 J=1,N
        P=0.0D0
        DO 170 I=J-2,J+2
          IF (I.LT.1.OR.I.GT.N) GO TO 170
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=DBLE(I+J)/1.0D1
          P=P+A*SIN(X(I))+B*COS(X(I))
  170   CONTINUE
        IF (J.GT.K) THEN
          I=J-K
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=DBLE(I+J)/1.0D1
          P=P+A*SIN(X(I))+B*COS(X(I))
        ELSE
          I=J+K
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=DBLE(I+J)/1.0D1
          P=P+A*SIN(X(I))+B*COS(X(I))
        END IF
        F=F+(DBLE(N+J)-P)**2/DBLE(N)
  180 CONTINUE
      RETURN
  190 K=N/2
      DO 210 J=1,N
        P=0.0D0
        DO 200 I=J-2,J+2
          IF (I.LT.1.OR.I.GT.N) GO TO 200
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=DBLE(I+J)/1.0D1
          P=P+A*SIN(X(I))+B*COS(X(I))
  200   CONTINUE
        IF (J.GT.K) THEN
          I=J-K
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=DBLE(I+J)/1.0D1
          P=P+A*SIN(X(I))+B*COS(X(I))
        ELSE
          I=J+K
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=DBLE(I+J)/1.0D1
          P=P+A*SIN(X(I))+B*COS(X(I))
        END IF
        F=F+(P+DBLE(J)*(1.0D0-COS(X(J))))/DBLE(N)
  210 CONTINUE
      RETURN
  220 K=N/2
      DO 240 J=1,N
        P=0.0D0
        Q=1.0D0+DBLE(J)/1.0D1
        DO 230 I=J-2,J+2
          IF (I.LT.1.OR.I.GT.N) GO TO 230
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=1.0D0+DBLE(I)/1.0D1
          C=DBLE(I+J)/1.0D1
          P=P+A*SIN(Q*X(J)+B*X(I)+C)
  230   CONTINUE
        IF (J.GT.K) THEN
          I=J-K
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=1.0D0+DBLE(I)/1.0D1
          C=DBLE(I+J)/1.0D1
          P=P+A*SIN(Q*X(J)+B*X(I)+C)
        ELSE
          I=J+K
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=1.0D0+DBLE(I)/1.0D1
          C=DBLE(I+J)/1.0D1
          P=P+A*SIN(Q*X(J)+B*X(I)+C)
        END IF
        F=F+P
  240 CONTINUE
      F=F/DBLE(N)
      RETURN
  250 P=-0.2008D-2
      Q=-0.1900D-2
      R=-0.0261D-2
      DO 270 I=0,N-5,5
        A=1.0D0
        B=0.0D0
        DO 260 J=1,5
          A=A*X(I+J)
          B=B+X(I+J)**2
  260   CONTINUE
        A=EXP(A)
        B=B-1.0D1-P
        C=X(I+2)*X(I+3)-5.0D0*X(I+4)*X(I+5)-Q
        D=X(I+1)**3+X(I+2)**3+1.0D0-R
        F=F+A+1.0D1*(B**2+C**2+D**2)
  270 CONTINUE
      RETURN
  280 DO 290 J=2,N
        A=X(J-1)-3.0D0
        B=X(J-1)-X(J)
        F=F+A**2+B**2+EXP(2.0D1*B)
  290 CONTINUE
      RETURN
  300 DO 310 J=2,N
        A=X(J)**2
        B=X(J-1)**2
        C=A+1.0D0
        D=B+1.0D0
        F=F+B**C+A**D
  310 CONTINUE
      RETURN
  320 P=1.0D0/DBLE(N+1)
      Q=0.5D0*P**2
      DO 330 J=1,N
        A=2.0D0*X(J)+Q*(X(J)+DBLE(J)*P+1.0D0)**3
        IF (J.GT.1) A=A-X(J-1)
        IF (J.LT.N) A=A-X(J+1)
        F=F+A**2
  330 CONTINUE
      RETURN
  340 P=1.0D0/DBLE(N+1)
      Q=2.0D0/P
      R=2.0D0*P
      DO 350 J=2,N
        A=X(J-1)-X(J)
        F=F+Q*X(J-1)*A
        IF (ABS(A).LE.1.0D-6) THEN
          F=F+R*EXP(X(J))*(1.0D0+A/2.0D0*(1.0D0+A/3.0D0*(1.0D0+A/4.0D0))
     &     )
        ELSE
          B=EXP(X(J-1))-EXP(X(J))
          F=F+R*B/A
        END IF
  350 CONTINUE
      F=F+Q*X(N)**2+R*(EXP(X(1))-1.0D0)/X(1)+R*(EXP(X(N))-1.0D0)/X(N)
      RETURN
  360 DO 370 J=1,N
        A=DBLE(J)*(1.0D0-COS(X(J)))
        IF (J.GT.1) A=A+DBLE(J)*SIN(X(J-1))
        IF (J.LT.N) A=A-DBLE(J)*SIN(X(J+1))
        F=F+A
  370 CONTINUE
      RETURN
  380 P=1.0D0/DBLE(N+1)
      DO 390 J=1,N
        IF (J.EQ.1) THEN
          F=F+0.25D0*X(J)**2/P+1.25D-1*X(J+1)**2/P+P*(EXP(X(J))-1.0D0)
        ELSE IF (J.EQ.N) THEN
          F=F+0.25D0*X(J)**2/P+1.25D-1*X(J-1)**2/P+P*(EXP(X(J))-1.0D0)
        ELSE
          F=F+1.25D-1*(X(J+1)-X(J-1))**2/P+P*(EXP(X(J))-1.0D0)
        END IF
  390 CONTINUE
      RETURN
  400 P=1.0D0/DBLE(N+1)
      DO 410 J=1,N
        Q=DBLE(J)*P
        IF (J.EQ.1) THEN
          F=F+0.5D0*X(J)**2/P+0.25D0*X(J+1)**2/P-P*(X(J)**2+2.0D0*X(J)*
     &     Q)
        ELSE IF (J.EQ.N) THEN
          F=F+0.5D0*X(J)**2/P+0.25D0*X(J-1)**2/P-P*(X(J)**2+2.0D0*X(J)*
     &     Q)
        ELSE
          F=F+2.5D-1*(X(J+1)-X(J-1))**2/P-P*(X(J)**2+2.0D0*X(J)*Q)
        END IF
  410 CONTINUE
      RETURN
  420 P=1.0D0/DBLE(N+1)
      DO 430 J=1,N
        Q=EXP(2.0D0*DBLE(J)*P)
        IF (J.EQ.1) THEN
          R=1.0D0/3.0D0
          F=F+0.5D0*(X(J)-R)**2/P+7.0D0*R**2+2.5D-1*(X(J+1)-R)**2/P+P*
     &     (X(J)**2+2.0D0*X(J)*Q)
        ELSE IF (J.EQ.N) THEN
          R=EXP(2.0D0)/3.0D0
          F=F+0.5D0*(X(J)-R)**2/P+7.0D0*R**2+2.5D-1*(X(J-1)-R)**2/P+P*
     &     (X(J)**2+2.0D0*X(J)*Q)
        ELSE
          F=F+2.5D-1*(X(J+1)-X(J-1))**2/P+P*(X(J)**2+2.0D0*X(J)*Q)
        END IF
  430 CONTINUE
      RETURN
  440 P=1.0D0/DBLE(N+1)
      DO 450 J=1,N
        IF (J.EQ.1) THEN
          F=F+(0.5D0*X(J)**2/P-P)+(2.5D-1*X(J+1)**2/P-P)*EXP(-2.0D0*X(J)
     &     **2)
        ELSE IF (J.EQ.N) THEN
          F=F+(0.5D0*X(J)**2/P-P)*EXP(-2.0D0)+(2.5D-1*X(J-1)**2/P-P)*
     &     EXP(-2.0D0*X(J)**2)
        ELSE
          F=F+(2.5D-1*(X(J+1)-X(J-1))**2/P-P)*EXP(-2.0D0*X(J)**2)
        END IF
  450 CONTINUE
      RETURN
  460 P=1.0D0/DBLE(N+1)
      DO 470 J=1,N
        IF (J.EQ.1) THEN
          A=0.5D0*(X(J+1)-1.0D0)/P
          B=(X(J)-1.0D0)/P
          F=F+P*(X(J)**2+A*ATAN(A)-LOG(SQRT(1.0D0+A**2)))+0.5D0*P*
     &     (1.0D0+B*ATAN(B)-LOG(SQRT(1.0D0+B**2)))
        ELSE IF (J.EQ.N) THEN
          A=0.5D0*(2.0D0-X(J-1))/P
          B=(2.0D0-X(J))/P
          F=F+P*(X(J)**2+A*ATAN(A)-LOG(SQRT(1.0D0+A**2)))+0.5D0*P*
     &     (4.0D0+B*ATAN(B)-LOG(SQRT(1.0D0+B**2)))
        ELSE
          A=0.5D0*(X(J+1)-X(J-1))/P
          F=F+P*(X(J)**2+A*ATAN(A)-LOG(SQRT(1.0D0+A**2)))
        END IF
  470 CONTINUE
      RETURN
  480 P=1.0D0/DBLE(N+1)
      DO 490 J=1,N
        IF (J.EQ.1) THEN
          A=0.5D0*X(J+1)/P
          B=X(J)/P
          F=F+P*(1.0D2*(X(J)-A**2)**2+(1.0D0-A)**2)+0.5D0*P*(1.0D2*B**4+
     &     (1.0D0-B)**2)
        ELSE IF (J.EQ.N) THEN
          A=-0.5D0*X(J-1)/P
          B=-X(J)/P
          F=F+P*(1.0D2*(X(J)-A**2)**2+(1.0D0-A)**2)+0.5D0*P*(1.0D2*B**4+
     &     (1.0D0-B)**2)
        ELSE
          A=0.5D0*(X(J+1)-X(J-1))/P
          F=F+P*(1.0D2*(X(J)-A**2)**2+(1.0D0-A)**2)
        END IF
  490 CONTINUE
      RETURN
      END
! SUBROUTINE TFGU14                ALL SYSTEMS                98/12/01
! PORTABILITY : ALL SYSTEMS
! 98/12/01 TU : ORIGINAL VERSION
!
! PURPOSE :
!  GRADIENTS OF MODEL FUNCTIONS FOR UNCONSTRAINED MINIMIZATION.
!  SPARSE VERSION - LU VERSION WITH MODIFIED TESTS NO 7-10,12.
!
! PARAMETERS :
!  II  N  NUMBER OF VARIABLES.
!  RI  X(N)  VECTOR OF VARIABLES.
!  RO  G(N)  GRADIENT OF THE MODEL FUNCTION.
!  II  NEXT  NUMBER OF THE TEST PROBLEM.
!
      SUBROUTINE TFGU14 (N, X, G, NEXT)
      INTEGER N,NEXT
      DOUBLE PRECISION X(*),G(*)
      DOUBLE PRECISION A,B,C,D,P,Q,R,U,V
      INTEGER I,J,K
      DO 10 I=1,N
        G(I)=0.0D0
   10 CONTINUE
      GO TO (20,40,60,80,100,120,160,180,220,250,280,310,330,350,370,
     &390,410,430,450,470,490,510),NEXT
   20 DO 30 J=2,N
        A=X(J-1)**2-X(J)
        B=X(J-1)-1.0D0
        G(J-1)=G(J-1)+4.0D2*X(J-1)*A+2.0D0*B
        G(J)=G(J)-2.0D2*A
   30 CONTINUE
      RETURN
   40 DO 50 J=2,N-2,2
        A=X(J-1)**2-X(J)
        B=X(J-1)-1.0D0
        C=X(J+1)**2-X(J+2)
        D=X(J+1)-1.0D0
        U=X(J)+X(J+2)-2.0D0
        V=X(J)-X(J+2)
        G(J-1)=G(J-1)+4.0D2*X(J-1)*A+2.0D0*B
        G(J)=G(J)-2.0D2*A+2.0D1*U+0.2D0*V
        G(J+1)=G(J+1)+3.6D2*X(J+1)*C+2.0D0*D
        G(J+2)=G(J+2)-1.8D2*C+2.0D1*U-0.2D0*V
   50 CONTINUE
      RETURN
   60 DO 70 J=2,N-2,2
        A=X(J-1)+1.0D1*X(J)
        B=X(J+1)-X(J+2)
        C=X(J)-2.0D0*X(J+1)
        D=X(J-1)-X(J+2)
        G(J-1)=G(J-1)+2.0D0*A+4.0D1*D**3
        G(J)=G(J)+2.0D1*A+4.0D0*C**3
        G(J+1)=G(J+1)-8.0D0*C**3+1.0D1*B
        G(J+2)=G(J+2)-4.0D1*D**3-1.0D1*B
   70 CONTINUE
      RETURN
   80 DO 90 J=2,N-2,2
        A=EXP(X(J-1))
        B=A-X(J)
        B=4.0D0*B**3
        D=X(J)-X(J+1)
        D=6.0D2*D**5
        P=X(J+1)-X(J+2)
        C=COS(P)
        Q=SIN(P)/COS(P)
        Q=4.0D0*Q**3/C**2
        U=X(J-1)
        V=X(J+2)-1.0D0
        G(J-1)=G(J-1)+A*B+8.0D0*U**7
        G(J)=G(J)+D-B
        G(J+1)=G(J+1)+Q-D
        G(J+2)=G(J+2)+2.0D0*V-Q
   90 CONTINUE
      RETURN
  100 P=7.0D0/3.0D0
      DO 110 J=1,N
        A=(3.0D0-2.0D0*X(J))*X(J)+1.0D0
        IF (J.GT.1) A=A-X(J-1)
        IF (J.LT.N) A=A-X(J+1)
        B=P*ABS(A)**(P-1.0D0)*SIGN(1.0D0,A)
        G(J)=G(J)+B*(3.0D0-4.0D0*X(J))
        IF (J.GT.1) G(J-1)=G(J-1)-B
        IF (J.LT.N) G(J+1)=G(J+1)-B
  110 CONTINUE
      RETURN
  120 P=7.0D0/3.0D0
      DO 150 J=1,N
        A=(2.0D0+5.0D0*X(J)**2)*X(J)+1.0D0
        DO 130 I=MAX(1,J-5),MIN(N,J+1)
          A=A+X(I)*(1.0D0+X(I))
  130   CONTINUE
        B=P*ABS(A)**(P-1.0D0)*SIGN(1.0D0,A)
        G(J)=G(J)+B*(2.0D0+1.5D1*X(J)**2)
        DO 140 I=MAX(1,J-5),MIN(N,J+1)
          G(I)=G(I)+B*(1.0D0+2.0D0*X(I))
  140   CONTINUE
  150 CONTINUE
      RETURN
  160 P=7.0D0/3.0D0
      K=N/2
      DO 170 J=1,N
        A=(3.0D0-2.0D0*X(J))*X(J)+1.0D0
        IF (J.GT.1) A=A-X(J-1)
        IF (J.LT.N) A=A-X(J+1)
        B=P*ABS(A)**(P-1.0D0)*SIGN(1.0D0,A)
        G(J)=G(J)+B*(3.0D0-4.0D0*X(J))
        IF (J.GT.1) G(J-1)=G(J-1)-B
        IF (J.LT.N) G(J+1)=G(J+1)-B
        IF (J.LE.K) THEN
          A=X(J)+X(J+K)
          B=P*ABS(A)**(P-1.0D0)*SIGN(1.0D0,A)
          G(J)=G(J)+B
          G(J+K)=G(J+K)+B
        END IF
  170 CONTINUE
      RETURN
  180 K=N/2
      DO 210 J=1,N
        P=0.0D0
        DO 190 I=J-2,J+2
          IF (I.LT.1.OR.I.GT.N) GO TO 190
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=DBLE(I+J)/1.0D1
          P=P+A*SIN(X(I))+B*COS(X(I))
  190   CONTINUE
        IF (J.GT.K) THEN
          I=J-K
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=DBLE(I+J)/1.0D1
          P=P+A*SIN(X(I))+B*COS(X(I))
        ELSE
          I=J+K
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=DBLE(I+J)/1.0D1
          P=P+A*SIN(X(I))+B*COS(X(I))
        END IF
        P=2.0D0*(DBLE(N+J)-P)/DBLE(N)
        DO 200 I=J-2,J+2
          IF (I.LT.1.OR.I.GT.N) GO TO 200
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=DBLE(I+J)/1.0D1
          G(I)=G(I)-P*(A*COS(X(I))-B*SIN(X(I)))
  200   CONTINUE
        IF (J.GT.K) THEN
          I=J-K
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=DBLE(I+J)/1.0D1
          G(I)=G(I)-P*(A*COS(X(I))-B*SIN(X(I)))
        ELSE
          I=J+K
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=DBLE(I+J)/1.0D1
          G(I)=G(I)-P*(A*COS(X(I))-B*SIN(X(I)))
        END IF
  210 CONTINUE
      RETURN
  220 K=N/2
      P=1.0D0/DBLE(N)
      DO 240 J=1,N
        DO 230 I=J-2,J+2
          IF (I.LT.1.OR.I.GT.N) GO TO 230
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=DBLE(I+J)/1.0D1
          G(I)=G(I)+P*(A*COS(X(I))-B*SIN(X(I)))
  230   CONTINUE
        IF (J.GT.K) THEN
          I=J-K
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=DBLE(I+J)/1.0D1
          G(I)=G(I)+P*(A*COS(X(I))-B*SIN(X(I)))
        ELSE
          I=J+K
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=DBLE(I+J)/1.0D1
          G(I)=G(I)+P*(A*COS(X(I))-B*SIN(X(I)))
        END IF
        G(J)=G(J)+P*DBLE(J)*SIN(X(J))
  240 CONTINUE
      RETURN
  250 K=N/2
      DO 270 J=1,N
        Q=1.0D0+DBLE(J)/1.0D1
        DO 260 I=J-2,J+2
          IF (I.LT.1.OR.I.GT.N) GO TO 260
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=1.0D0+DBLE(I)/1.0D1
          C=DBLE(I+J)/1.0D1
          P=A*COS(Q*X(J)+B*X(I)+C)/DBLE(N)
          G(J)=G(J)+P*Q
          G(I)=G(I)+P*B
  260   CONTINUE
        IF (J.GT.K) THEN
          I=J-K
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=1.0D0+DBLE(I)/1.0D1
          C=DBLE(I+J)/1.0D1
          P=A*COS(Q*X(J)+B*X(I)+C)/DBLE(N)
          G(J)=G(J)+P*Q
          G(I)=G(I)+P*B
        ELSE
          I=J+K
          A=5.0D0*(1.0D0+MOD(I,5)+MOD(J,5))
          B=1.0D0+DBLE(I)/1.0D1
          C=DBLE(I+J)/1.0D1
          P=A*COS(Q*X(J)+B*X(I)+C)/DBLE(N)
          G(J)=G(J)+P*Q
          G(I)=G(I)+P*B
        END IF
  270 CONTINUE
      RETURN
  280 P=-0.2008D-2
      Q=-0.1900D-2
      R=-0.0261D-2
      DO 300 I=0,N-5,5
        A=1.0D0
        B=0.0D0
        DO 290 J=1,5
          A=A*X(I+J)
          B=B+X(I+J)**2
  290   CONTINUE
        A=A*EXP(A)
        B=B-1.0D1-P
        C=X(I+2)*X(I+3)-5.0D0*X(I+4)*X(I+5)-Q
        D=X(I+1)**3+X(I+2)**3+1.0D0-R
        G(I+1)=G(I+1)+A/X(I+1)+2.0D1*(2.0D0*B*X(I+1)+3.0D0*D*X(I+1)**2)
        G(I+2)=G(I+2)+A/X(I+2)+2.0D1*(2.0D0*B*X(I+2)+C*X(I+3)+3.0D0*D*
     &   X(I+2)**2)
        G(I+3)=G(I+3)+A/X(I+3)+2.0D1*(2.0D0*B*X(I+3)+C*X(I+2))
        G(I+4)=G(I+4)+A/X(I+4)+2.0D1*(2.0D0*B*X(I+4)-5.0D0*C*X(I+5))
        G(I+5)=G(I+5)+A/X(I+5)+2.0D1*(2.0D0*B*X(I+5)-5.0D0*C*X(I+4))
  300 CONTINUE
      RETURN
  310 DO 320 J=2,N
        A=X(J-1)-3.0D0
        B=X(J-1)-X(J)
        G(J-1)=G(J-1)+2.0D0*A+2.0D0*B+2.0D1*EXP(2.0D1*B)
        G(J)=G(J)-2.0D0*B-2.0D1*EXP(2.0D1*B)
  320 CONTINUE
      RETURN
  330 DO 340 J=2,N
        A=X(J)**2
        B=X(J-1)**2
        C=A+1.0D0
        D=B+1.0D0
        P=0.0D0
        IF (A.GT.P) P=LOG(A)
        Q=0.0D0
        IF (B.GT.Q) Q=LOG(B)
        G(J-1)=G(J-1)+2.0D0*X(J-1)*(C*B**A+P*A**D)
        G(J)=G(J)+2.0D0*X(J)*(D*A**B+Q*B**C)
  340 CONTINUE
      RETURN
  350 P=1.0D0/DBLE(N+1)
      Q=0.5D0*P**2
      DO 360 J=1,N
        A=2.0D0*X(J)+Q*(X(J)+DBLE(J)*P+1.0D0)**3
        IF (J.GT.1) A=A-X(J-1)
        IF (J.LT.N) A=A-X(J+1)
        G(J)=G(J)+A*(4.0D0+6.0D0*Q*(X(J)+DBLE(J)*P+1.0D0)**2.0D0)
        IF (J.GT.1) G(J-1)=G(J-1)-2.0D0*A
        IF (J.LT.N) G(J+1)=G(J+1)-2.0D0*A
  360 CONTINUE
      RETURN
  370 P=1.0D0/DBLE(N+1)
      Q=2.0D0/P
      R=2.0D0*P
      DO 380 J=2,N
        A=X(J-1)-X(J)
        G(J-1)=G(J-1)+Q*(2.0D0*X(J-1)-X(J))
        G(J)=G(J)-Q*X(J-1)
        IF (ABS(A).LE.1.0D-6) THEN
          G(J-1)=G(J-1)+R*EXP(X(J))*(1.0D0/2.0D0+A*(1.0D0/3.0D0+A/8.0D0)
     &     )
          G(J)=G(J)+R*EXP(X(J))*(1.0D0/2.0D0+A*(1.0D0/6.0D0+A/24.0D0))
        ELSE
          B=EXP(X(J-1))-EXP(X(J))
          G(J-1)=G(J-1)+R*(EXP(X(J-1))*A-B)/A**2
          G(J)=G(J)-R*(EXP(X(J))*A-B)/A**2
        END IF
  380 CONTINUE
      G(1)=G(1)+R*(EXP(X(1))*(X(1)-1.0D0)+1.0D0)/X(1)**2
      G(N)=G(N)+2.0D0*Q*X(N)+R*(EXP(X(N))*(X(N)-1.0D0)+1.0D0)/X(N)**2
      RETURN
  390 DO 400 J=1,N
        A=DBLE(J)*SIN(X(J))
        G(J)=G(J)+A
        IF (J.GT.1) G(J-1)=G(J-1)+DBLE(J)*COS(X(J-1))
        IF (J.LT.N) G(J+1)=G(J+1)-DBLE(J)*COS(X(J+1))
  400 CONTINUE
      RETURN
  410 P=1.0D0/DBLE(N+1)
      DO 420 J=1,N
        IF (J.EQ.1) THEN
          G(J)=G(J)+0.5D0*X(J)/P+P*EXP(X(J))
          G(J+1)=G(J+1)+0.25D0*X(J+1)/P
        ELSE IF (J.EQ.N) THEN
          G(J)=G(J)+0.5D0*X(J)/P+P*EXP(X(J))
          G(J-1)=G(J-1)+0.25D0*X(J-1)/P
        ELSE
          A=0.25D0*(X(J+1)-X(J-1))/P
          G(J)=G(J)+P*EXP(X(J))
          G(J-1)=G(J-1)-A
          G(J+1)=G(J+1)+A
        END IF
  420 CONTINUE
      RETURN
  430 P=1.0D0/DBLE(N+1)
      DO 440 J=1,N
        Q=DBLE(J)*P
        IF (J.EQ.1) THEN
          G(J)=G(J)+X(J)/P-2.0D0*P*(X(J)+Q)
          G(J+1)=G(J+1)+0.5D0*X(J+1)/P
        ELSE IF (J.EQ.N) THEN
          G(J)=G(J)+X(J)/P-2.0D0*P*(X(J)+Q)
          G(J-1)=G(J-1)+0.5D0*X(J-1)/P
        ELSE
          A=0.5D0*(X(J+1)-X(J-1))/P
          G(J)=G(J)-2.0D0*P*(X(J)+Q)
          G(J-1)=G(J-1)-A
          G(J+1)=G(J+1)+A
        END IF
  440 CONTINUE
      RETURN
  450 P=1.0D0/DBLE(N+1)
      DO 460 J=1,N
        Q=EXP(2.0D0*DBLE(J)*P)
        IF (J.EQ.1) THEN
          R=1.0D0/3.0D0
          A=0.5D0*(X(J+1)-R)/P
          G(J)=G(J)+2.0D0*P*(X(J)+Q)+(X(J)-R)/P
          G(J+1)=G(J+1)+A
        ELSE IF (J.EQ.N) THEN
          R=EXP(2.0D0)/3.0D0
          A=0.5D0*(X(J-1)-R)/P
          G(J)=G(J)+2.0D0*P*(X(J)+Q)+(X(J)-R)/P
          G(J-1)=G(J-1)+A
        ELSE
          A=0.5D0*(X(J+1)-X(J-1))/P
          G(J)=G(J)+2.0D0*P*(X(J)+Q)
          G(J-1)=G(J-1)-A
          G(J+1)=G(J+1)+A
        END IF
  460 CONTINUE
      RETURN
  470 P=1.0D0/DBLE(N+1)
      DO 480 J=1,N
        A=EXP(-2.0D0*X(J)**2)
        IF (J.EQ.1) THEN
          B=0.5D0*X(J+1)/P
          G(J)=G(J)+X(J)/P-4.0D0*X(J)*A*P*(B**2-1.0D0)
          G(J+1)=G(J+1)+A*B
        ELSE IF (J.EQ.N) THEN
          B=0.5D0*X(J-1)/P
          G(J)=G(J)+X(J)/P*EXP(-2.0D0)-4.0D0*X(J)*A*P*(B**2-1.0D0)
          G(J-1)=G(J-1)+A*B
        ELSE
          B=0.5D0*(X(J+1)-X(J-1))/P
          G(J)=G(J)-4.0D0*X(J)*A*P*(B**2-1.0D0)
          G(J-1)=G(J-1)-A*B
          G(J+1)=G(J+1)+A*B
        END IF
  480 CONTINUE
      RETURN
  490 P=1.0D0/DBLE(N+1)
      DO 500 J=1,N
        IF (J.EQ.1) THEN
          A=0.5D0*(X(J+1)-1.0D0)/P
          B=(X(J)-1.0D0)/P
          U=0.5D0*ATAN(A)
          V=0.5D0*ATAN(B)
          G(J)=G(J)+2.0D0*P*X(J)+V
          G(J+1)=G(J+1)+U
        ELSE IF (J.EQ.N) THEN
          A=0.5D0*(2.0D0-X(J-1))/P
          B=(2.0D0-X(J))/P
          U=0.5D0*ATAN(A)
          V=0.5D0*ATAN(B)
          G(J)=G(J)+2.0D0*P*X(J)-V
          G(J-1)=G(J-1)-U
        ELSE
          A=0.5D0*(X(J+1)-X(J-1))/P
          U=0.5D0*ATAN(A)
          G(J)=G(J)+2.0D0*P*X(J)
          G(J-1)=G(J-1)-U
          G(J+1)=G(J+1)+U
        END IF
  500 CONTINUE
      RETURN
  510 P=1.0D0/DBLE(N+1)
      DO 520 J=1,N
        IF (J.EQ.1) THEN
          A=0.5D0*X(J+1)/P
          B=X(J)/P
          G(J)=G(J)+2.0D2*P*(X(J)-A**2)+2.0D2*B**3-(1.0D0-B)
          G(J+1)=G(J+1)-2.0D2*(X(J)-A**2)*A-(1.0D0-A)
        ELSE IF (J.EQ.N) THEN
          A=-0.5D0*X(J-1)/P
          B=-X(J)/P
          G(J)=G(J)+2.0D2*P*(X(J)-A**2)-2.0D2*B**3+(1.0D0-B)
          G(J-1)=G(J-1)+2.0D2*(X(J)-A**2)*A+(1.0D0-A)
        ELSE
          A=0.5D0*(X(J+1)-X(J-1))/P
          G(J)=G(J)+2.0D2*P*(X(J)-A**2)
          G(J-1)=G(J-1)+2.0D2*(X(J)-A**2)*A+(1.0D0-A)
          G(J+1)=G(J+1)-2.0D2*(X(J)-A**2)*A-(1.0D0-A)
        END IF
  520 CONTINUE
      RETURN
      END
! SUBROUTINE TIUB15                ALL SYSTEMS                92/12/01
! PORTABILITY : ALL SYSTEMS
! 92/12/01 RA : ORIGINAL VERSION
!
! PURPOSE :
!  INITIATION OF VARIABLES AND DEFINITION OF STRUCTURE OF THE SPARSE
!  JACOBIAN MATRIX FOR THE SUM OF SQUARES.
!
! PARAMETERS :
!  IU  N  NUMBER OF VARIABLES.
!  IO  NA  NUMBER OF PARTIAL FUNCTIONS.
!  IO  MA  NUMBER OF NONZERO ELEMENTS IN THE SPARSE JACOBIAN MATRIX.
!  RO  X(N)  VECTOR OF VARIABLES.
!  IO  IA(NA+1)  POINTERS OF FIRST IN THE ROW ELEMENTS IN THE
!         SPARSE JACOBIAN MATRIX.
!  IO  JA(MA)  COLUMN INDICES OF NONZERO ELEMENTS IN THE SPARSE
!         JACOBIAN MATRIX.
!  RO  FMIN  LOWER BOUND FOR THE OBJECTIVE FUNCTION VALUE.
!  RO  XMAX  MAXIMUM ALLOWED STEPSIZE.
!  II  NEXT  NUMBER OF THE SELECTED TEST PROBLEM.
!  IO  IERR  ERROR INDICATOR. IERR=0 FOR CORRECT OUTPUT.
!
      SUBROUTINE TIUB15 (N, NA, MA, X, IA, JA, FMIN, XMAX, NEXT, IERR)
      INTEGER N,NA,MA,NEXT,IERR
      INTEGER IA(*),JA(*)
      DOUBLE PRECISION X(*),FMIN,XMAX
      INTEGER I,J,K,L,II,JJ,KK,LL,MM,KA
      DOUBLE PRECISION Y(20)
      COMMON /EMPR15/ Y
      DOUBLE PRECISION ETA9
      PARAMETER  (ETA9=1.0D60)
      FMIN=0.0D0
      XMAX=1.0D3
      IERR=0
      GO TO (10,40,70,100,130,160,230,250,280,320,350,380,410,440,500,
     &520,540,560,580,600,620,650),NEXT
   10 IF (N.LT.2) GO TO 680
      N=N-MOD(N,2)
      DO 20 I=1,N
        IF (MOD(I,2).EQ.1) THEN
          X(I)=-1.2D0
        ELSE
          X(I)=1.0D0
        END IF
   20 CONTINUE
      DO 30 I=1,N-1
        L=3*(I-1)+1
        JA(L)=I
        JA(L+1)=I+1
        JA(L+2)=I
        K=2*(I-1)+1
        IA(K)=L
        IA(K+1)=L+2
   30 CONTINUE
      IA(K+2)=L+3
      NA=2*(N-1)
      MA=3*(N-1)
      RETURN
   40 IF (N.LT.4) GO TO 680
      N=N-MOD(N,2)
      DO 50 I=1,N
        IF (MOD(I,2).EQ.1) THEN
          X(I)=-2.0D0
          IF (I.LE.4) X(I)=-3.0D0
        ELSE
          X(I)=0.0D0
          IF (I.LE.4) X(I)=-1.0D0
        END IF
   50 CONTINUE
      DO 60 I=2,N-2,2
        L=5*(I-2)+1
        JA(L)=I-1
        JA(L+1)=I
        JA(L+2)=I-1
        JA(L+3)=I+1
        JA(L+4)=I+2
        JA(L+5)=I+1
        JA(L+6)=I
        JA(L+7)=I+2
        JA(L+8)=I
        JA(L+9)=I+2
        K=3*(I-2)+1
        IA(K)=L
        IA(K+1)=L+2
        IA(K+2)=L+3
        IA(K+3)=L+5
        IA(K+4)=L+6
        IA(K+5)=L+8
   60 CONTINUE
      IA(K+6)=L+10
      NA=K+5
      MA=L+9
      RETURN
   70 IF (N.LT.4) GO TO 680
      N=N-MOD(N,2)
      DO 80 I=1,N
        IF (MOD(I,4).EQ.1) THEN
          X(I)=3.0D0
        ELSE IF (MOD(I,4).EQ.2) THEN
          X(I)=-1.0D0
        ELSE IF (MOD(I,4).EQ.3) THEN
          X(I)=0.0D0
        ELSE
          X(I)=1.0D0
        END IF
   80 CONTINUE
      DO 90 I=2,N-2,2
        L=4*(I-2)+1
        JA(L)=I-1
        JA(L+1)=I
        JA(L+2)=I+1
        JA(L+3)=I+2
        JA(L+4)=I
        JA(L+5)=I+1
        JA(L+6)=I-1
        JA(L+7)=I+2
        K=2*(I-2)+1
        IA(K)=L
        IA(K+1)=L+2
        IA(K+2)=L+4
        IA(K+3)=L+6
   90 CONTINUE
      IA(K+4)=L+8
      NA=K+3
      MA=L+7
      RETURN
  100 IF (N.LT.4) GO TO 680
      N=N-MOD(N,2)
      DO 110 I=1,N
        X(I)=2.0D0
  110 CONTINUE
      X(1)=1.0D0
      DO 120 I=2,N-2,2
        L=4*(I-2)+1
        JA(L)=I-1
        JA(L+1)=I
        JA(L+2)=I
        JA(L+3)=I+1
        JA(L+4)=I+1
        JA(L+5)=I+2
        JA(L+6)=I-1
        JA(L+7)=I+2
        K=5*(I-2)/2+1
        IA(K)=L
        IA(K+1)=L+2
        IA(K+2)=L+4
        IA(K+3)=L+6
        IA(K+4)=L+7
  120 CONTINUE
      IA(K+5)=L+8
      NA=K+4
      MA=L+7
      XMAX=1.0D1
      RETURN
  130 IF (N.LT.3) GO TO 680
      DO 140 I=1,N
        X(I)=-1.0D0
  140 CONTINUE
      JA(1)=1
      JA(2)=2
      IA(1)=1
      DO 150 I=2,N-1
        K=3*(I-2)+3
        JA(K)=I-1
        JA(K+1)=I
        JA(K+2)=I+1
        IA(I)=K
  150 CONTINUE
      JA(K+3)=N-1
      JA(K+4)=N
      IA(N)=K+3
      IA(N+1)=IA(N)+2
      NA=N
      MA=3*N-2
      RETURN
  160 IF (N.LT.6) GO TO 680
      DO 170 I=1,N
        X(I)=-1.0D0
  170 CONTINUE
      L=1
      DO 190 I=1,5
        IA(I)=L
        DO 180 K=1,I+1
          JA(L)=K
          L=L+1
  180   CONTINUE
  190 CONTINUE
      DO 210 I=6,N-1
        IA(I)=L
        DO 200 K=I-5,I+1
          JA(L)=K
          L=L+1
  200   CONTINUE
  210 CONTINUE
      IA(N)=L
      DO 220 K=N-5,N
        JA(L)=K
        L=L+1
  220 CONTINUE
      IA(N+1)=L
      NA=N
      MA=L-1
      RETURN
  230 IF (N.LT.2) GO TO 680
      DO 240 I=1,N-1
        X(I)=0.5D0
        K=4*(I-1)+1
        JA(K)=I
        JA(K+1)=I+1
        JA(K+2)=I
        JA(K+3)=I+1
        L=2*(I-1)+1
        IA(L)=K
        IA(L+1)=K+2
  240 CONTINUE
      X(N)=-2.0D0
      IA(L+2)=IA(L+1)+2
      NA=2*(N-1)
      MA=2*NA
      RETURN
  250 IF (N.LT.4) GO TO 680
      N=N-MOD(N,4)
      DO 260 I=1,N
        X(I)=SIN(DBLE(I))**2
  260 CONTINUE
      MM=5*N
      K=1
      DO 270 KA=1,MM
        I=MOD(KA,N/2)+1
        J=I+N/2
        K=2*(KA-1)+1
        JA(K)=I
        JA(K+1)=J
        IA(KA)=K
  270 CONTINUE
      IA(MM+1)=K+2
      NA=MM
      MA=2*MM
      RETURN
  280 IF (N.LT.4) GO TO 680
      N=N-MOD(N,2)
      DO 290 I=1,N
        X(I)=5.0D0
  290 CONTINUE
      DO 310 I=2,N-2,2
        L=12*(I-2)+1
        DO 300 J=1,6
          K=(J-1)*4
          JA(L+K)=I-1
          JA(L+K+1)=I
          JA(L+K+2)=I+1
          JA(L+K+3)=I+2
  300   CONTINUE
        K=3*(I-2)+1
        IA(K)=L
        IA(K+1)=L+4
        IA(K+2)=L+8
        IA(K+3)=L+12
        IA(K+4)=L+16
        IA(K+5)=L+20
  310 CONTINUE
      IA(K+6)=L+24
      NA=K+5
      MA=L+23
      RETURN
  320 IF (N.LT.2) GO TO 680
      DO 330 I=1,N
        X(I)=0.2D0
  330 CONTINUE
      MM=2*N-2
      JA(1)=1
      JA(2)=2
      JA(3)=1
      JA(4)=2
      IA(1)=1
      IA(2)=3
      K=5
      L=3
      DO 340 KA=3,MM,2
        I=(KA+1)/2-1
        JA(K)=I
        JA(K+1)=I+1
        JA(K+2)=I+2
        JA(K+3)=I+1
        JA(K+4)=I+2
        IA(L)=K
        IA(L+1)=K+3
        K=K+5
        L=L+2
  340 CONTINUE
      JA(K)=N-1
      JA(K+1)=N
      IA(L)=K
      IA(L+1)=K+2
      NA=L
      MA=K+1
      RETURN
  350 CONTINUE
      IF (N.LT.2) GO TO 680
      N=N-MOD(N,2)
      DO 360 I=1,N
        IF (MOD(I,2).EQ.1) THEN
          X(I)=-0.8D0
        ELSE
          X(I)=-0.8D0
        END IF
  360 CONTINUE
      IA(1)=1
      NA=2*(N-1)
      MA=1
      L=1
      DO 370 K=1,N-1
        JA(MA)=K
        MA=MA+1
        JA(MA)=K+1
        MA=MA+1
        L=L+1
        IA(L)=MA
        JA(MA)=K
        MA=MA+1
        L=L+1
        IA(L)=MA
  370 CONTINUE
      MA=IA(NA+1)-1
      RETURN
  380 CONTINUE
      IF (N.LT.5) GO TO 680
      IF (MOD(N-5,3).NE.0) N=N-MOD(N-5,3)
      DO 390 I=1,N
        X(I)=-1.0D0
  390 CONTINUE
      IA(1)=1
      KK=(N-5)/3+1
      NA=6*KK
      MA=1
      L=1
      DO 400 K=1,KK
        I=3*(K-1)+1
        JA(MA)=I
        MA=MA+1
        JA(MA)=I+1
        MA=MA+1
        L=L+1
        IA(L)=MA
        JA(MA)=I+2
        MA=MA+1
        L=L+1
        IA(L)=MA
        JA(MA)=I+3
        MA=MA+1
        L=L+1
        IA(L)=MA
        JA(MA)=I+4
        MA=MA+1
        L=L+1
        IA(L)=MA
        JA(MA)=I
        MA=MA+1
        JA(MA)=I+3
        MA=MA+1
        JA(MA)=I+4
        MA=MA+1
        L=L+1
        IA(L)=MA
        JA(MA)=I+1
        MA=MA+1
        JA(MA)=I+2
        MA=MA+1
        JA(MA)=I+3
        MA=MA+1
        L=L+1
        IA(L)=MA
  400 CONTINUE
      MA=IA(NA+1)-1
      RETURN
  410 CONTINUE
      IF (N.LT.5) GO TO 680
      IF (MOD(N-5,3).NE.0) N=N-MOD(N-5,3)
      DO 420 I=1,N
        X(I)=-1.0D0
  420 CONTINUE
      IA(1)=1
      KK=(N-5)/3+1
      NA=7*KK
      MA=1
      L=1
      DO 430 K=1,KK
        I=3*(K-1)+1
        JA(MA)=I
        MA=MA+1
        JA(MA)=I+1
        MA=MA+1
        L=L+1
        IA(L)=MA
        JA(MA)=I+1
        MA=MA+1
        JA(MA)=I+2
        MA=MA+1
        L=L+1
        IA(L)=MA
        JA(MA)=I+2
        MA=MA+1
        JA(MA)=I+3
        MA=MA+1
        L=L+1
        IA(L)=MA
        JA(MA)=I+3
        MA=MA+1
        JA(MA)=I+4
        MA=MA+1
        L=L+1
        IA(L)=MA
        JA(MA)=I
        MA=MA+1
        JA(MA)=I+1
        MA=MA+1
        JA(MA)=I+2
        MA=MA+1
        L=L+1
        IA(L)=MA
        JA(MA)=I+1
        MA=MA+1
        JA(MA)=I+2
        MA=MA+1
        JA(MA)=I+3
        MA=MA+1
        L=L+1
        IA(L)=MA
        JA(MA)=I
        MA=MA+1
        JA(MA)=I+4
        MA=MA+1
        L=L+1
        IA(L)=MA
  430 CONTINUE
      MA=IA(NA+1)-1
      RETURN
  440 CONTINUE
      IF (N.LT.4) GO TO 680
      DO 450 I=1,N
        IF (MOD(I,4).EQ.1) THEN
          X(I)=-0.8D0
        ELSE IF (MOD(I,4).EQ.2) THEN
          X(I)=1.2D0
        ELSE IF (MOD(I,4).EQ.3) THEN
          X(I)=-1.2D0
        ELSE
          X(I)=0.8D0
        END IF
  450 CONTINUE
      Y(1)=14.4D0
      Y(2)=6.8D0
      Y(3)=4.2D0
      Y(4)=3.2D0
      II=4
      JJ=4
      LL=2
      IA(1)=1
  460 IF (MOD(N-II,LL).NE.0) N=N-MOD(N-II,LL)
      KK=(N-II)/LL+1
      NA=JJ*KK
      MA=1
      L=1
      DO 490 K=1,KK
        MM=(K-1)*LL
        DO 480 J=1,JJ
          DO 470 I=1,II
            JA(MA)=MM+I
            MA=MA+1
  470     CONTINUE
          L=L+1
          IA(L)=MA
  480   CONTINUE
  490 CONTINUE
      MA=IA(NA+1)-1
      RETURN
  500 CONTINUE
      IF (N.LT.4) GO TO 680
      DO 510 I=1,N
        IF (MOD(I,4).EQ.1) THEN
          X(I)=-0.8D0
        ELSE IF (MOD(I,4).EQ.2) THEN
          X(I)=1.2D0
        ELSE IF (MOD(I,4).EQ.3) THEN
          X(I)=-1.2D0
        ELSE
          X(I)=0.8D0
        END IF
  510 CONTINUE
      Y(1)=35.8D0
      Y(2)=11.2D0
      Y(3)=6.2D0
      Y(4)=4.4D0
      II=4
      JJ=4
      LL=2
      IA(1)=1
      GO TO 460
  520 CONTINUE
      IF (N.LT.4) GO TO 680
      DO 530 I=1,N
        IF (MOD(I,4).EQ.1) THEN
          X(I)=-0.8D0
        ELSE IF (MOD(I,4).EQ.2) THEN
          X(I)=1.2D0
        ELSE IF (MOD(I,4).EQ.3) THEN
          X(I)=-1.2D0
        ELSE
          X(I)=0.8D0
        END IF
  530 CONTINUE
      Y(1)=30.6D0
      Y(2)=72.2D0
      Y(3)=124.4D0
      Y(4)=187.4D0
      II=4
      JJ=4
      LL=2
      IA(1)=1
      GO TO 460
  540 IF (N.LT.4) GO TO 680
      N=N-MOD(N,2)
      NA=N
      MA=0
      IA(1)=1
      DO 550 I=1,N
        IF (MOD(I,8).EQ.1) X(I)=1.0D-1
        IF (MOD(I,8).EQ.2.OR.MOD(I,8).EQ.0) X(I)=2.0D-1
        IF (MOD(I,8).EQ.3.OR.MOD(I,8).EQ.7) X(I)=3.0D-1
        IF (MOD(I,8).EQ.4.OR.MOD(I,8).EQ.6) X(I)=4.0D-1
        IF (MOD(I,8).EQ.5) X(I)=5.0D-1
        IF (I.GT.2) THEN
          MA=MA+1
          JA(MA)=I-2
        END IF
        MA=MA+2
        IF (MOD(I,2).EQ.1) THEN
          JA(MA-1)=I
          JA(MA)=I+1
        ELSE
          JA(MA-1)=I-1
          JA(MA)=I
        END IF
        IF (I.LT.N-1) THEN
          MA=MA+1
          JA(MA)=I+2
        END IF
        IA(I+1)=MA+1
  550 CONTINUE
      RETURN
  560 IF (N.LT.3) GO TO 680
      NA=N
      MA=0
      IA(1)=1
      DO 570 I=1,N
        X(I)=1.2D1
        IF (I.GT.1) THEN
          MA=MA+1
          JA(MA)=I-1
        END IF
        MA=MA+1
        JA(MA)=I
        IF (I.LT.N) THEN
          MA=MA+1
          JA(MA)=I+1
        END IF
        IA(I+1)=MA+1
  570 CONTINUE
      RETURN
  580 IF (N.LT.7) GO TO 680
      NA=N
      MA=0
      IA(1)=1
      DO 590 I=1,N
        X(I)=-1.0D0
        IF (I.GT.1.AND.I.LT.N-3) THEN
          MA=MA+1
          JA(MA)=I-1
        END IF
        IF (I.LT.N-4) THEN
          MA=MA+1
          JA(MA)=I
        END IF
        IF (I.LT.N-5) THEN
          MA=MA+1
          JA(MA)=I+1
        END IF
        MA=MA+5
        JA(MA-4)=N-4
        JA(MA-3)=N-3
        JA(MA-2)=N-2
        JA(MA-1)=N-1
        JA(MA)=N
        IA(I+1)=MA+1
  590 CONTINUE
      RETURN
  600 IF (N.LT.3) GO TO 680
      NA=N
      MA=0
      IA(1)=1
      DO 610 I=1,N
        X(I)=DBLE(I)/DBLE(N+1)
        X(I)=X(I)*(X(I)-1.0D0)
        IF (I.GT.1) THEN
          MA=MA+1
          JA(MA)=I-1
        END IF
        MA=MA+1
        JA(MA)=I
        IF (I.LT.N) THEN
          MA=MA+1
          JA(MA)=I+1
        END IF
        IA(I+1)=MA+1
  610 CONTINUE
      RETURN
  620 CONTINUE
      IF (N.LT.5) GO TO 680
      IF (MOD(N-5,3).NE.0) N=N-MOD(N-5,3)
      DO 630 I=1,N
        X(I)=-1.0D0
  630 CONTINUE
      IA(1)=1
      KK=(N-5)/3+1
      NA=7*KK
      MA=1
      L=1
      DO 640 K=1,KK
        I=3*(K-1)+1
        JA(MA)=I
        MA=MA+1
        JA(MA)=I+1
        MA=MA+1
        L=L+1
        IA(L)=MA
        JA(MA)=I+1
        MA=MA+1
        JA(MA)=I+2
        MA=MA+1
        L=L+1
        IA(L)=MA
        JA(MA)=I+3
        MA=MA+1
        L=L+1
        IA(L)=MA
        JA(MA)=I+4
        MA=MA+1
        L=L+1
        IA(L)=MA
        JA(MA)=I
        MA=MA+1
        JA(MA)=I+1
        MA=MA+1
        L=L+1
        IA(L)=MA
        JA(MA)=I+2
        MA=MA+1
        JA(MA)=I+3
        MA=MA+1
        JA(MA)=I+4
        MA=MA+1
        L=L+1
        IA(L)=MA
        JA(MA)=I+1
        MA=MA+1
        JA(MA)=I+4
        MA=MA+1
        L=L+1
        IA(L)=MA
  640 CONTINUE
      MA=IA(NA+1)-1
      RETURN
  650 IF (N.LT.3) GO TO 680
      DO 660 I=1,N
        IF (MOD(I,2).EQ.1) THEN
          X(I)=-1.2D0
        ELSE
          X(I)=1.0D0
        END IF
  660 CONTINUE
      NA=2*(N-1)
      MA=2
      IA(1)=1
      IA(2)=2
      JA(1)=1
      DO 670 KA=1,NA-1
        I=(KA+1)/2
        IF (MOD(KA,2).EQ.1) THEN
          JA(MA)=I
          MA=MA+1
          JA(MA)=I+1
          MA=MA+1
        ELSE
          JA(MA)=I
          MA=MA+1
          JA(MA)=I+1
          MA=MA+1
          JA(MA)=I+2
          MA=MA+1
        END IF
        IA(KA+2)=MA
  670 CONTINUE
      MA=MA-1
      RETURN
  680 IERR=1
      RETURN
      END
! SUBROUTINE TAFU15             ALL SYSTEMS                92/12/01
! PORTABILITY : ALL SYSTEMS
! 92/12/01 RA : ORIGINAL VERSION
!
! PURPOSE :
!  VALUES OF PARTIAL FUNCTIONS IN THE SUM OF SQUARES.
!
! PARAMETERS :
!  II  N  NUMBER OF VARIABLES.
!  II  KA  INDEX OF THE GIVEN PARTIAL FUNCTION.
!  RI  X(N)  VECTOR OF VARIABLES.
!  RO  FA  VALUE OF THE KA-TH PARTIAL FUNCTION AT THE POINT X.
!  II  NEXT  NUMBER OF THE SELECTED TEST PROBLEM.
!
      SUBROUTINE TAFU15 (N, KA, X, FA, NEXT)
      INTEGER N,KA,NEXT
      DOUBLE PRECISION X(*),FA
      DOUBLE PRECISION A,C,D,P,ALFA,U,V
      INTEGER I,J,K,L,M,IA,IB,IC
      DOUBLE PRECISION Y(20)
      COMMON /EMPR15/ Y
      GO TO (10,20,30,40,50,60,80,90,100,110,120,130,140,150,180,210,
     &230,240,250,260,270,280),NEXT
   10 I=(KA+1)/2
      IF (MOD(KA,2).EQ.1) THEN
        FA=1.0D1*(X(I)**2-X(I+1))
      ELSE
        FA=X(I)-1.0D0
      END IF
      RETURN
   20 A=SQRT(10.0D0)
      D=SQRT(90.0D0)
      I=2*((KA+5)/6)
      IF (MOD(KA,6).EQ.1) THEN
        FA=1.0D1*(X(I-1)**2-X(I))
      ELSE IF (MOD(KA,6).EQ.2) THEN
        FA=X(I-1)-1.0D0
      ELSE IF (MOD(KA,6).EQ.3) THEN
        FA=D*(X(I+1)**2-X(I+2))
      ELSE IF (MOD(KA,6).EQ.4) THEN
        FA=X(I+1)-1.0D0
      ELSE IF (MOD(KA,6).EQ.5) THEN
        FA=A*(X(I)+X(I+2)-2.0D0)
      ELSE
        FA=(X(I)-X(I+2))/A
      END IF
      RETURN
   30 A=SQRT(1.0D1)
      C=SQRT(5.0D0)
      I=2*((KA+3)/4)
      IF (MOD(KA,4).EQ.1) THEN
        FA=X(I-1)+1.0D1*X(I)
      ELSE IF (MOD(KA,4).EQ.2) THEN
        FA=C*(X(I+1)-X(I+2))
      ELSE IF (MOD(KA,4).EQ.3) THEN
        FA=(X(I)-2.0D0*X(I+1))**2
      ELSE
        FA=A*(X(I-1)-X(I+2))**2
      END IF
      RETURN
   40 I=2*((KA+4)/5)
      IF (MOD(KA,5).EQ.1) THEN
        FA=(EXP(X(I-1))-X(I))**2
      ELSE IF (MOD(KA,5).EQ.2) THEN
        FA=1.0D1*(X(I)-X(I+1))**3
      ELSE IF (MOD(KA,5).EQ.3) THEN
        P=X(I+1)-X(I+2)
        FA=(SIN(P)/COS(P))**2
      ELSE IF (MOD(KA,5).EQ.4) THEN
        FA=X(I-1)**4
      ELSE
        FA=X(I+2)-1.0D0
      END IF
      RETURN
   50 I=KA
      FA=(3.0D0-2.0D0*X(I))*X(I)+1.0D0
      IF (I.GT.1) FA=FA-X(I-1)
      IF (I.LT.N) FA=FA-X(I+1)
      RETURN
   60 I=KA
      FA=(2.0D0+5.0D0*X(I)**2)*X(I)+1.0D0
      DO 70 J=MAX(1,I-5),MIN(N,I+1)
        FA=FA+X(J)*(1.0D0+X(J))
   70 CONTINUE
      RETURN
   80 I=(KA+1)/2
      IF (MOD(KA,2).EQ.1) THEN
        FA=X(I)+X(I+1)*((5.0D0-X(I+1))*X(I+1)-2.0D0)-1.3D1
      ELSE
        FA=X(I)+X(I+1)*((1.0D0+X(I+1))*X(I+1)-1.4D1)-2.9D1
      END IF
      RETURN
   90 I=MOD(KA,N/2)+1
      J=I+N/2
      M=5*N
      IF (KA.LE.M/2) THEN
        IA=1
      ELSE
        IA=2
      END IF
      IB=5-KA/(M/4)
      IC=MOD(KA,5)+1
      FA=(X(I)**IA-X(J)**IB)**IC
      RETURN
  100 I=2*((KA+5)/6)-1
      IF (MOD(KA,6).EQ.1) THEN
        FA=X(I)+3.0D0*X(I+1)*(X(I+2)-1.0D0)+X(I+3)**2-1.0D0
      ELSE IF (MOD(KA,6).EQ.2) THEN
        FA=(X(I)+X(I+1))**2+(X(I+2)-1.0D0)**2-X(I+3)-3.0D0
      ELSE IF (MOD(KA,6).EQ.3) THEN
        FA=X(I)*X(I+1)-X(I+2)*X(I+3)
      ELSE IF (MOD(KA,6).EQ.4) THEN
        FA=2.0D0*X(I)*X(I+2)+X(I+1)*X(I+3)-3.0D0
      ELSE IF (MOD(KA,6).EQ.5) THEN
        FA=(X(I)+X(I+1)+X(I+2)+X(I+3))**2+(X(I)-1.0D0)**2
      ELSE
        FA=X(I)*X(I+1)*X(I+2)*X(I+3)+(X(I+3)-1.0D0)**2-1.0D0
      END IF
      RETURN
  110 I=(KA+1)/2
      J=MOD(KA,2)
      IF (J.EQ.0) THEN
        FA=6.0D0-EXP(2.0D0*X(I))-EXP(2.0D0*X(I+1))
      ELSE IF (I.EQ.1) THEN
        FA=4.0D0-EXP(X(I))-EXP(X(I+1))
      ELSE IF (I.EQ.N) THEN
        FA=8.0D0-EXP(3.0D0*X(I-1))-EXP(3.0D0*X(I))
      ELSE
        FA=8.0D0-EXP(3.0D0*X(I-1))-EXP(3.0D0*X(I))+4.0D0-EXP(X(I))-
     &   EXP(X(I+1))
      END IF
      RETURN
  120 I=(KA+1)/2
      IF (MOD(KA,2).EQ.1) THEN
        FA=1.0D1*(2.0D0*X(I)/(1.0D0+X(I)**2)-X(I+1))
      ELSE
        FA=X(I)-1.0D0
      END IF
      RETURN
  130 I=3*((KA+5)/6)-2
      IF (MOD(KA,6).EQ.1) THEN
        FA=1.0D1*(X(I)**2-X(I+1))
      ELSE IF (MOD(KA,6).EQ.2) THEN
        FA=X(I+2)-1.0D0
      ELSE IF (MOD(KA,6).EQ.3) THEN
        FA=(X(I+3)-1.0D0)**2
      ELSE IF (MOD(KA,6).EQ.4) THEN
        FA=(X(I+4)-1.0D0)**3
      ELSE IF (MOD(KA,6).EQ.5) THEN
        FA=X(I)**2*X(I+3)+SIN(X(I+3)-X(I+4))-1.0D1
      ELSE
        FA=X(I+1)+(X(I+2)**2*X(I+3))**2-2.0D1
      END IF
      RETURN
  140 I=3*((KA+6)/7)-2
      IF (MOD(KA,7).EQ.1) THEN
        FA=1.0D1*(X(I)**2-X(I+1))
      ELSE IF (MOD(KA,7).EQ.2) THEN
        FA=1.0D1*(X(I+1)**2-X(I+2))
      ELSE IF (MOD(KA,7).EQ.3) THEN
        FA=(X(I+2)-X(I+3))**2
      ELSE IF (MOD(KA,7).EQ.4) THEN
        FA=(X(I+3)-X(I+4))**2
      ELSE IF (MOD(KA,7).EQ.5) THEN
        FA=X(I)+X(I+1)**2+X(I+2)-3.0D1
      ELSE IF (MOD(KA,7).EQ.6) THEN
        FA=X(I+1)-X(I+2)**2+X(I+3)-1.0D1
      ELSE
        FA=X(I)*X(I+4)-1.0D1
      END IF
      RETURN
  150 I=2*((KA+3)/4)-2
      L=MOD((KA-1),4)+1
      FA=-Y(L)
      DO 170 K=1,3
        A=DBLE(K*K)/DBLE(L)
        DO 160 J=1,4
          IF (X(I+J).EQ.0) X(I+J)=1.0D-16
          A=A*SIGN(1.0D0,X(I+J))*ABS(X(I+J))**(DBLE(J)/DBLE(K*L))
  160   CONTINUE
        FA=FA+A
  170 CONTINUE
      RETURN
  180 I=2*((KA+3)/4)-2
      L=MOD((KA-1),4)+1
      FA=-Y(L)
      DO 200 K=1,3
        A=0.0D0
        DO 190 J=1,4
          A=A+X(I+J)*(DBLE(J)/DBLE(K*L))
  190   CONTINUE
        FA=FA+EXP(A)*DBLE(K*K)/DBLE(L)
  200 CONTINUE
      RETURN
  210 I=2*((KA+3)/4)-2
      L=MOD((KA-1),4)+1
      FA=-Y(L)
      DO 220 J=1,4
        FA=FA+DBLE((1-2*MOD(J,2))*L*J*J)*SIN(X(I+J))+DBLE(L*L*J)*
     &   COS(X(I+J))
  220 CONTINUE
      RETURN
  230 ALFA=0.5D0
      IF (KA.EQ.1) THEN
        FA=ALFA-(1.0D0-ALFA)*X(3)-X(1)*(1.0D0+4.0D0*X(2))
      ELSE IF (KA.EQ.2) THEN
        FA=-(2.0D0-ALFA)*X(4)-X(2)*(1.0D0+4.0D0*X(1))
      ELSE IF (KA.EQ.N-1) THEN
        FA=ALFA*X(N-3)-X(N-1)*(1.0D0+4.0D0*X(N))
      ELSE IF (KA.EQ.N) THEN
        FA=ALFA*X(N-2)-(2.0D0-ALFA)-X(N)*(1.0D0+4.0D0*X(N-1))
      ELSE IF (MOD(KA,2).EQ.1) THEN
        FA=ALFA*X(KA-2)-(1.0D0-ALFA)*X(KA+2)-X(KA)*(1.0D0+4.0D0*X(KA+1))
      ELSE
        FA=ALFA*X(KA-2)-(2.0D0-ALFA)*X(KA+2)-X(KA)*(1.0D0+4.0D0*X(KA-1))
      END IF
      RETURN
  240 IF (KA.LT.2) THEN
        FA=4.0D0*(X(KA)-X(KA+1)**2)
      ELSE IF (KA.LT.N) THEN
        FA=8.0D0*X(KA)*(X(KA)**2-X(KA-1))-2.0D0*(1.0D0-X(KA))+4.0D0*
     &   (X(KA)-X(KA+1)**2)
      ELSE
        FA=8.0D0*X(KA)*(X(KA)**2-X(KA-1))-2.0D0*(1.0D0-X(KA))
      END IF
      RETURN
  250 IF (KA.EQ.1) THEN
        FA=-2.0D0*X(KA)**2+3.0D0*X(KA)-2.0D0*X(KA+1)+3.0D0*X(N-4)-X(N-3)
     &   -X(N-2)+0.5D0*X(N-1)-X(N)+1.0D0
      ELSE IF (KA.LE.N-1) THEN
        FA=-2.0D0*X(KA)**2+3.0D0*X(KA)-X(KA-1)-2.0D0*X(KA+1)+3.0D0*X(N-
     &   4)-X(N-3)-X(N-2)+0.5D0*X(N-1)-X(N)+1.0D0
      ELSE
        FA=-2.0D0*X(N)**2+3.0D0*X(N)-X(N-1)+3.0D0*X(N-4)-X(N-3)-X(N-2)+
     &   0.5D0*X(N-1)-X(N)+1.0D0
      END IF
      RETURN
  260 U=1.0D0/DBLE(N+1)
      V=DBLE(KA)*U
      FA=2.0D0*X(KA)+0.5D0*U*U*(X(KA)+V+1.0D0)**3+1.0D0
      IF (KA.GT.1) FA=FA-X(KA-1)
      IF (KA.LT.N) FA=FA-X(KA+1)
      RETURN
  270 I=3*((KA+6)/7)-2
      IF (MOD(KA,7).EQ.1) THEN
        FA=1.0D1*(X(I)**2-X(I+1))
      ELSE IF (MOD(KA,7).EQ.2) THEN
        FA=X(I+1)+X(I+2)-2.0D0
      ELSE IF (MOD(KA,7).EQ.3) THEN
        FA=X(I+3)-1.0D0
      ELSE IF (MOD(KA,7).EQ.4) THEN
        FA=X(I+4)-1.0D0
      ELSE IF (MOD(KA,7).EQ.5) THEN
        FA=X(I)+3.0D0*X(I+1)
      ELSE IF (MOD(KA,7).EQ.6) THEN
        FA=X(I+2)+X(I+3)-2.0D0*X(I+4)
      ELSE
        FA=1.0D1*(X(I+1)**2-X(I+4))
      END IF
      RETURN
  280 I=KA/2
      IF (KA.EQ.1) THEN
        FA=X(KA)-1.0D0
      ELSE IF (MOD(KA,2).EQ.0) THEN
        FA=1.0D1*(X(I)**2-X(I+1))
      ELSE
        FA=2.0D0*EXP(-(X(I)-X(I+1))**2)+EXP(-2.0D0*(X(I+1)-X(I+2))**2)
      END IF
      RETURN
      END
! SUBROUTINE TAGU15                ALL SYSTEMS                92/12/01
! PORTABILITY : ALL SYSTEMS
! 92/12/01 RA : ORIGINAL VERSION
!
! PURPOSE :
!  GRADIENTS OF PARTIAL FUNCTIONS IN THE SUM OF SQUARES.
!
! PARAMETERS :
!  II  N  NUMBER OF VARIABLES.
!  II  KA  INDEX OF THE GIVEN PARTIAL FUNCTION.
!  RI  X(N)  VECTOR OF VARIABLES.
!  RO  GA(N)  GRADIENT OF THE KA-TH PARTIAL FUNCTION AT THE POINT X.
!  II  NEXT  NUMBER OF THE SELECTED TEST PROBLEM.
!
      SUBROUTINE TAGU15 (N, KA, X, GA, NEXT)
      INTEGER N,KA,NEXT
      DOUBLE PRECISION X(*),GA(*)
      DOUBLE PRECISION A,B,C,D,Q,R,E,ALFA,U,V
      INTEGER I,J,K,L,M,IA,IB,IC
      DOUBLE PRECISION Y(20)
      COMMON /EMPR15/ Y
      GO TO (10,20,30,40,50,60,80,90,100,110,120,130,140,150,200,250,
     &270,280,290,300,310,320,330,340),NEXT
   10 I=(KA+1)/2
      IF (MOD(KA,2).EQ.1) THEN
        GA(I)=2.0D1*X(I)
        GA(I+1)=-1.0D1
      ELSE
        GA(I)=1.0D0
      END IF
      RETURN
   20 I=2*((KA+5)/6)
      A=SQRT(9.0D1)
      B=SQRT(1.0D1)
      IF (MOD(KA,6).EQ.1) THEN
        GA(I-1)=2.0D1*X(I-1)
        GA(I)=-1.0D1
      ELSE IF (MOD(KA,6).EQ.2) THEN
        GA(I-1)=1.0D0
      ELSE IF (MOD(KA,6).EQ.3) THEN
        GA(I+1)=2.0D0*A*X(I+1)
        GA(I+2)=-A
      ELSE IF (MOD(KA,6).EQ.4) THEN
        GA(I+1)=1.0D0
      ELSE IF (MOD(KA,6).EQ.5) THEN
        GA(I)=B
        GA(I+2)=B
      ELSE
        GA(I)=1.0D0/B
        GA(I+2)=-1.0D0/B
      END IF
      RETURN
   30 I=2*((KA+3)/4)
      A=SQRT(5.0D0)
      B=SQRT(1.0D1)
      IF (MOD(KA,4).EQ.1) THEN
        GA(I-1)=1.0D0
        GA(I)=1.0D1
      ELSE IF (MOD(KA,4).EQ.2) THEN
        GA(I+1)=A
        GA(I+2)=-A
      ELSE IF (MOD(KA,4).EQ.3) THEN
        C=X(I)-2.0D0*X(I+1)
        GA(I)=2.0D0*C
        GA(I+1)=-4.0D0*C
      ELSE
        C=X(I-1)-X(I+2)
        D=2.0D0*B*C
        GA(I-1)=D
        GA(I+2)=-D
      END IF
      RETURN
   40 I=2*((KA+4)/5)
      IF (MOD(KA,5).EQ.1) THEN
        A=EXP(X(I-1))
        B=A-X(I)
        C=2.0D0*B
        GA(I-1)=C*A
        GA(I)=-C
      ELSE IF (MOD(KA,5).EQ.2) THEN
        A=X(I)-X(I+1)
        B=3.0D1*A**2
        GA(I)=B
        GA(I+1)=-B
      ELSE IF (MOD(KA,5).EQ.3) THEN
        C=X(I+1)-X(I+2)
        Q=SIN(C)/COS(C)
        R=COS(C)
        D=2.0D0*Q/R**2
        GA(I+1)=D
        GA(I+2)=-D
      ELSE IF (MOD(KA,5).EQ.4) THEN
        GA(I-1)=4.0D0*X(I-1)**3
      ELSE
        GA(I+2)=1.0D0
      END IF
      RETURN
   50 I=KA
      GA(I)=3.0D0-4.0D0*X(I)
      IF (I.GT.1) GA(I-1)=-1.0D0
      IF (I.LT.N) GA(I+1)=-1.0D0
      RETURN
   60 I=KA
      DO 70 J=MAX(1,I-5),MIN(N,I+1)
        GA(J)=1.0D0+2.0D0*X(J)
   70 CONTINUE
      GA(I)=GA(I)+2.0D0+1.5D1*X(I)**2
      RETURN
   80 I=(KA+1)/2
      IF (MOD(KA,2).EQ.1) THEN
        GA(I)=1.0D0
        GA(I+1)=1.0D1*X(I+1)-3.0D0*X(I+1)**2-2.0D0
      ELSE
        GA(I)=1.0D0
        GA(I+1)=2.0D0*X(I+1)+3.0D0*X(I+1)**2-1.4D1
      END IF
      RETURN
   90 I=MOD(KA,N/2)+1
      J=I+N/2
      M=5*N
      IF (KA.LE.M/2) THEN
        IA=1
      ELSE
        IA=2
      END IF
      IB=5-KA/(M/4)
      IC=MOD(KA,5)+1
      A=DBLE(IA)
      B=DBLE(IB)
      C=DBLE(IC)
      D=X(I)**IA-X(J)**IB
      IF (D.EQ.0.0D0) THEN
        GA(I)=0.0D0
        GA(J)=0.0D0
      ELSE
        E=C*D**(IC-1)
        IF (X(I).EQ.0.0D0.AND.IA.LE.1) THEN
          GA(I)=0.0D0
        ELSE
          GA(I)=E*A*X(I)**(IA-1)
        END IF
        IF (X(J).EQ.0.0D0.AND.IB.LE.1) THEN
          GA(J)=0.0D0
        ELSE
          GA(J)=-E*B*X(J)**(IB-1)
        END IF
      END IF
      RETURN
  100 I=2*((KA+5)/6)-1
      IF (MOD(KA,6).EQ.1) THEN
        GA(I)=1.0D0
        GA(I+1)=3.0D0*(X(I+2)-1.0D0)
        GA(I+2)=3.0D0*X(I+1)
        GA(I+3)=2.0D0*X(I+3)
      ELSE IF (MOD(KA,6).EQ.2) THEN
        GA(I)=2.0D0*(X(I)+X(I+1))
        GA(I+1)=2.0D0*(X(I)+X(I+1))
        GA(I+2)=2.0D0*(X(I+2)-1.0D0)
        GA(I+3)=-1.0D0
      ELSE IF (MOD(KA,6).EQ.3) THEN
        GA(I)=X(I+1)
        GA(I+1)=X(I)
        GA(I+2)=-X(I+3)
        GA(I+3)=-X(I+2)
      ELSE IF (MOD(KA,6).EQ.4) THEN
        GA(I)=2.0D0*X(I+2)
        GA(I+1)=X(I+3)
        GA(I+2)=2.0D0*X(I)
        GA(I+3)=X(I+1)
      ELSE IF (MOD(KA,6).EQ.5) THEN
        GA(I)=2.0D0*(X(I)+X(I+1)+X(I+2)+X(I+3))+2.0D0*(X(I)-1.0D0)
        GA(I+1)=2.0D0*(X(I)+X(I+1)+X(I+2)+X(I+3))
        GA(I+2)=2.0D0*(X(I)+X(I+1)+X(I+2)+X(I+3))
        GA(I+3)=2.0D0*(X(I)+X(I+1)+X(I+2)+X(I+3))
      ELSE
        GA(I)=X(I+1)*X(I+2)*X(I+3)
        GA(I+1)=X(I)*X(I+2)*X(I+3)
        GA(I+2)=X(I)*X(I+1)*X(I+3)
        GA(I+3)=X(I)*X(I+1)*X(I+2)+2.0D0*(X(I+3)-1.0D0)
      END IF
      RETURN
  110 IF (N.GE.2) THEN
        I=(KA+1)/2
        J=MOD(KA,2)
        IF (J.EQ.0) THEN
          GA(I)=-2.0D0*EXP(2.0D0*X(I))
          GA(I+1)=-2.0D0*EXP(2.0D0*X(I+1))
        ELSE IF (I.EQ.1) THEN
          GA(I)=-EXP(X(I))
          GA(I+1)=-EXP(X(I+1))
        ELSE IF (I.EQ.N) THEN
          GA(I-1)=-3.0D0*EXP(3.0D0*X(I-1))
          GA(I)=-3.0D0*EXP(3.0D0*X(I))
        ELSE
          GA(I-1)=-3.0D0*EXP(3.0D0*X(I-1))
          GA(I)=-3.0D0*EXP(3.0D0*X(I))-EXP(X(I))
          GA(I+1)=-EXP(X(I+1))
        END IF
      END IF
      RETURN
  120 I=(KA+1)/2
      IF (MOD(KA,2).EQ.1) THEN
        GA(I)=2.0D1*(1.0D0-X(I)**2)/(1.0D0+X(I)**2)**2
        GA(I+1)=-1.0D1
      ELSE
        GA(I)=1.0D0
      END IF
      RETURN
  130 I=3*((KA+5)/6)-2
      IF (MOD(KA,6).EQ.1) THEN
        GA(I)=2.0D1*X(I)
        GA(I+1)=-1.0D1
      ELSE IF (MOD(KA,6).EQ.2) THEN
        GA(I+2)=1.0D0
      ELSE IF (MOD(KA,6).EQ.3) THEN
        GA(I+3)=2.0D0*(X(I+3)-1.0D0)
      ELSE IF (MOD(KA,6).EQ.4) THEN
        GA(I+4)=3.0D0*(X(I+4)-1.0D0)**2
      ELSE IF (MOD(KA,6).EQ.5) THEN
        GA(I)=2.0D0*X(I)*X(I+3)
        GA(I+3)=X(I)**2+COS(X(I+3)-X(I+4))
        GA(I+4)=-COS(X(I+3)-X(I+4))
      ELSE
        GA(I+1)=1.0D0
        GA(I+2)=4.0D0*X(I+2)*(X(I+2)*X(I+3))**2
        GA(I+3)=2.0D0*X(I+2)**4*X(I+3)
      END IF
      RETURN
  140 I=3*((KA+6)/7)-2
      IF (MOD(KA,7).EQ.1) THEN
        GA(I)=2.0D1*X(I)
        GA(I+1)=-1.0D1
      ELSE IF (MOD(KA,7).EQ.2) THEN
        GA(I+1)=2.0D1*X(I+1)
        GA(I+2)=-1.0D1
      ELSE IF (MOD(KA,7).EQ.3) THEN
        GA(I+2)=2.0D0*(X(I+2)-X(I+3))
        GA(I+3)=-2.0D0*(X(I+2)-X(I+3))
      ELSE IF (MOD(KA,7).EQ.4) THEN
        GA(I+3)=2.0D0*(X(I+3)-X(I+4))
        GA(I+4)=-2.0D0*(X(I+3)-X(I+4))
      ELSE IF (MOD(KA,7).EQ.5) THEN
        GA(I)=1.0D0
        GA(I+1)=2.0D0*X(I+1)
        GA(I+2)=1.0D0
      ELSE IF (MOD(KA,7).EQ.6) THEN
        GA(I+1)=1.0D0
        GA(I+2)=-2.0D0*X(I+2)
        GA(I+3)=1.0D0
      ELSE
        GA(I)=X(I+4)
        GA(I+4)=X(I)
      END IF
      RETURN
  150 I=2*((KA+3)/4)-2
      L=MOD((KA-1),4)+1
      DO 160 J=1,4
        GA(I+J)=0.0D0
  160 CONTINUE
      DO 190 K=1,3
        A=DBLE(K*K)/DBLE(L)
        DO 170 J=1,4
          A=A*SIGN(1.0D0,X(I+J))*ABS(X(I+J))**(DBLE(J)/DBLE(K*L))
  170   CONTINUE
        DO 180 J=1,4
          GA(I+J)=GA(I+J)+(DBLE(J)/DBLE(K*L))*A/X(I+J)
  180   CONTINUE
  190 CONTINUE
      RETURN
  200 I=2*((KA+3)/4)-2
      L=MOD((KA-1),4)+1
      DO 210 J=1,4
        GA(I+J)=0.0D0
  210 CONTINUE
      DO 240 K=1,3
        A=0.0D0
        DO 220 J=1,4
          A=A+X(I+J)*(DBLE(J)/DBLE(K*L))
  220   CONTINUE
        A=EXP(A)*DBLE(K*K)/DBLE(L)
        DO 230 J=1,4
          GA(I+J)=GA(I+J)+A*(DBLE(J)/DBLE(K*L))
  230   CONTINUE
  240 CONTINUE
      RETURN
  250 I=2*((KA+3)/4)-2
      L=MOD((KA-1),4)+1
      DO 260 J=1,4
        GA(I+J)=DBLE((1-2*MOD(J,2))*L*J*J)*COS(X(I+J))-DBLE(L*L*J)*
     &   SIN(X(I+J))
  260 CONTINUE
      RETURN
  270 ALFA=0.5D0
      IF (KA.EQ.1) THEN
        GA(1)=-1.0D0-4.0D0*X(2)
        GA(2)=-4.0D0*X(1)
        GA(3)=ALFA-1.0D0
      ELSE IF (KA.EQ.2) THEN
        GA(1)=-4.0D0*X(2)
        GA(2)=-1.0D0-4.0D0*X(1)
        GA(4)=-2.0D0+ALFA
      ELSE IF (KA.EQ.N-1) THEN
        GA(N-3)=ALFA
        GA(N-1)=-1.0D0-4.0D0*X(N)
        GA(N)=-4.0D0*X(N-1)
      ELSE IF (KA.EQ.N) THEN
        GA(N-2)=ALFA
        GA(N-1)=-4.0D0*X(N)
        GA(N)=-1.0D0-4.0D0*X(N-1)
      ELSE IF (MOD(KA,2).EQ.1) THEN
        GA(KA-2)=ALFA
        GA(KA)=-1.0D0-4.0D0*X(KA+1)
        GA(KA+1)=-4.0D0*X(KA)
        GA(KA+2)=-1.0D0+ALFA
      ELSE
        GA(KA-2)=ALFA
        GA(KA-1)=-4.0D0*X(KA)
        GA(KA)=-1.0D0-4.0D0*X(KA-1)
        GA(KA+2)=-2.0D0+ALFA
      END IF
      RETURN
  280 IF (KA.LT.2) THEN
        GA(KA)=4.0D0
        GA(KA+1)=-8.0D0*X(KA+1)
      ELSE IF (KA.LT.N) THEN
        GA(KA-1)=-8.0D0*X(KA)
        GA(KA)=24.0D0*X(KA)**2-8.0D0*X(KA-1)+6.0D0
        GA(KA+1)=-8.0D0*X(KA+1)
      ELSE
        GA(KA-1)=-8.0D0*X(KA)
        GA(KA)=24.0D0*X(KA)**2-8.0D0*X(KA-1)+2.0D0
      END IF
      RETURN
  290 IF (KA.EQ.1) THEN
        GA(N-4)=3.0D0
        GA(N-3)=-1.0D0
        GA(N-2)=-1.0D0
        GA(N-1)=0.50D0
        GA(N)=-1.0D0
        GA(KA)=-4.0D0*X(KA)+3.0D0
        GA(KA+1)=-2.0D0
      ELSE IF (KA.LE.N-1) THEN
        GA(KA-1)=0.0D0
        GA(KA)=0.0D0
        GA(KA+1)=0.0D0
        GA(N-4)=3.0D0
        GA(N-3)=-1.0D0
        GA(N-2)=-1.0D0
        GA(N-1)=0.50D0
        GA(N)=-1.0D0
        GA(KA-1)=GA(KA-1)-1.0D0
        GA(KA)=GA(KA)-4.0D0*X(KA)+3.0D0
        GA(KA+1)=GA(KA+1)-2.0D0
      ELSE
        GA(N-4)=3.0D0
        GA(N-3)=-1.0D0
        GA(N-2)=-1.0D0
        GA(N-1)=0.50D0
        GA(N)=-4.0D0*X(N)+2.0D0
      END IF
      RETURN
  300 U=1.0D0/DBLE(N+1)
      V=DBLE(KA)*U
      GA(KA)=2.0D0+1.5D0*U**2*(X(KA)+V+1.0D0)**2
      IF (KA.GT.1) GA(KA-1)=-1.0D0
      IF (KA.LT.N) GA(KA+1)=-1.0D0
      RETURN
  310 I=3*((KA+6)/7)-2
      IF (MOD(KA,7).EQ.1) THEN
        GA(I)=2.0D1*X(I)
        GA(I+1)=-1.0D1
      ELSE IF (MOD(KA,7).EQ.2) THEN
        GA(I+1)=1.0D0
        GA(I+2)=1.0D0
      ELSE IF (MOD(KA,7).EQ.3) THEN
        GA(I+3)=1.0D0
      ELSE IF (MOD(KA,7).EQ.4) THEN
        GA(I+4)=1.0D0
      ELSE IF (MOD(KA,7).EQ.5) THEN
        GA(I)=1.0D0
        GA(I+1)=3.0D0
      ELSE IF (MOD(KA,7).EQ.6) THEN
        GA(I+2)=1.0D0
        GA(I+3)=1.0D0
        GA(I+4)=-2.0D0
      ELSE
        GA(I+1)=2.0D1*X(I+1)
        GA(I+4)=-1.0D1
      END IF
      RETURN
  320 I=KA/2
      IF (KA.EQ.1) THEN
        GA(KA)=1.0D0
      ELSE IF (MOD(KA,2).EQ.0) THEN
        GA(I)=2.0D1*X(I)
        GA(I+1)=-1.0D1
      ELSE
        A=2.0D0*EXP(-(X(I)-X(I+1))**2)
        B=EXP(-2.0D0*(X(I+1)-X(I+2))**2)
        GA(I)=-2.0D0*(X(I)-X(I+1))*A
        GA(I+1)=2.0D0*(X(I)-X(I+1))*A-4.0D0*(X(I+1)-X(I+2))*B
        GA(I+2)=4.0D0*(X(I+1)-X(I+2))*B
      END IF
      RETURN
  330 I=MOD(KA,N/2)+1
      J=I+N/2
      M=5*N
      IA=KA/(M/4)+1
      A=DBLE(IA)
      B=DBLE(KA/(M/5)+1)
      GA(I)=A*X(I)**(IA-1)*EXP(B*X(J))
      GA(J)=B*X(I)**IA*EXP(B*X(J))+1.0D0
      RETURN
  340 IA=MIN(MAX(MOD(KA,13)-2,1),7)
      IB=(KA+12)/13
      I=IA+IB-1
      IF (IA.EQ.7) THEN
        J=IB
      ELSE
        J=IA+IB
      END IF
      C=3.0D0*DBLE(IA)/1.0D1
      A=COS(C)
      B=EXP(SIN(C*X(J)))
      GA(I)=-COS(X(I))*(1.0D0+A)+5.0D0*B
      GA(J)=(1.0D0+A)+5.0D0*(X(I)-2.0D0)*C*COS(C*X(J))*B
      DO 350 L=0,6
        IF (IB+L.NE.I.AND.IB+L.NE.J) GA(IB+L)=0.5D0*COS(X(IB+L))
  350 CONTINUE
      RETURN
      END
! SUBROUTINE TIUB18             ALL SYSTEMS                 97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
!  INITIATION OF VARIABLES AND DEFINITION OF STRUCTURE OF THE SPARSE
!  JACOBIAN MATRIX FOR THE SYSTEM OF NONLINEAR EQUATIONS.
!
! PARAMETERS :
!  IU  N  NUMBER OF VARIABLES.
!  IO  NA  NUMBER OF EQUATIONS.
!  IO  MA  NUMBER OF NONZERO ELEMENTS IN THE SPARSE JACOBIAN MATRIX.
!  RO  X(N)  VECTOR OF VARIABLES.
!  IO  IA(NA+1)  POINTERS OF FIRST IN THE ROW ELEMENTS IN THE
!         SPARSE JACOBIAN MATRIX.
!  IO  JA(MA)  COLUMN INDICES OF NONZERO ELEMENTS IN THE SPARSE
!         JACOBIAN MATRIX.
!  RO  FMIN  LOWER BOUND FOR THE OBJECTIVE FUNCTION VALUE.
!  RO  XMAX  MAXIMUM ALLOWED STEPSIZE.
!  II  NEXT  NUMBER OF THE TEST PROBLEM.
!  IO  IERR  ERROR INDICATOR.
!
      SUBROUTINE TIUB18 (N, NA, MA, X, IA, JA, FMIN, XMAX, NEXT, IERR)
      INTEGER N,NA,MA,IA(*),JA(*),NEXT,IERR
      DOUBLE PRECISION X(*),FMIN,XMAX,PAR
      INTEGER I,J,K,M
      COMMON /EMPR18/ PAR,M
      FMIN=0.0D0
      XMAX=1.0D4
      IERR=0
      GO TO (10,30,50,70,110,130,150,170,190,210,230,250,270,290,310,
     &370,410,580,430,450,560,780,800,470,660,720,630,690,600,500),NEXT
   10 IF (N.LT.4) GO TO 830
      N=N-MOD(N,2)
      NA=N
      MA=0
      IA(1)=1
      DO 20 I=1,N
        IF (MOD(I,8).EQ.1) X(I)=1.0D-1
        IF (MOD(I,8).EQ.2.OR.MOD(I,8).EQ.0) X(I)=2.0D-1
        IF (MOD(I,8).EQ.3.OR.MOD(I,8).EQ.7) X(I)=3.0D-1
        IF (MOD(I,8).EQ.4.OR.MOD(I,8).EQ.6) X(I)=4.0D-1
        IF (MOD(I,8).EQ.5) X(I)=5.0D-1
        IF (I.GT.2) THEN
          MA=MA+1
          JA(MA)=I-2
        END IF
        MA=MA+2
        IF (MOD(I,2).EQ.1) THEN
          JA(MA-1)=I
          JA(MA)=I+1
        ELSE
          JA(MA-1)=I-1
          JA(MA)=I
        END IF
        IF (I.LT.N-1) THEN
          MA=MA+1
          JA(MA)=I+2
        END IF
        IA(I+1)=MA+1
   20 CONTINUE
      RETURN
   30 IF (N.LT.5) GO TO 830
      N=N-MOD(N,2)
      NA=N
      MA=0
      IA(1)=1
      DO 40 I=1,N
        IF (MOD(I,8).EQ.1) X(I)=1.0D-1
        IF (MOD(I,8).EQ.2.OR.MOD(I,8).EQ.0) X(I)=2.0D-1
        IF (MOD(I,8).EQ.3.OR.MOD(I,8).EQ.7) X(I)=3.0D-1
        IF (MOD(I,8).EQ.4.OR.MOD(I,8).EQ.6) X(I)=4.0D-1
        IF (MOD(I,8).EQ.5) X(I)=5.0D-1
        MA=MA+1
        JA(MA)=1
        IF (I.GT.3) THEN
          MA=MA+1
          JA(MA)=I-2
        END IF
        IF (I.GT.2) THEN
          MA=MA+1
          JA(MA)=I-1
        END IF
        IF (I.GT.1) THEN
          MA=MA+1
          JA(MA)=I
        END IF
        IF (I.EQ.1) THEN
          MA=MA+1
          JA(MA)=I+1
        END IF
        IF (I.LT.N-1) THEN
          MA=MA+1
          JA(MA)=I+2
        END IF
        IA(I+1)=MA+1
   40 CONTINUE
      RETURN
   50 IF (N.LT.5) GO TO 830
      N=N-MOD(N,5)
      NA=N
      MA=0
      IA(1)=1
      DO 60 I=1,N
        X(I)=1.0D0/DBLE(N)
        J=(I-1)/5*5
        MA=MA+5
        JA(MA-4)=J+1
        JA(MA-3)=J+2
        JA(MA-2)=J+3
        JA(MA-1)=J+4
        JA(MA)=J+5
        IA(I+1)=MA+1
   60 CONTINUE
      RETURN
   70 IF (N.LT.3) GO TO 830
      DO 80 I=1,N
        X(I)=0.0D0
   80 CONTINUE
   90 NA=N
      MA=0
      IA(1)=1
      DO 100 I=1,N
        IF (I.GT.1) THEN
          MA=MA+1
          JA(MA)=I-1
        END IF
        MA=MA+1
        JA(MA)=I
        IF (I.LT.N) THEN
          MA=MA+1
          JA(MA)=I+1
        END IF
        IA(I+1)=MA+1
  100 CONTINUE
      RETURN
  110 IF (N.LT.3) GO TO 830
      IF (MOD(N,2).NE.1) N=N-1
      NA=N
      MA=0
      IA(1)=1
      DO 120 I=1,N
        X(I)=1.0D0
        IF (MOD(I,2).EQ.1) THEN
          IF (I.GT.1) THEN
            MA=MA+2
            JA(MA-1)=I-2
            JA(MA)=I-1
          END IF
          MA=MA+1
          JA(MA)=I
          IF (I.LT.N) THEN
            MA=MA+2
            JA(MA-1)=I+1
            JA(MA)=I+2
          END IF
        ELSE
          MA=MA+3
          JA(MA-2)=I-1
          JA(MA-1)=I
          JA(MA)=I+1
        END IF
        IA(I+1)=MA+1
  120 CONTINUE
      RETURN
  130 IF (N.LT.3) GO TO 830
      DO 140 I=1,N
        X(I)=-1.0D0
  140 CONTINUE
      GO TO 90
  150 IF (N.LT.3) GO TO 830
      DO 160 I=1,N
        X(I)=1.2D1
  160 CONTINUE
      XMAX=1.0D2
      GO TO 90
  170 IF (N.LT.5) GO TO 830
      NA=N
      MA=0
      IA(1)=1
      DO 180 I=1,N
        X(I)=-2.0D0
        IF (I.GT.2) THEN
          MA=MA+1
          JA(MA)=I-2
        END IF
        IF (I.GT.1) THEN
          MA=MA+1
          JA(MA)=I-1
        END IF
        MA=MA+1
        JA(MA)=I
        IF (I.LT.N) THEN
          MA=MA+1
          JA(MA)=I+1
        END IF
        IF (I.LT.N-1) THEN
          MA=MA+1
          JA(MA)=I+2
        END IF
        IA(I+1)=MA+1
  180 CONTINUE
      RETURN
  190 IF (N.LT.7) GO TO 830
      NA=N
      MA=0
      IA(1)=1
      DO 200 I=1,N
        X(I)=-3.0D0
        IF (I.GT.3) THEN
          MA=MA+1
          JA(MA)=I-3
        END IF
        IF (I.GT.2) THEN
          MA=MA+1
          JA(MA)=I-2
        END IF
        IF (I.GT.1) THEN
          MA=MA+1
          JA(MA)=I-1
        END IF
        MA=MA+1
        JA(MA)=I
        IF (I.LT.N) THEN
          MA=MA+1
          JA(MA)=I+1
        END IF
        IF (I.LT.N-1) THEN
          MA=MA+1
          JA(MA)=I+2
        END IF
        IF (I.LT.N-2) THEN
          MA=MA+1
          JA(MA)=I+3
        END IF
        IA(I+1)=MA+1
  200 CONTINUE
      XMAX=1.0D1
      RETURN
  210 IF (N.LT.7) GO TO 830
      NA=N
      MA=0
      IA(1)=1
      DO 220 I=1,N
        X(I)=-1.0D0
        IF (I.GT.1.AND.I.LT.N-3) THEN
          MA=MA+1
          JA(MA)=I-1
        END IF
        IF (I.LT.N-4) THEN
          MA=MA+1
          JA(MA)=I
        END IF
        IF (I.LT.N-5) THEN
          MA=MA+1
          JA(MA)=I+1
        END IF
        MA=MA+5
        JA(MA-4)=N-4
        JA(MA-3)=N-3
        JA(MA-2)=N-2
        JA(MA-1)=N-1
        JA(MA)=N
        IA(I+1)=MA+1
  220 CONTINUE
      RETURN
  230 IF (N.LT.2) GO TO 830
      N=N-MOD(N,2)
      NA=N
      MA=0
      IA(1)=1
      DO 240 I=1,N
        IF (MOD(I,2).EQ.1) THEN
          X(I)=9.0D1
          MA=MA+2
          JA(MA-1)=I
          JA(MA)=I+1
        ELSE
          X(I)=6.0D1
          MA=MA+2
          JA(MA-1)=I-1
          JA(MA)=I
        END IF
        IA(I+1)=MA+1
  240 CONTINUE
      RETURN
  250 IF (N.LT.4) GO TO 830
      N=N-MOD(N,4)
      NA=N
      MA=0
      IA(1)=1
      DO 260 I=1,N
        IF (MOD(I,4).EQ.1) THEN
          X(I)=3.0D0
          MA=MA+2
          JA(MA-1)=I
          JA(MA)=I+1
        ELSE IF (MOD(I,4).EQ.2) THEN
          X(I)=-1.0D0
          MA=MA+2
          JA(MA-1)=I+1
          JA(MA)=I+2
        ELSE IF (MOD(I,4).EQ.3) THEN
          X(I)=0.0D0
          MA=MA+2
          JA(MA-1)=I-1
          JA(MA)=I
        ELSE
          X(I)=1.0D0
          MA=MA+2
          JA(MA-1)=I-3
          JA(MA)=I
        END IF
        IA(I+1)=MA+1
  260 CONTINUE
      XMAX=1.0D2
      RETURN
  270 IF (N.LT.4) GO TO 830
      N=N-MOD(N,4)
      NA=N
      MA=0
      IA(1)=1
      DO 280 I=1,N
        IF (MOD(I,4).EQ.1) THEN
          X(I)=1.0D0
          MA=MA+2
          JA(MA-1)=I
          JA(MA)=I+1
        ELSE IF (MOD(I,4).EQ.2) THEN
          X(I)=2.0D0
          MA=MA+2
          JA(MA-1)=I
          JA(MA)=I+1
        ELSE IF (MOD(I,4).EQ.3) THEN
          X(I)=2.0D0
          MA=MA+2
          JA(MA-1)=I
          JA(MA)=I+1
        ELSE
          X(I)=2.0D0
          MA=MA+1
          JA(MA)=I
        END IF
        IA(I+1)=MA+1
  280 CONTINUE
      XMAX=1.0D1
      RETURN
  290 IF (N.LT.3) GO TO 830
      DO 300 I=1,N
        X(I)=-1.0D0
  300 CONTINUE
      GO TO 90
  310 IF (N.LT.7) GO TO 830
      NA=N
      MA=0
      IA(1)=1
      DO 330 I=1,5
        X(I)=-1.0D0
        DO 320 K=1,I+1
          MA=MA+1
          JA(MA)=K
  320   CONTINUE
        IA(I+1)=MA+1
  330 CONTINUE
      DO 350 I=6,N-1
        X(I)=-1.0D0
        DO 340 K=I-5,I+1
          MA=MA+1
          JA(MA)=K
  340   CONTINUE
        IA(I+1)=MA+1
  350 CONTINUE
      X(N)=-1.0D0
      DO 360 K=N-5,N
        MA=MA+1
        JA(MA)=K
  360 CONTINUE
      IA(N+1)=MA+1
      RETURN
  370 IF (N.LT.2) GO TO 830
      N=N-MOD(N,2)
      NA=N
      DO 380 I=1,N
        IF (MOD(I,2).EQ.1) THEN
          X(I)=0.0D0
        ELSE
          X(I)=1.0D0
        END IF
  380 CONTINUE
      MA=1
      DO 400 I=1,N-1,2
        IA(I)=2*I-1
        IA(I+1)=2*I+1
        DO 390 J=1,2
          JA(MA)=I
          MA=MA+1
          JA(MA)=I+1
          MA=MA+1
  390   CONTINUE
  400 CONTINUE
      IA(N+1)=MA
      MA=MA-1
      RETURN
  410 IF (N.LT.4) GO TO 830
      N=N-MOD(N,4)
      NA=N
      IA(1)=1
      J=0
      DO 420 I=2,N,2
        X(I-1)=-3.0D0
        X(I)=-1.0D0
        IA(I)=IA(I-1)+2
        IA(I+1)=IA(I)+3
        MA=IA(I-1)+4
        JA(MA-4)=I-1
        JA(MA-3)=I
        IF (J.EQ.0) THEN
          JA(MA-2)=I-1
          JA(MA-1)=I
          JA(MA)=I+2
          J=1
        ELSE
          JA(MA-2)=I-2
          JA(MA-1)=I-1
          JA(MA)=I
          J=0
        END IF
  420 CONTINUE
      RETURN
  430 IF (N.LT.3) GO TO 830
      DO 440 I=1,N
        X(I)=DBLE(I)/DBLE(N+1)
        X(I)=X(I)*(X(I)-1.0D0)
  440 CONTINUE
      GO TO 90
  450 IF (N.LT.3) GO TO 830
      DO 460 I=1,N
        X(I)=1.0D1
  460 CONTINUE
      GO TO 90
  470 IF (N.LT.16) GO TO 830
      M=INT(SQRT(DBLE(N)))
      PAR=6.7D0/DBLE(M+1)**2
      N=M*M
      K=0
      DO 490 J=1,M
        DO 480 I=1,M
          K=K+1
          X(K)=0.0D0
  480   CONTINUE
  490 CONTINUE
      GO TO 750
  500 IF (N.LT.16) GO TO 830
      M=INT(SQRT(DBLE(N)))
      N=M*M
      PAR=500.0D0
      K=0
      DO 520 J=1,M
        DO 510 I=1,M
          K=K+1
          X(K)=0.0D0
  510   CONTINUE
  520 CONTINUE
  530 NA=N
      MA=0
      IA(1)=1
      K=0
      DO 550 J=1,M
        DO 540 I=1,M
          K=K+1
          IF (J.GT.2) THEN
            MA=MA+1
            JA(MA)=K-M-M
          END IF
          IF (J.GT.1) THEN
            IF (I.GT.1) THEN
              MA=MA+1
              JA(MA)=K-M-1
            END IF
            MA=MA+1
            JA(MA)=K-M
            IF (I.LT.M) THEN
              MA=MA+1
              JA(MA)=K-M+1
            END IF
          END IF
          IF (I.GT.1) THEN
            IF (I.GT.2) THEN
              MA=MA+1
              JA(MA)=K-2
            END IF
            MA=MA+1
            JA(MA)=K-1
          END IF
          MA=MA+1
          JA(MA)=K
          IF (I.LT.M) THEN
            MA=MA+1
            JA(MA)=K+1
            IF (I.LT.M-1) THEN
              MA=MA+1
              JA(MA)=K+2
            END IF
          END IF
          IF (J.LT.M) THEN
            IF (I.GT.1) THEN
              MA=MA+1
              JA(MA)=K+M-1
            END IF
            MA=MA+1
            JA(MA)=K+M
            IF (I.LT.M) THEN
              MA=MA+1
              JA(MA)=K+M+1
            END IF
          END IF
          IF (J.LT.M-1) THEN
            MA=MA+1
            JA(MA)=K+M+M
          END IF
          IA(K+1)=MA+1
  540   CONTINUE
  550 CONTINUE
      RETURN
  560 IF (N.LT.3) GO TO 830
      DO 570 I=1,N
        X(I)=1.0D0
  570 CONTINUE
      PAR=1.0D1
      GO TO 90
  580 IF (N.LT.3) GO TO 830
      DO 590 I=1,N
        X(I)=1.5D0
  590 CONTINUE
      XMAX=1.0D0
      GO TO 90
  600 IF (N.LT.16) GO TO 830
      M=INT(SQRT(DBLE(N)))
      N=M*M
      PAR=500.0D0/DBLE(M+2)**4
      K=0
      DO 620 J=1,M
        DO 610 I=1,M
          K=K+1
          X(K)=0.0D0
  610   CONTINUE
  620 CONTINUE
      GO TO 530
  630 IF (N.LT.16) GO TO 830
      M=INT(SQRT(DBLE(N)))
      PAR=5.0D1/DBLE(M+1)
      N=M*M
      K=0
      DO 650 J=1,M
        DO 640 I=1,M
          K=K+1
          X(K)=1.0D0-DBLE(I)*DBLE(J)/DBLE(M+1)**2
  640   CONTINUE
  650 CONTINUE
      GO TO 750
  660 IF (N.LT.16) GO TO 830
      M=INT(SQRT(DBLE(N)))
      PAR=1.0D0/DBLE(M+1)**2
      N=M*M
      K=0
      DO 680 J=1,M
        DO 670 I=1,M
          K=K+1
          X(K)=-1.0D0
  670   CONTINUE
  680 CONTINUE
      GO TO 750
  690 IF (N.LT.16) GO TO 830
      M=INT(SQRT(DBLE(N)))
      PAR=1.0D0/DBLE(M+1)
      N=M*M
      K=0
      DO 710 J=1,M
        DO 700 I=1,M
          K=K+1
          X(K)=0.0D0
  700   CONTINUE
  710 CONTINUE
      GO TO 750
  720 IF (N.LT.16) GO TO 830
      M=INT(SQRT(DBLE(N)))
      PAR=1.0D0/DBLE(M+1)**2
      N=M*M
      K=0
      DO 740 J=1,M
        DO 730 I=1,M
          K=K+1
          X(K)=0.0D0
  730   CONTINUE
  740 CONTINUE
  750 NA=N
      MA=0
      IA(1)=1
      K=0
      DO 770 J=1,M
        DO 760 I=1,M
          K=K+1
          IF (J.GT.1) THEN
            MA=MA+1
            JA(MA)=K-M
          END IF
          IF (I.GT.1) THEN
            MA=MA+1
            JA(MA)=K-1
          END IF
          MA=MA+1
          JA(MA)=K
          IF (I.LT.M) THEN
            MA=MA+1
            JA(MA)=K+1
          END IF
          IF (J.LT.M) THEN
            MA=MA+1
            JA(MA)=K+M
          END IF
          IA(K+1)=MA+1
  760   CONTINUE
  770 CONTINUE
      RETURN
  780 IF (N.LT.5) GO TO 830
      PAR=5.0D2/DBLE(N+2)
      NA=N
      MA=0
      IA(1)=1
      DO 790 I=1,N
        X(I)=((DBLE(I)+0.5D0)/DBLE(N+2)-0.5D0)**2
        IF (I.GT.2) THEN
          MA=MA+1
          JA(MA)=I-2
        END IF
        IF (I.GT.1) THEN
          MA=MA+1
          JA(MA)=I-1
        END IF
        MA=MA+1
        JA(MA)=I
        IF (I.LT.N) THEN
          MA=MA+1
          JA(MA)=I+1
        END IF
        IF (I.LT.N-1) THEN
          MA=MA+1
          JA(MA)=I+2
        END IF
        IA(I+1)=MA+1
  790 CONTINUE
      RETURN
  800 IF (N.LT.10) GO TO 830
      N=N-MOD(N,2)
      M=N/2
      PAR=5.0D2
      NA=N
      MA=0
      IA(1)=1
      DO 810 I=1,M
        K=I+M
        X(I)=(DBLE(I)/DBLE(M+1)-0.5D0)**2
        IF (I.GT.2) THEN
          MA=MA+1
          JA(MA)=I-2
        END IF
        IF (I.GT.1) THEN
          MA=MA+1
          JA(MA)=I-1
        END IF
        MA=MA+1
        JA(MA)=I
        IF (I.LT.M) THEN
          MA=MA+1
          JA(MA)=I+1
        END IF
        IF (I.LT.M-1) THEN
          MA=MA+1
          JA(MA)=I+2
        END IF
        IF (I.GT.1) THEN
          MA=MA+1
          JA(MA)=K-1
        END IF
        MA=MA+1
        JA(MA)=K
        IF (I.LT.M) THEN
          MA=MA+1
          JA(MA)=K+1
        END IF
        IA(I+1)=MA+1
  810 CONTINUE
      DO 820 I=M+1,N
        K=I-M
        X(I)=DBLE(K)/DBLE(M+1)-0.5D0
        IF (K.GT.1) THEN
          MA=MA+1
          JA(MA)=K-1
        END IF
        MA=MA+1
        JA(MA)=K
        IF (K.LT.M) THEN
          MA=MA+1
          JA(MA)=K+1
        END IF
        IF (K.GT.1) THEN
          MA=MA+1
          JA(MA)=I-1
        END IF
        MA=MA+1
        JA(MA)=I
        IF (K.LT.M) THEN
          MA=MA+1
          JA(MA)=I+1
        END IF
        IA(I+1)=MA+1
  820 CONTINUE
      RETURN
  830 IERR=1
      RETURN
      END
! SUBROUTINE TAFU18             ALL SYSTEMS                 97/12/01
! PORTABILITY : ALL SYSTEMS
! 97/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
!  VALUES OF EQUATION FUNCTIONS (LEFT HAND SIDES) IN THE SYSTEM OF
!  NONLINEAR EQUATIONS.
!
! PARAMETERS :
!  II  N  NUMBER OF VARIABLES.
!  II  KA  INDEX OF THE GIVEN EQUATION.
!  RI  X(N)  VECTOR OF VARIABLES.
!  RO  FA  VALUE OF THE KA-TH EQUATION FUNCTION AT THE POINT X.
!  II  NEXT  NUMBER OF THE SELECTED TEST PROBLEM.
!
      SUBROUTINE TAFU18 (N, KA, X, FA, NEXT)
      INTEGER N,KA,NEXT
      DOUBLE PRECISION X(*),FA
      INTEGER I,J,M
      DOUBLE PRECISION A1,A2,A3,A4,H,PAR,PI
      DATA PI/3.14159265358979323D0/
      COMMON /EMPR18/ PAR,M
      GO TO (10,20,30,50,60,70,80,90,100,110,120,130,140,150,160,180,
     &190,200,210,220,230,240,250,260,270,280,290,300,310,320),NEXT
   10 A2=0.5D0
      IF (KA.EQ.1) THEN
        FA=A2-(1.0D0-A2)*X(3)-X(1)*(1.0D0+4.0D0*X(2))
      ELSE IF (KA.EQ.2) THEN
        FA=-(2.0D0-A2)*X(4)-X(2)*(1.0D0+4.0D0*X(1))
      ELSE IF (KA.EQ.N-1) THEN
        FA=A2*X(N-3)-X(N-1)*(1.0D0+4.0D0*X(N))
      ELSE IF (KA.EQ.N) THEN
        FA=A2*X(N-2)-(2.0D0-A2)-X(N)*(1.0D0+4.0D0*X(N-1))
      ELSE IF (MOD(KA,2).EQ.1) THEN
        FA=A2*X(KA-2)-(1.0D0-A2)*X(KA+2)-X(KA)*(1.0D0+4.0D0*X(KA+1))
      ELSE
        FA=A2*X(KA-2)-(2.0D0-A2)*X(KA+2)-X(KA)*(1.0D0+4.0D0*X(KA-1))
      END IF
      RETURN
   20 A1=0.414214D0
      IF (KA.EQ.1) THEN
        FA=X(1)-(1.0D0-X(1))*X(3)-A1*(1.0D0+4.0D0*X(2))
      ELSE IF (KA.EQ.2) THEN
        FA=-(1.0D0-X(1))*X(4)-A1*(1.0D0+4.0D0*X(2))
      ELSE IF (KA.EQ.3) THEN
        FA=A1*X(1)-(1.0D0-X(1))*X(5)-X(3)*(1.0D0+4.0D0*X(2))
      ELSE IF (KA.LE.N-2) THEN
        FA=X(1)*X(KA-2)-(1.0D0-X(1))*X(KA+2)-X(KA)*(1.0D0+4.0D0*X(KA-1))
      ELSE IF (KA.EQ.N-1) THEN
        FA=X(1)*X(N-3)-X(N-1)*(1.0D0+4.0D0*X(N-2))
      ELSE
        FA=X(1)*X(N-2)-(1.0D0-X(1))-X(N)*(1.0D0+4.0D0*X(N-1))
      END IF
      RETURN
   30 J=(KA-1)/5
      FA=5.0D0-DBLE(J+1)*(1.0D0-COS(X(KA)))-SIN(X(KA))
      J=J*5
      DO 40 I=J+1,J+5
        FA=FA-COS(X(I))
   40 CONTINUE
      RETURN
   50 IF (KA.LT.2) THEN
        FA=3.0D0*X(KA)**3+2.0D0*X(KA+1)-5.0D0+SIN(X(KA)-X(KA+1))*
     &   SIN(X(KA)+X(KA+1))
      ELSE IF (KA.LT.N) THEN
        FA=3.0D0*X(KA)**3+2.0D0*X(KA+1)-5.0D0+SIN(X(KA)-X(KA+1))*
     &   SIN(X(KA)+X(KA+1))+4.0D0*X(KA)-X(KA-1)*EXP(X(KA-1)-X(KA))-
     &   3.0D0
      ELSE
        FA=4.0D0*X(KA)-X(KA-1)*EXP(X(KA-1)-X(KA))-3.0D0
      END IF
      RETURN
   60 IF (MOD(KA,2).EQ.1) THEN
        FA=0.0D0
        IF (KA.NE.1) FA=FA-6.0D0*(X(KA-2)-X(KA))**3+1.0D1-4.0D0*X(KA-1)-
     &   2.0D0*SIN(X(KA-2)-X(KA-1)-X(KA))*SIN(X(KA-2)+X(KA-1)-X(KA))
        IF (KA.NE.N) FA=FA+3.0D0*(X(KA)-X(KA+2))**3-5.0D0+2.0D0*X(KA+1)+
     &   SIN(X(KA)-X(KA+1)-X(KA+2))*SIN(X(KA)+X(KA+1)-X(KA+2))
      ELSE
        FA=4.0D0*X(KA)-(X(KA-1)-X(KA+1))*EXP(X(KA-1)-X(KA)-X(KA+1))-
     &   3.0D0
      END IF
      RETURN
   70 H=2.0D0
      IF (KA.EQ.1) THEN
        FA=((3.0D0-H*X(1))*X(1)-2.0D0*X(2)+1.0D0)**2
      ELSE IF (KA.LE.N-1) THEN
        FA=((3.0D0-H*X(KA))*X(KA)-X(KA-1)-2.0D0*X(KA+1)+1.0D0)**2
      ELSE
        FA=((3.0D0-H*X(N))*X(N)-X(N-1)+1.0D0)**2
      END IF
      RETURN
   80 IF (KA.LT.2) THEN
        FA=4.0D0*(X(KA)-X(KA+1)**2)
      ELSE IF (KA.LT.N) THEN
        FA=8.0D0*X(KA)*(X(KA)**2-X(KA-1))-2.0D0*(1.0D0-X(KA))+4.0D0*
     &   (X(KA)-X(KA+1)**2)
      ELSE
        FA=8.0D0*X(KA)*(X(KA)**2-X(KA-1))-2.0D0*(1.0D0-X(KA))
      END IF
      RETURN
   90 IF (KA.LT.2) THEN
        FA=4.0D0*(X(KA)-X(KA+1)**2)+X(KA+1)-X(KA+2)**2
      ELSE IF (KA.LT.3) THEN
        FA=8.0D0*X(KA)*(X(KA)**2-X(KA-1))-2.0D0*(1.0D0-X(KA))+4.0D0*
     &   (X(KA)-X(KA+1)**2)+X(KA+1)-X(KA+2)**2
      ELSE IF (KA.LT.N-1) THEN
        FA=8.0D0*X(KA)*(X(KA)**2-X(KA-1))-2.0D0*(1.0D0-X(KA))+4.0D0*
     &   (X(KA)-X(KA+1)**2)+X(KA-1)**2-X(KA-2)+X(KA+1)-X(KA+2)**2
      ELSE IF (KA.LT.N) THEN
        FA=8.0D0*X(KA)*(X(KA)**2-X(KA-1))-2.0D0*(1.0D0-X(KA))+4.0D0*
     &   (X(KA)-X(KA+1)**2)+X(KA-1)**2-X(KA-2)
      ELSE
        FA=8.0D0*X(KA)*(X(KA)**2-X(KA-1))-2.0D0*(1.0D0-X(KA))+X(KA-1)**
     &   2-X(KA-2)
      END IF
      RETURN
  100 IF (KA.LT.2) THEN
        FA=4.0D0*(X(KA)-X(KA+1)**2)+X(KA+1)-X(KA+2)**2+X(KA+2)-X(KA+3)**
     &   2
      ELSE IF (KA.LT.3) THEN
        FA=8.0D0*X(KA)*(X(KA)**2-X(KA-1))-2.0D0*(1.0D0-X(KA))+4.0D0*
     &   (X(KA)-X(KA+1)**2)+X(KA-1)**2+X(KA+1)-X(KA+2)**2+X(KA+2)-X(KA+
     &   3)**2
      ELSE IF (KA.LT.4) THEN
        FA=8.0D0*X(KA)*(X(KA)**2-X(KA-1))-2.0D0*(1.0D0-X(KA))+4.0D0*
     &   (X(KA)-X(KA+1)**2)+X(KA-1)**2-X(KA-2)+X(KA+1)-X(KA+2)**2+X(KA-
     &   2)**2+X(KA+2)-X(KA+3)**2
      ELSE IF (KA.LT.N-2) THEN
        FA=8.0D0*X(KA)*(X(KA)**2-X(KA-1))-2.0D0*(1.0D0-X(KA))+4.0D0*
     &   (X(KA)-X(KA+1)**2)+X(KA-1)**2-X(KA-2)+X(KA+1)-X(KA+2)**2+X(KA-
     &   2)**2-X(KA-3)+X(KA+2)-X(KA+3)**2
      ELSE IF (KA.LT.N-1) THEN
        FA=8.0D0*X(KA)*(X(KA)**2-X(KA-1))-2.0D0*(1.0D0-X(KA))+4.0D0*
     &   (X(KA)-X(KA+1)**2)+X(KA-1)**2-X(KA-2)+X(KA+1)-X(KA+2)**2+X(KA-
     &   2)**2-X(KA-3)+X(KA+2)
      ELSE IF (KA.LT.N) THEN
        FA=8.0D0*X(KA)*(X(KA)**2-X(KA-1))-2.0D0*(1.0D0-X(KA))+4.0D0*
     &   (X(KA)-X(KA+1)**2)+X(KA-1)**2-X(KA-2)+X(KA+1)+X(KA-2)**2-X(KA-
     &   3)
      ELSE
        FA=8.0D0*X(KA)*(X(KA)**2-X(KA-1))-2.0D0*(1.0D0-X(KA))+X(KA-1)**
     &   2-X(KA-2)+X(KA-2)**2-X(KA-3)
      END IF
      RETURN
  110 IF (KA.EQ.1) THEN
        FA=-2.0D0*X(KA)**2+3.0D0*X(KA)-2.0D0*X(KA+1)+3.0D0*X(N-4)-X(N-3)
     &   -X(N-2)+0.5D0*X(N-1)-X(N)+1.0D0
      ELSE IF (KA.LE.N-1) THEN
        FA=-2.0D0*X(KA)**2+3.0D0*X(KA)-X(KA-1)-2.0D0*X(KA+1)+3.0D0*X(N-
     &   4)-X(N-3)-X(N-2)+0.5D0*X(N-1)-X(N)+1.0D0
      ELSE
        FA=-2.0D0*X(N)**2+3.0D0*X(N)-X(N-1)+3.0D0*X(N-4)-X(N-3)-X(N-2)+
     &   0.5D0*X(N-1)-X(N)+1.0D0
      END IF
      RETURN
  120 IF (MOD(KA,2).EQ.1) THEN
        FA=X(KA)+((5.0D0-X(KA+1))*X(KA+1)-2.0D0)*X(KA+1)-1.3D1
      ELSE
        FA=X(KA-1)+((X(KA)+1.0D0)*X(KA)-1.4D1)*X(KA)-2.9D1
      END IF
      RETURN
  130 IF (MOD(KA,4).EQ.1) THEN
        FA=X(KA)+1.0D1*X(KA+1)
      ELSE IF (MOD(KA,4).EQ.2) THEN
        FA=2.23606797749979D0*(X(KA+1)-X(KA+2))
      ELSE IF (MOD(KA,4).EQ.3) THEN
        FA=(X(KA-1)-2.0D0*X(KA))**2
      ELSE
        FA=3.16227766016838D0*(X(KA-3)-X(KA))**2
      END IF
      RETURN
  140 IF (MOD(KA,4).EQ.1) THEN
        FA=(EXP(X(KA))-X(KA+1))**2
      ELSE IF (MOD(KA,4).EQ.2) THEN
        FA=1.0D1*(X(KA)-X(KA+1))**3
      ELSE IF (MOD(KA,4).EQ.3) THEN
        FA=X(KA)-X(KA+1)
        FA=(SIN(FA)/COS(FA))**2
      ELSE
        FA=X(KA)-1.0D0
      END IF
      RETURN
  150 IF (KA.LT.2) THEN
        FA=X(KA)*(0.5D0*X(KA)-3.0D0)-1.0D0+2.0D0*X(KA+1)
      ELSE IF (KA.LT.N) THEN
        FA=X(KA-1)+X(KA)*(0.5D0*X(KA)-3.0D0)-1.0D0+2.0D0*X(KA+1)
      ELSE
        FA=X(KA-1)+X(KA)*(0.5D0*X(KA)-3.0D0)-1.0D0
      END IF
      RETURN
  160 FA=(2.0D0+5.0D0*X(KA)**2)*X(KA)+1.0D0
      DO 170 I=MAX(1,KA-5),MIN(N,KA+1)
        FA=FA+X(I)*(1.0D0+X(I))
  170 CONTINUE
      RETURN
  180 IF (MOD(KA,2).EQ.1) THEN
        FA=1.0D4*X(KA)*X(KA+1)-1.0D0
      ELSE
        FA=EXP(-X(KA-1))+EXP(-X(KA))-1.0001D0
      END IF
      RETURN
  190 IF (MOD(KA,4).EQ.1) THEN
        FA=-2.0D2*X(KA)*(X(KA+1)-X(KA)**2)-(1.0D0-X(KA))
      ELSE IF (MOD(KA,4).EQ.2) THEN
        FA=2.0D2*(X(KA)-X(KA-1)**2)+2.02D1*(X(KA)-1.0D0)+1.98D1*(X(KA+2)
     &   -1.0D0)
      ELSE IF (MOD(KA,4).EQ.3) THEN
        FA=-1.8D2*X(KA)*(X(KA+1)-X(KA)**2)-(1.0D0-X(KA))
      ELSE
        FA=1.8D2*(X(KA)-X(KA-1)**2)+2.02D1*(X(KA)-1.0D0)+1.98D1*(X(KA-2)
     &   -1.0D0)
      END IF
      RETURN
  200 IF (KA.LT.2) THEN
        FA=X(KA)-EXP(COS(DBLE(KA)*(X(KA)+X(KA+1))))
      ELSE IF (KA.LT.N) THEN
        FA=X(KA)-EXP(COS(DBLE(KA)*(X(KA-1)+X(KA)+X(KA+1))))
      ELSE
        FA=X(KA)-EXP(COS(DBLE(KA)*(X(KA-1)+X(KA))))
      END IF
      RETURN
  210 A3=1.0D0/DBLE(N+1)
      A4=DBLE(KA)*A3
      FA=2.0D0*X(KA)+0.5D0*A3*A3*(X(KA)+A4+1.0D0)**3
      IF (KA.GT.1) FA=FA-X(KA-1)
      IF (KA.LT.N) FA=FA-X(KA+1)
      RETURN
  220 IF (KA.EQ.1) THEN
        FA=3.0D0*X(KA)*(X(KA+1)-2.0D0*X(KA))+0.25D0*X(KA+1)**2
      ELSE IF (KA.EQ.N) THEN
        FA=3.0D0*X(KA)*(2.0D1-2.0D0*X(KA)+X(KA-1))+0.25D0*(2.0D1-X(KA-1)
     &   )**2
      ELSE
        FA=3.0D0*X(KA)*(X(KA+1)-2.0D0*X(KA)+X(KA-1))+0.25D0*(X(KA+1)-
     &   X(KA-1))**2
      END IF
      RETURN
  230 H=1.0D0/DBLE(N+1)
      IF (KA.LT.2) THEN
        FA=2.0D0*X(KA)+PAR*H**2*SINH(PAR*X(KA))-X(KA+1)
      ELSE IF (KA.LT.N) THEN
        FA=2.0D0*X(KA)+PAR*H**2*SINH(PAR*X(KA))-X(KA-1)-X(KA+1)
      ELSE
        FA=2.0D0*X(KA)+PAR*H**2*SINH(PAR*X(KA))-X(KA-1)-1.0D0
      END IF
      RETURN
  240 FA=6.0D0*X(KA)
      A1=0.0D0
      A2=0.0D0
      A3=0.0D0
      IF (KA.GT.1) THEN
        FA=FA-4.0D0*X(KA-1)
        A1=A1-X(KA-1)
        A2=A2+X(KA-1)
        A3=A3+2.0D0*X(KA-1)
      END IF
      IF (KA.GT.2) THEN
        FA=FA+X(KA-2)
        A3=A3-X(KA-2)
      END IF
      IF (KA.LT.N-1) THEN
        FA=FA+X(KA+2)
        A3=A3+X(KA+2)
      END IF
      IF (KA.LT.N) THEN
        FA=FA-4.0D0*X(KA+1)
        A1=A1+X(KA+1)
        A2=A2+X(KA+1)
        A3=A3-2.0D0*X(KA+1)
      END IF
      IF (KA.GE.N-1) THEN
        FA=FA+1.0D0
        A3=A3+1.0D0
      END IF
      IF (KA.GE.N) THEN
        FA=FA-4.0D0
        A1=A1+1.0D0
        A2=A2+1.0D0
        A3=A3-2.0D0
      END IF
      FA=FA-0.5D0*PAR*(A1*A2-X(KA)*A3)
      RETURN
  250 H=1.0D0/(M+1)
      IF (KA.LE.M) THEN
        J=KA+M
        FA=6.0D0*X(KA)
        A1=0.0D0
        A2=0.0D0
        IF (KA.EQ.1) THEN
          A1=A1+1.0D0
        END IF
        IF (KA.GT.1) THEN
          FA=FA-4.0D0*X(KA-1)
          A1=A1-X(J-1)
          A2=A2+2.0D0*X(KA-1)
        END IF
        IF (KA.GT.2) THEN
          FA=FA+X(KA-2)
          A2=A2-X(KA-2)
        END IF
        IF (KA.LT.M-1) THEN
          FA=FA+X(KA+2)
          A2=A2+X(KA+2)
        END IF
        IF (KA.LT.M) THEN
          FA=FA-4.0D0*X(KA+1)
          A1=A1+X(J+1)
          A2=A2-2.0D0*X(KA+1)
        END IF
        IF (KA.EQ.M) THEN
          A1=A1+1.0D0
        END IF
        FA=FA+0.5D0*PAR*H*(X(KA)*A2+X(J)*A1*H**2)
      ELSE
        J=KA-M
        FA=-2.0D0*X(KA)
        A1=0.0D0
        A2=0.0D0
        IF (J.EQ.1) THEN
          A2=A2+1.0D0
        END IF
        IF (J.GT.1) THEN
          FA=FA+X(KA-1)
          A1=A1-X(J-1)
          A2=A2-X(KA-1)
        END IF
        IF (J.LT.M) THEN
          FA=FA+X(KA+1)
          A1=A1+X(J+1)
          A2=A2+X(KA+1)
        END IF
        IF (J.EQ.M) THEN
          A2=A2+1.0D0
        END IF
        FA=FA+0.5D0*PAR*H*(X(KA)*A1+X(J)*A2)
      END IF
      RETURN
  260 FA=4.0D0*X(KA)-PAR*EXP(X(KA))
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      IF (I.GT.1) FA=FA-X(KA-1)
      IF (I.LT.M) FA=FA-X(KA+1)
      IF (J.GT.1) FA=FA-X(KA-M)
      IF (J.LT.M) FA=FA-X(KA+M)
      RETURN
  270 FA=4.0D0*X(KA)
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      FA=FA+PAR*X(KA)**3/(1.0D0+PAR*DBLE(I)**2+PAR*DBLE(J)**2)
      IF (I.EQ.1) FA=FA-1.0D0
      IF (I.GT.1) FA=FA-X(KA-1)
      IF (I.LT.M) FA=FA-X(KA+1)
      IF (I.EQ.M) FA=FA-2.0D0+EXP(DBLE(J)/DBLE(M+1))
      IF (J.EQ.1) FA=FA-1.0D0
      IF (J.GT.1) FA=FA-X(KA-M)
      IF (J.LT.M) FA=FA-X(KA+M)
      IF (J.EQ.M) FA=FA-2.0D0+EXP(DBLE(I)/DBLE(M+1))
      RETURN
  280 FA=4.0D0*X(KA)-PAR*SIN(2.0D0*PI*X(KA))
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      A1=DBLE(I)/DBLE(M+1)
      A2=DBLE(J)/DBLE(M+1)
      FA=FA-1.0D4*((A1-0.25D0)**2+(A2-0.75D0)**2)*PAR
      IF (I.EQ.1) FA=FA-X(KA+1)-PAR*SIN(PI*X(KA+1)*DBLE(M+1))
      IF (I.GT.1.AND.I.LT.M) FA=FA-X(KA+1)-X(KA-1)-PAR*SIN(PI*(X(KA+1)-
     &X(KA-1))*DBLE(M+1))
      IF (I.EQ.M) FA=FA-X(KA-1)+PAR*SIN(PI*X(KA-1)*DBLE(M+1))
      IF (J.EQ.1) FA=FA-X(KA+M)-PAR*SIN(PI*X(KA+M)*DBLE(M+1))
      IF (J.GT.1.AND.J.LT.M) FA=FA-X(KA+M)-X(KA-M)-PAR*SIN(PI*(X(KA+M)-
     &X(KA-M))*DBLE(M+1))
      IF (J.EQ.M) FA=FA-X(KA-M)+PAR*SIN(PI*X(KA-M)*DBLE(M+1))
      RETURN
  290 FA=8.0D0*X(KA)**2
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      IF (I.EQ.1) FA=FA-2.0D0*X(KA)*(X(KA+1)+1.0D0)-0.5D0*(X(KA+1)-
     &1.0D0)**2-1.5D0*X(KA)**2*(X(KA+1)-1.0D0)*PAR
      IF (I.GT.1.AND.I.LT.M) FA=FA-2.0D0*X(KA)*(X(KA+1)+X(KA-1))-0.5D0*
     &(X(KA+1)-X(KA-1))**2-1.5D0*X(KA)**2*(X(KA+1)-X(KA-1))*PAR
      IF (I.EQ.M) FA=FA-2.0D0*X(KA)*X(KA-1)-0.5D0*X(KA-1)**2+1.5D0*X(KA)
     &**2*X(KA-1)*PAR
      IF (J.EQ.1) FA=FA-2.0D0*X(KA)*(X(KA+M)+1.0D0)-0.5D0*(X(KA+M)-
     &1.0D0)**2
      IF (J.GT.1.AND.J.LT.M) FA=FA-2.0D0*X(KA)*(X(KA+M)+X(KA-M))-0.5D0*
     &(X(KA+M)-X(KA-M))**2
      IF (J.EQ.M) FA=FA-2.0D0*X(KA)*X(KA-M)-0.5D0*X(KA-M)**2
      IF (I.EQ.1.AND.J.EQ.1) FA=FA-PAR/DBLE(M+1)
      RETURN
  300 FA=4.0D0*X(KA)
      A3=0.0D0
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      A1=PAR*DBLE(I)
      A2=PAR*DBLE(J)
      FA=FA-2.0D3*A1*A2*(1.0D0-A1)*(1.0D0-A2)*PAR**2
      IF (I.GT.1) THEN
        FA=FA-X(KA-1)
        A3=A3-X(KA-1)
      END IF
      IF (I.LT.M) THEN
        FA=FA-X(KA+1)
        A3=A3+X(KA+1)
      END IF
      IF (J.GT.1) THEN
        FA=FA-X(KA-M)
        A3=A3-X(KA-M)
      END IF
      IF (J.LT.M) THEN
        FA=FA-X(KA+M)
        A3=A3+X(KA+M)
      END IF
      FA=FA+2.0D1*PAR*A3*X(KA)
      RETURN
  310 FA=2.0D1*X(KA)-PAR*MAX(0.0D0,X(KA))
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      FA=FA-SIGN(PAR,(DBLE(I)/DBLE(M+2)-0.5D0))
      IF (J.GT.2) THEN
        FA=FA+X(KA-M-M)
      END IF
      IF (J.GT.1) THEN
        IF (I.GT.1) THEN
          FA=FA+2.0D0*X(KA-M-1)
        END IF
        FA=FA-8.0D0*X(KA-M)
        IF (I.LT.M) THEN
          FA=FA+2.0D0*X(KA-M+1)
        END IF
      END IF
      IF (I.GT.1) THEN
        IF (I.GT.2) THEN
          FA=FA+X(KA-2)
        END IF
        FA=FA-8.0D0*X(KA-1)
      END IF
      IF (I.LT.M) THEN
        FA=FA-8.0D0*X(KA+1)
        IF (I.LT.M-1) THEN
          FA=FA+X(KA+2)
        END IF
      END IF
      IF (J.LT.M) THEN
        IF (I.GT.1) THEN
          FA=FA+2.0D0*X(KA+M-1)
        END IF
        FA=FA-8.0D0*X(KA+M)
        IF (I.LT.M) THEN
          FA=FA+2.0D0*X(KA+M+1)
        END IF
      END IF
      IF (J.LT.M-1) THEN
        FA=FA+X(KA+M+M)
      END IF
      RETURN
  320 H=0.5D0/DBLE(M+2)
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      FA=2.0D1*X(KA)
      A1=0.0D0
      A2=0.0D0
      A3=0.0D0
      A4=0.0D0
      IF (J.GT.2) THEN
        FA=FA+X(KA-M-M)
        A4=A4+X(KA-M-M)
      END IF
      IF (J.GT.1) THEN
        IF (I.GT.1) THEN
          FA=FA+2.0D0*X(KA-M-1)
          A3=A3+X(KA-M-1)
          A4=A4+X(KA-M-1)
        END IF
        FA=FA-8.0D0*X(KA-M)
        A1=A1-X(KA-M)
        A4=A4-4.0D0*X(KA-M)
        IF (I.LT.M) THEN
          FA=FA+2.0D0*X(KA-M+1)
          A3=A3-X(KA-M+1)
          A4=A4+X(KA-M+1)
        END IF
      END IF
      IF (I.GT.1) THEN
        IF (I.GT.2) THEN
          FA=FA+X(KA-2)
          A3=A3+X(KA-2)
        END IF
        FA=FA-8.0D0*X(KA-1)
        A2=A2-X(KA-1)
        A3=A3-4.0D0*X(KA-1)
      END IF
      IF (I.LT.M) THEN
        FA=FA-8.0D0*X(KA+1)
        A2=A2+X(KA+1)
        A3=A3+4.0D0*X(KA+1)
        IF (I.LT.M-1) THEN
          FA=FA+X(KA+2)
          A3=A3-X(KA+2)
        END IF
      END IF
      IF (J.LT.M) THEN
        IF (I.GT.1) THEN
          FA=FA+2.0D0*X(KA+M-1)
          A3=A3+X(KA+M-1)
          A4=A4-X(KA+M-1)
        END IF
        FA=FA-8.0D0*X(KA+M)
        A1=A1+X(KA+M)
        A4=A4+4.0D0*X(KA+M)
        IF (I.LT.M) THEN
          FA=FA+2.0D0*X(KA+M+1)
          A3=A3-X(KA+M+1)
          A4=A4-X(KA+M+1)
        END IF
      END IF
      IF (J.LT.M-1) THEN
        FA=FA+X(KA+M+M)
        A4=A4-X(KA+M+M)
      END IF
      IF (J.EQ.M) THEN
        IF (I.GT.1) THEN
          FA=FA-H-H
          A3=A3-H
          A4=A4+H
        END IF
        FA=FA+8.0D0*H
        A1=A1-H
        A4=A4-4.0D0*H
        IF (I.LT.M) THEN
          FA=FA-2.0D0*H
          A3=A3+H
          A4=A4+H
        END IF
        FA=FA+H
        A4=A4-H
      END IF
      IF (J.EQ.M-1) THEN
        FA=FA-H
        A4=A4+H
      END IF
      FA=FA+0.25D0*PAR*(A1*A3-A2*A4)
      RETURN
      END
! SUBROUTINE TAGU18             ALL SYSTEMS                 92/12/01
! PORTABILITY : ALL SYSTEMS
! 92/12/01 LU : ORIGINAL VERSION
!
! PURPOSE :
!  GRADIENTS OF TEST FUNCTIONS FOR NONLINEAR EQUATIONS.
!  UNIVERSAL VERSION.
!
! PARAMETERS :
!  II  N  NUMBER OF VARIABLES.
!  II  KA  INDEX OF THE APPROXIMATED FUNCTION.
!  RI  X(N)  VECTOR OF VARIABLES.
!  RO  GA(N)  GRADIENT OF THE APPROXIMATED FUNCTION AT THE
!          SELECTED POINT.
!  II  NEXT  NUMBER OF THE TEST PROBLEM.
!
      SUBROUTINE TAGU18 (N, KA, X, GA, NEXT)
      INTEGER I,J,K,M,N,KA,NEXT
      DOUBLE PRECISION X(*),GA(*)
      DOUBLE PRECISION GA1(2),GA2(2),GA3(6),GA4(6)
      DOUBLE PRECISION U,V,W,EX,PI,D1S,D2S,H,ALFA,A1,A2,A3,A4,PAR
      DATA PI/3.14159265358979323D0/
      COMMON /EMPR18/ PAR,M
      GO TO (10,20,30,50,60,70,80,90,100,110,120,130,140,150,160,180,
     &190,200,210,220,230,240,250,260,270,280,290,300,310,320,350,340),
     &NEXT
   10 ALFA=0.5D0
      IF (KA.EQ.1) THEN
        GA(1)=-1.0D0-4.0D0*X(2)
        GA(2)=-4.0D0*X(1)
        GA(3)=ALFA-1.0D0
      ELSE IF (KA.EQ.2) THEN
        GA(1)=-4.0D0*X(2)
        GA(2)=-1.0D0-4.0D0*X(1)
        GA(4)=-2.0D0+ALFA
      ELSE IF (KA.EQ.N-1) THEN
        GA(N-3)=ALFA
        GA(N-1)=-1.0D0-4.0D0*X(N)
        GA(N)=-4.0D0*X(N-1)
      ELSE IF (KA.EQ.N) THEN
        GA(N-2)=ALFA
        GA(N-1)=-4.0D0*X(N)
        GA(N)=-1.0D0-4.0D0*X(N-1)
      ELSE IF (MOD(KA,2).EQ.1) THEN
        GA(KA-2)=ALFA
        GA(KA)=-1.0D0-4.0D0*X(KA+1)
        GA(KA+1)=-4.0D0*X(KA)
        GA(KA+2)=-1.0D0+ALFA
      ELSE
        GA(KA-2)=ALFA
        GA(KA-1)=-4.0D0*X(KA)
        GA(KA)=-1.0D0-4.0D0*X(KA-1)
        GA(KA+2)=-2.0D0+ALFA
      END IF
      RETURN
   20 A1=0.414214D0
      IF (KA.EQ.1) THEN
        GA(1)=1.0D0+X(3)
        GA(2)=-4.0D0*A1
        GA(3)=-1.0D0+X(1)
      ELSE IF (KA.EQ.2) THEN
        GA(1)=X(4)
        GA(2)=-4.0D0*A1
        GA(4)=-1.0D0+X(1)
      ELSE IF (KA.EQ.3) THEN
        GA(1)=A1+X(5)
        GA(2)=-4.0D0*X(3)
        GA(3)=-1.0D0-4.0D0*X(2)
        GA(5)=-1.0D0+X(1)
      ELSE IF (KA.LE.N-2) THEN
        GA(1)=X(KA-2)+X(KA+2)
        GA(KA-2)=X(1)
        GA(KA-1)=-4.0D0*X(KA)
        GA(KA)=-1.0D0-4.0D0*X(KA-1)
        GA(KA+2)=-1.0D0+X(1)
      ELSE IF (KA.EQ.N-1) THEN
        GA(1)=X(N-3)
        GA(N-3)=X(1)
        GA(N-2)=-4.0D0*X(N-1)
        GA(N-1)=-1.0D0-4.0D0*X(N-2)
      ELSE
        GA(1)=X(N-2)+1.0D0
        GA(N-2)=X(1)
        GA(N-1)=-4.0D0*X(N)
        GA(N)=-1.0D0-4.0D0*X(N-1)
      END IF
      RETURN
   30 J=(KA-1)/5
      GA(KA)=-DBLE(J+1)*SIN(X(KA))-COS(X(KA))
      J=J*5
      DO 40 I=J+1,J+5
        IF (I.EQ.KA) THEN
          GA(I)=GA(I)+SIN(X(I))
        ELSE
          GA(I)=SIN(X(I))
        END IF
   40 CONTINUE
      RETURN
   50 IF (KA.LT.2) THEN
        D1S=COS(X(KA)-X(KA+1))*SIN(X(KA)+X(KA+1))
        D2S=SIN(X(KA)-X(KA+1))*COS(X(KA)+X(KA+1))
        GA(KA)=9.0D0*X(KA)**2+D1S+D2S
        GA(KA+1)=2.0D0-D1S+D2S
      ELSE IF (KA.LT.N) THEN
        D1S=COS(X(KA)-X(KA+1))*SIN(X(KA)+X(KA+1))
        D2S=SIN(X(KA)-X(KA+1))*COS(X(KA)+X(KA+1))
        EX=EXP(X(KA-1)-X(KA))
        GA(KA-1)=-EX-X(KA-1)*EX
        GA(KA)=9.0D0*X(KA)**2+D1S+D2S+4.0D0+X(KA-1)*EX
        GA(KA+1)=2.0D0-D1S+D2S
      ELSE
        EX=EXP(X(KA-1)-X(KA))
        GA(KA-1)=-EX-X(KA-1)*EX
        GA(KA)=4.0D0+X(KA-1)*EX
      END IF
      RETURN
   60 IF (MOD(KA,2).EQ.1) THEN
        GA(KA)=0.0D0
        IF (KA.NE.1) THEN
          D1S=COS(X(KA-2)-X(KA-1)-X(KA))*SIN(X(KA-2)+X(KA-1)-X(KA))
          D2S=SIN(X(KA-2)-X(KA-1)-X(KA))*COS(X(KA-2)+X(KA-1)-X(KA))
          GA(KA-2)=-18.0D0*(X(KA-2)-X(KA))**2-2.0D0*(D1S+D2S)
          GA(KA-1)=-4.0D0+2.0D0*(D1S-D2S)
          GA(KA)=GA(KA)+18.0D0*(X(KA-2)-X(KA))**2+2.0D0*(D1S+D2S)
        END IF
        IF (KA.NE.N) THEN
          D1S=COS(X(KA)-X(KA+1)-X(KA+2))*SIN(X(KA)+X(KA+1)-X(KA+2))
          D2S=SIN(X(KA)-X(KA+1)-X(KA+2))*COS(X(KA)+X(KA+1)-X(KA+2))
          GA(KA)=GA(KA)+9.0D0*(X(KA)-X(KA+2))**2+D1S+D2S
          GA(KA+1)=2.0D0-D1S+D2S
          GA(KA+2)=-9.0D0*(X(KA)-X(KA+2))**2-D1S-D2S
        END IF
      ELSE
        EX=EXP(X(KA-1)-X(KA)-X(KA+1))
        W=X(KA-1)-X(KA+1)
        GA(KA-1)=-EX-W*EX
        GA(KA)=4.0D0+W*EX
        GA(KA+1)=EX+W*EX
      END IF
      RETURN
   70 H=2.0D0
      IF (KA.EQ.1) THEN
        GA(1)=2.0D0*((3.0D0-H*X(1))*X(1)-2.0D0*X(2)+1.0D0)*(3.0D0-2.0D0*
     &   H*X(1))
        GA(2)=-4.0D0*((3.0D0-H*X(1))*X(1)-2.0D0*X(2)+1.0D0)
      ELSE IF (KA.LE.N-1) THEN
        GA(KA-1)=-2.0D0*((3.0D0-H*X(KA))*X(KA)-X(KA-1)-2.0D0*X(KA+1)+
     &   1.0D0)
        GA(KA)=2.0D0*((3.0D0-H*X(KA))*X(KA)-X(KA-1)-2.0D0*X(KA+1)+1.0D0)
     &   *(3.0D0-2.0D0*H*X(KA))
        GA(KA+1)=-4.0D0*((3.0D0-H*X(KA))*X(KA)-X(KA-1)-2.0D0*X(KA+1)+
     &   1.0D0)
      ELSE
        GA(N-1)=-2.0D0*((3.0D0-H*X(N))*X(N)-X(N-1)+1.0D0)
        GA(N)=2.0D0*((3.0D0-H*X(N))*X(N)-X(N-1)+1.0D0)*(3.0D0-2.0D0*H*
     &   X(N))
      END IF
      RETURN
   80 IF (KA.LT.2) THEN
        GA(KA)=4.0D0
        GA(KA+1)=-8.0D0*X(KA+1)
      ELSE IF (KA.LT.N) THEN
        GA(KA-1)=-8.0D0*X(KA)
        GA(KA)=24.0D0*X(KA)**2-8.0D0*X(KA-1)+6.0D0
        GA(KA+1)=-8.0D0*X(KA+1)
      ELSE
        GA(KA-1)=-8.0D0*X(KA)
        GA(KA)=24.0D0*X(KA)**2-8.0D0*X(KA-1)+2.0D0
      END IF
      RETURN
   90 IF (KA.LT.2) THEN
        GA(KA)=4.0D0
        GA(KA+1)=-8.0D0*X(KA+1)+1.0D0
        GA(KA+2)=-2.0D0*X(KA+2)
      ELSE IF (KA.LT.3) THEN
        GA(KA-1)=-8.0D0*X(KA)
        GA(KA)=24.0D0*X(KA)**2-8.0D0*X(KA-1)+6.0D0
        GA(KA+1)=-8.0D0*X(KA+1)+1.0D0
        GA(KA+2)=-2.0D0*X(KA+2)
      ELSE IF (KA.LT.N-1) THEN
        GA(KA-2)=-1.0D0
        GA(KA-1)=-8.0D0*X(KA)+2.0D0*X(KA-1)
        GA(KA)=24.0D0*X(KA)**2-8.0D0*X(KA-1)+6.0D0
        GA(KA+1)=-8.0D0*X(KA+1)+1.0D0
        GA(KA+2)=-2.0D0*X(KA+2)
      ELSE IF (KA.LT.N) THEN
        GA(KA-2)=-1.0D0
        GA(KA-1)=-8.0D0*X(KA)+2.0D0*X(KA-1)
        GA(KA)=24.0D0*X(KA)**2-8.0D0*X(KA-1)+6.0D0
        GA(KA+1)=-8.0D0*X(KA+1)
      ELSE
        GA(KA-2)=-1.0D0
        GA(KA-1)=-8.0D0*X(KA)+2.0D0*X(KA-1)
        GA(KA)=24.0D0*X(KA)**2-8.0D0*X(KA-1)+2.0D0
      END IF
      RETURN
  100 IF (KA.LT.2) THEN
        GA(KA)=4.0D0
        GA(KA+1)=-8.0D0*X(KA+1)+1.0D0
        GA(KA+2)=-2.0D0*X(KA+2)+1.0D0
        GA(KA+3)=-2.0D0*X(KA+3)
      ELSE IF (KA.LT.3) THEN
        GA(KA-1)=-8.0D0*X(KA)+2.0D0*X(KA-1)
        GA(KA)=24.0D0*X(KA)**2-8.0D0*X(KA-1)+6.0D0
        GA(KA+1)=-8.0D0*X(KA+1)+1.0D0
        GA(KA+2)=-2.0D0*X(KA+2)+1.0D0
        GA(KA+3)=-2.0D0*X(KA+3)
      ELSE IF (KA.LT.4) THEN
        GA(KA-2)=2.0D0*X(KA-2)-1.0D0
        GA(KA-1)=-8.0D0*X(KA)+2.0D0*X(KA-1)
        GA(KA)=24.0D0*X(KA)**2-8.0D0*X(KA-1)+6.0D0
        GA(KA+1)=-8.0D0*X(KA+1)+1.0D0
        GA(KA+2)=-2.0D0*X(KA+2)+1.0D0
        GA(KA+3)=-2.0D0*X(KA+3)
      ELSE IF (KA.LT.N-2) THEN
        GA(KA-3)=-1.0D0
        GA(KA-2)=2.0D0*X(KA-2)-1.0D0
        GA(KA-1)=-8.0D0*X(KA)+2.0D0*X(KA-1)
        GA(KA)=24.0D0*X(KA)**2-8.0D0*X(KA-1)+6.0D0
        GA(KA+1)=-8.0D0*X(KA+1)+1.0D0
        GA(KA+2)=-2.0D0*X(KA+2)+1.0D0
        GA(KA+3)=-2.0D0*X(KA+3)
      ELSE IF (KA.LT.N-1) THEN
        GA(KA-3)=-1.0D0
        GA(KA-2)=2.0D0*X(KA-2)-1.0D0
        GA(KA-1)=-8.0D0*X(KA)+2.0D0*X(KA-1)
        GA(KA)=24.0D0*X(KA)**2-8.0D0*X(KA-1)+6.0D0
        GA(KA+1)=-8.0D0*X(KA+1)+1.0D0
        GA(KA+2)=-2.0D0*X(KA+2)+1.0D0
      ELSE IF (KA.LT.N) THEN
        GA(KA-3)=-1.0D0
        GA(KA-2)=2.0D0*X(KA-2)-1.0D0
        GA(KA-1)=-8.0D0*X(KA)+2.0D0*X(KA-1)
        GA(KA)=24.0D0*X(KA)**2-8.0D0*X(KA-1)+6.0D0
        GA(KA+1)=-8.0D0*X(KA+1)+1.0D0
      ELSE
        GA(KA-3)=-1.0D0
        GA(KA-2)=2.0D0*X(KA-2)-1.0D0
        GA(KA-1)=-8.0D0*X(KA)+2.0D0*X(KA-1)
        GA(KA)=24.0D0*X(KA)**2-8.0D0*X(KA-1)+2.0D0
      END IF
      RETURN
  110 IF (KA.EQ.1) THEN
        GA(N-4)=3.0D0
        GA(N-3)=-1.0D0
        GA(N-2)=-1.0D0
        GA(N-1)=0.50D0
        GA(N)=-1.0D0
        GA(KA)=-4.0D0*X(KA)+3.0D0
        GA(KA+1)=-2.0D0
      ELSE IF (KA.LE.N-1) THEN
        GA(KA-1)=0.0D0
        GA(KA)=0.0D0
        GA(KA+1)=0.0D0
        GA(N-4)=3.0D0
        GA(N-3)=-1.0D0
        GA(N-2)=-1.0D0
        GA(N-1)=0.50D0
        GA(N)=-1.0D0
        GA(KA-1)=GA(KA-1)-1.0D0
        GA(KA)=GA(KA)-4.0D0*X(KA)+3.0D0
        GA(KA+1)=GA(KA+1)-2.0D0
      ELSE
        GA(N-4)=3.0D0
        GA(N-3)=-1.0D0
        GA(N-2)=-1.0D0
        GA(N-1)=-0.5D0
        GA(N)=-4.0D0*X(N)+2.0D0
      END IF
      RETURN
  120 IF (MOD(KA,2).EQ.1) THEN
        GA(KA)=1.0D0
        GA(KA+1)=10.0D0*X(KA+1)-3.0D0*X(KA+1)**2-2.0D0
      ELSE
        GA(KA-1)=1.0D0
        GA(KA)=3.0D0*X(KA)**2+2.0D0*X(KA)-1.4D1
      END IF
      RETURN
  130 IF (MOD(KA,4).EQ.1) THEN
        GA(KA)=1.0D0
        GA(KA+1)=1.0D1
      ELSE IF (MOD(KA,4).EQ.2) THEN
        GA(KA+1)=2.23606797749979D0
        GA(KA+2)=-GA(KA+1)
      ELSE IF (MOD(KA,4).EQ.3) THEN
        GA(KA-1)=2.0D0*(X(KA-1)-2.0D0*X(KA))
        GA(KA)=-4.0D0*(X(KA-1)-2.0D0*X(KA))
      ELSE
        GA(KA-3)=3.16227766016838D0*2.0D0*(X(KA-3)-X(KA))
        GA(KA)=-3.16227766016838D0*2.0D0*(X(KA-3)-X(KA))
      END IF
      RETURN
  140 IF (MOD(KA,4).EQ.1) THEN
        GA(KA)=2.0D0*(EXP(X(KA))-X(KA+1))*EXP(X(KA))
        GA(KA+1)=-2.0D0*(EXP(X(KA))-X(KA+1))
      ELSE IF (MOD(KA,4).EQ.2) THEN
        GA(KA)=3.0D1*(X(KA)-X(KA+1))**2
        GA(KA+1)=-GA(KA)
      ELSE IF (MOD(KA,4).EQ.3) THEN
        GA(KA)=2.0D0*SIN(X(KA)-X(KA+1))/(COS(X(KA)-X(KA+1)))**3
        GA(KA+1)=-GA(KA)
      ELSE
        GA(KA)=1.0D0
      END IF
      RETURN
  150 IF (KA.LT.2) THEN
        GA(KA)=X(KA)-3.0D0
        GA(KA+1)=2.0D0
      ELSE IF (KA.LT.N) THEN
        GA(KA-1)=1.0D0
        GA(KA)=X(KA)-3.0D0
        GA(KA+1)=2.0D0
      ELSE
        GA(KA-1)=1.0D0
        GA(KA)=X(KA)-3.0D0
      END IF
      RETURN
  160 DO 170 J=MAX(1,KA-5),MIN(N,KA+1)
        GA(J)=1.0D0+2.0D0*X(J)
  170 CONTINUE
      GA(KA)=GA(KA)+2.0D0+15.0D0*X(KA)**2
      RETURN
  180 IF (MOD(KA,2).EQ.1) THEN
        GA(KA)=1.0D4*X(KA+1)
        GA(KA+1)=1.0D4*X(KA)
      ELSE
        GA(KA-1)=-EXP(-X(KA-1))
        GA(KA)=-EXP(-X(KA))
      END IF
      RETURN
  190 IF (MOD(KA,4).EQ.1) THEN
        GA(KA)=-2.0D2*(X(KA+1)-3.0D0*X(KA)**2)+1.0D0
        GA(KA+1)=-2.0D2*X(KA)
      ELSE IF (MOD(KA,4).EQ.2) THEN
        GA(KA-1)=-4.0D2*X(KA-1)
        GA(KA)=2.202D2
        GA(KA+2)=1.98D1
      ELSE IF (MOD(KA,4).EQ.3) THEN
        GA(KA)=-1.8D2*(X(KA+1)-3.0D0*X(KA)**2)+1.0D0
        GA(KA+1)=-1.8D2*X(KA)
      ELSE
        GA(KA-2)=1.98D1
        GA(KA-1)=-3.6D2*X(KA-1)
        GA(KA)=2.002D2
      END IF
      RETURN
  200 IF (KA.LT.2) THEN
        GA(KA+1)=EXP(COS(DBLE(KA)*(X(KA)+X(KA+1))))*DBLE(KA)*
     &   SIN(DBLE(KA)*(X(KA)+X(KA+1)))
        GA(KA)=GA(KA+1)+1.0D0
      ELSE IF (KA.LT.N) THEN
        GA(KA-1)=EXP(COS(DBLE(KA)*(X(KA-1)+X(KA)+X(KA+1))))*SIN(DBLE(KA)
     &   *(X(KA-1)+X(KA)+X(KA+1)))*DBLE(KA)
        GA(KA+1)=GA(KA-1)
        GA(KA)=1.0D0+GA(KA-1)
      ELSE
        GA(KA-1)=EXP(COS(DBLE(KA)*(X(KA-1)+X(KA))))*SIN(DBLE(KA)*(X(KA-
     &   1)+X(KA)))*DBLE(KA)
        GA(KA)=1.0D0+GA(KA-1)
      END IF
      RETURN
  210 U=1.0D0/DBLE(N+1)
      V=DBLE(KA)*U
      GA(KA)=2.0D0+1.5D0*U**2*(X(KA)+V+1.0D0)**2
      IF (KA.GT.1) GA(KA-1)=-1.0D0
      IF (KA.LT.N) GA(KA+1)=-1.0D0
      RETURN
  220 IF (KA.EQ.1) THEN
        GA(KA)=3.0D0*(X(KA+1)-4.0D0*X(KA))
        GA(KA+1)=3.0D0*X(KA)+0.5D0*X(KA+1)
      ELSE IF (KA.EQ.N) THEN
        GA(KA-1)=3.0D0*X(KA)-0.5D0*(2.0D1-X(KA-1))
        GA(KA)=3.0D0*(2.0D1-4.0D0*X(KA)+X(KA-1))
      ELSE
        GA(KA-1)=3.0D0*X(KA)-0.5D0*(X(KA+1)-X(KA-1))
        GA(KA)=3.0D0*(X(KA+1)-4.0D0*X(KA)+X(KA-1))
        GA(KA+1)=3.0D0*X(KA)+0.5D0*(X(KA+1)-X(KA-1))
      END IF
      RETURN
  230 H=1.0D0/DBLE(N+1)
      IF (KA.LT.2) THEN
        GA(KA)=2.0D0+PAR**2*H**2*COSH(PAR*X(KA))
        GA(KA+1)=-1.0D0
      ELSE IF (KA.LT.N) THEN
        GA(KA-1)=-1.0D0
        GA(KA)=2.0D0+PAR**2*H**2*COSH(PAR*X(KA))
        GA(KA+1)=-1.0D0
      ELSE
        GA(KA)=2.0D0+PAR**2*H**2*COSH(PAR*X(KA))
        GA(KA-1)=-1.0D0
      END IF
      RETURN
  240 GA(KA)=6.0D0
!      FA=6.0D 0*X(KA)
      A1=0.0D0
      A2=0.0D0
      A3=0.0D0
      GA1(1)=0.0D0
      GA1(2)=0.0D0
      GA2(1)=0.0D0
      GA2(2)=0.0D0
      IF (KA.GT.1) THEN
        GA(KA-1)=-4.0D0+PAR*X(KA)
        A1=A1-X(KA-1)
        A2=A2+X(KA-1)
        A3=A3+2.0D0*X(KA-1)
        GA1(1)=-1.0D0
        GA2(1)=1.0D0
      END IF
      IF (KA.GT.2) THEN
        GA(KA-2)=1.0D0-0.5D0*PAR*X(KA)
        A3=A3-X(KA-2)
      END IF
      IF (KA.LT.N-1) THEN
        GA(KA+2)=1.0D0+0.5D0*PAR*X(KA)
        A3=A3+X(KA+2)
      END IF
      IF (KA.LT.N) THEN
        GA(KA+1)=-4.0D0-PAR*X(KA)
        A1=A1+X(KA+1)
        A2=A2+X(KA+1)
        A3=A3-2.0D0*X(KA+1)
        GA1(2)=1.0D0
        GA2(2)=1.0D0
      END IF
      IF (KA.GE.N-1) THEN
        A3=A3+1.0D0
      END IF
      IF (KA.GE.N) THEN
        A1=A1+1.0D0
        A2=A2+1.0D0
        A3=A3-2.0D0
      END IF
      GA(KA)=GA(KA)+0.5D0*PAR*A3
      IF (KA.GT.1) GA(KA-1)=GA(KA-1)-0.5D0*PAR*(GA1(1)*A2+A1*GA2(1))
      IF (KA.LT.N) GA(KA+1)=GA(KA+1)-0.5D0*PAR*(GA1(2)*A2+A1*GA2(2))
      RETURN
  250 H=1.0D0/(M+1)
      IF (KA.LE.M) THEN
        J=KA+M
        GA(KA)=6.0D0
        GA(J)=0.0D0
        A1=0.0D0
        A2=0.0D0
        IF (KA.EQ.1) THEN
          A1=A1+1.0D0
        END IF
        IF (KA.GT.1) THEN
          A1=A1-X(J-1)
          A2=A2+2.0D0*X(KA-1)
          GA(KA-1)=-4.0D0+PAR*H*X(KA)
          GA(J-1)=-0.5D0*PAR*H**3*X(J)
        END IF
        IF (KA.GT.2) THEN
          A2=A2-X(KA-2)
          GA(KA-2)=1.0D0-0.5D0*PAR*H*X(KA)
        END IF
        IF (KA.LT.M-1) THEN
          A2=A2+X(KA+2)
          GA(KA+2)=1.0D0+0.5D0*PAR*H*X(KA)
        END IF
        IF (KA.LT.M) THEN
          A1=A1+X(J+1)
          A2=A2-2.0D0*X(KA+1)
          GA(KA+1)=-4.0D0-PAR*H*X(KA)
          GA(J+1)=0.5D0*PAR*H**3*X(J)
        END IF
        IF (KA.EQ.M) THEN
          A1=A1+1.0D0
        END IF
        GA(KA)=GA(KA)+0.5D0*PAR*H*A2
        GA(J)=GA(J)+0.5D0*PAR*H**3*A1
      ELSE
        J=KA-M
        GA(KA)=-2.0D0
        GA(J)=0.0D0
        A1=0.0D0
        A2=0.0D0
        IF (J.EQ.1) THEN
          A2=A2+1.0D0
        END IF
        IF (J.GT.1) THEN
          A1=A1-X(J-1)
          A2=A2-X(KA-1)
          GA(KA-1)=1.0D0-0.5D0*PAR*H*X(J)
          GA(J-1)=-0.5D0*PAR*H*X(KA)
        END IF
        IF (J.LT.M) THEN
          A1=A1+X(J+1)
          A2=A2+X(KA+1)
          GA(KA+1)=1.0D0+0.5D0*PAR*H*X(J)
          GA(J+1)=0.5D0*PAR*H*X(KA)
        END IF
        IF (J.EQ.M) THEN
          A2=A2+1.0D0
        END IF
        GA(KA)=GA(KA)+0.5D0*PAR*H*A1
        GA(J)=GA(J)+0.5D0*PAR*H*A2
      END IF
      RETURN
  260 GA(KA)=4.0D0-PAR*EXP(X(KA))
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      IF (J.GT.1) GA(KA-M)=-1.0D0
      IF (I.GT.1) GA(KA-1)=-1.0D0
      IF (I.LT.M) GA(KA+1)=-1.0D0
      IF (J.LT.M) GA(KA+M)=-1.0D0
      RETURN
  270 GA(KA)=4.0D0
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      GA(KA)=GA(KA)+3.0D0*PAR*X(KA)**2/(1.0D0+PAR*DBLE(I)**2+PAR*DBLE(J)
     &**2)
      IF (I.GT.1) GA(KA-1)=-1.0D0
      IF (I.LT.M) GA(KA+1)=-1.0D0
      IF (J.GT.1) GA(KA-M)=-1.0D0
      IF (J.LT.M) GA(KA+M)=-1.0D0
      RETURN
  280 GA(KA)=4.0D0-2.0D0*PI*PAR*COS(2.0D0*PI*X(KA))
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      A1=DBLE(I)/DBLE(M+1)
      A2=DBLE(J)/DBLE(M+1)
      IF (I.EQ.1) GA(KA+1)=-1.0D0-PI*DBLE(M+1)*PAR*COS(PI*X(KA+1)*
     &DBLE(M+1))
      IF (I.GT.1.AND.I.LT.M) THEN
        GA(KA-1)=-1.0D0+PI*DBLE(M+1)*PAR*COS(PI*(X(KA+1)-X(KA-1))*
     &   DBLE(M+1))
        GA(KA+1)=-1.0D0-PI*DBLE(M+1)*PAR*COS(PI*(X(KA+1)-X(KA-1))*
     &   DBLE(M+1))
      END IF
      IF (I.EQ.M) GA(KA-1)=-1.0D0+PI*DBLE(M+1)*PAR*COS(PI*X(KA-1)*
     &DBLE(M+1))
      IF (J.EQ.1) GA(KA+M)=-1.0D0-PI*DBLE(M+1)*PAR*COS(PI*X(KA+M)*
     &DBLE(M+1))
      IF (J.GT.1.AND.J.LT.M) THEN
        GA(KA-M)=-1.0D0+PI*DBLE(M+1)*PAR*COS(PI*(X(KA+M)-X(KA-M))*
     &   DBLE(M+1))
        GA(KA+M)=-1.0D0-PI*DBLE(M+1)*PAR*COS(PI*(X(KA+M)-X(KA-M))*
     &   DBLE(M+1))
      END IF
      IF (J.EQ.M) GA(KA-M)=-1.0D0+PI*DBLE(M+1)*PAR*COS(PI*X(KA-M)*
     &DBLE(M+1))
      RETURN
  290 GA(KA)=1.6D1*X(KA)
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      IF (I.EQ.1) THEN
        GA(KA)=GA(KA)-2.0D0*(X(KA+1)+1.0D0)-3.0D0*X(KA)*(X(KA+1)-1.0D0)*
     &   PAR
        GA(KA+1)=-2.0D0*X(KA)-(X(KA+1)-1.0D0)-1.5D0*X(KA)**2*PAR
      END IF
      IF (I.GT.1.AND.I.LT.M) THEN
        GA(KA)=GA(KA)-2.0D0*(X(KA+1)+X(KA-1))-3.0D0*X(KA)*(X(KA+1)-X(KA-
     &   1))*PAR
        GA(KA-1)=-2.0D0*X(KA)+(X(KA+1)-X(KA-1))+1.5D0*X(KA)**2*PAR
        GA(KA+1)=-2.0D0*X(KA)-(X(KA+1)-X(KA-1))-1.5D0*X(KA)**2*PAR
      END IF
      IF (I.EQ.M) THEN
        GA(KA)=GA(KA)-2.0D0*X(KA-1)+3.0D0*X(KA)*X(KA-1)*PAR
        GA(KA-1)=-2.0D0*X(KA)-X(KA-1)+1.5D0*X(KA)**2*PAR
      END IF
      IF (J.EQ.1) THEN
        GA(KA)=GA(KA)-2.0D0*(X(KA+M)+1.0D0)
        GA(KA+M)=-2.0D0*X(KA)-(X(KA+M)-1.0D0)
      END IF
      IF (J.GT.1.AND.J.LT.M) THEN
        GA(KA)=GA(KA)-2.0D0*(X(KA+M)+X(KA-M))
        GA(KA-M)=-2.0D0*X(KA)+(X(KA+M)-X(KA-M))
        GA(KA+M)=-2.0D0*X(KA)-(X(KA+M)-X(KA-M))
      END IF
      IF (J.EQ.M) THEN
        GA(KA)=GA(KA)-2.0D0*X(KA-M)
        GA(KA-M)=-2.0D0*X(KA)-X(KA-M)
      END IF
      RETURN
  300 GA(KA)=4.0D0
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      A1=PAR*DBLE(I)
      A2=PAR*DBLE(J)
      A3=0.0D0
      IF (I.GT.1) THEN
        GA(KA-1)=-1.0D0-2.0D1*PAR*X(KA)
        A3=A3-X(KA-1)
      END IF
      IF (I.LT.M) THEN
        GA(KA+1)=-1.0D0+2.0D1*PAR*X(KA)
        A3=A3+X(KA+1)
      END IF
      IF (J.GT.1) THEN
        GA(KA-M)=-1.0D0-2.0D1*PAR*X(KA)
        A3=A3-X(KA-M)
      END IF
      IF (J.LT.M) THEN
        GA(KA+M)=-1.0D0+2.0D1*PAR*X(KA)
        A3=A3+X(KA+M)
      END IF
      GA(KA)=GA(KA)+2.0D1*PAR*A3
      RETURN
  310 GA(KA)=2.0D1-PAR
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      IF (J.GT.2) THEN
        GA(KA-M-M)=1.0D0
      END IF
      IF (J.GT.1) THEN
        IF (I.GT.1) THEN
          GA(KA-M-1)=2.0D0
        END IF
        GA(KA-M)=-8.0D0
        IF (I.LT.M) THEN
          GA(KA-M+1)=2.0D0
        END IF
      END IF
      IF (I.GT.1) THEN
        IF (I.GT.2) THEN
          GA(KA-2)=1.0D0
        END IF
        GA(KA-1)=-8.0D0
      END IF
      IF (I.LT.M) THEN
        GA(KA+1)=-8.0D0
        IF (I.LT.M-1) THEN
          GA(KA+2)=1.0D0
        END IF
      END IF
      IF (J.LT.M) THEN
        IF (I.GT.1) THEN
          GA(KA+M-1)=2.0D0
        END IF
        GA(KA+M)=-8.0D0
        IF (I.LT.M) THEN
          GA(KA+M+1)=2.0D0
        END IF
      END IF
      IF (J.LT.M-1) THEN
        GA(KA+M+M)=1.0D0
      END IF
      RETURN
  320 H=0.5D0/DBLE(M+2)
      J=(KA-1)/M+1
      I=KA-(J-1)*M
      GA(KA)=2.0D1
      A1=0.0D0
      A2=0.0D0
      A3=0.0D0
      A4=0.0D0
      GA1(1)=0.0D0
      GA1(2)=0.0D0
      GA2(1)=0.0D0
      GA2(2)=0.0D0
      DO 330 K=1,6
        GA3(K)=0.0D0
        GA4(K)=0.0D0
  330 CONTINUE
      IF (J.GT.2) THEN
        GA(KA-M-M)=1.0D0
        GA4(1)=GA4(1)+1.0D0
        A4=A4+X(KA-M-M)
      END IF
      IF (J.GT.1) THEN
        IF (I.GT.1) THEN
          GA(KA-M-1)=2.0D0
          GA3(1)=GA3(1)+1.0D0
          GA4(2)=GA4(2)+1.0D0
          A3=A3+X(KA-M-1)
          A4=A4+X(KA-M-1)
        END IF
        GA(KA-M)=-8.0D0
        GA1(1)=GA1(1)-1.0D0
        A1=A1-X(KA-M)
        IF (I.LT.M) THEN
          GA(KA-M+1)=2.0D0
          GA3(2)=GA3(2)-1.0D0
          GA4(3)=GA4(3)+1.0D0
          A3=A3-X(KA-M+1)
          A4=A4+X(KA-M+1)
        END IF
      END IF
      IF (I.GT.1) THEN
        IF (I.GT.2) THEN
          GA(KA-2)=1.0D0
          GA3(3)=GA3(3)+1.0D0
          A3=A3+X(KA-2)
        END IF
        GA(KA-1)=-8.0D0
        GA2(1)=GA2(1)-1.0D0
        A2=A2-X(KA-1)
      END IF
      IF (I.LT.M) THEN
        GA(KA+1)=-8.0D0
        GA2(2)=GA2(2)+1.0D0
        A2=A2+X(KA+1)
        IF (I.LT.M-1) THEN
          GA(KA+2)=1.0D0
          GA3(4)=GA3(4)-1.0D0
          A3=A3-X(KA+2)
        END IF
      END IF
      IF (J.LT.M) THEN
        IF (I.GT.1) THEN
          GA(KA+M-1)=2.0D0
          GA3(5)=GA3(5)+1.0D0
          GA4(4)=GA4(4)-1.0D0
          A3=A3+X(KA+M-1)
          A4=A4-X(KA+M-1)
        END IF
        GA(KA+M)=-8.0D0
        GA1(2)=GA1(2)+1.0D0
        A1=A1+X(KA+M)
        IF (I.LT.M) THEN
          GA(KA+M+1)=2.0D0
          GA3(6)=GA3(6)-1.0D0
          GA4(5)=GA4(5)-1.0D0
          A3=A3-X(KA+M+1)
          A4=A4-X(KA+M+1)
        END IF
      END IF
      IF (J.LT.M-1) THEN
        GA(KA+M+M)=1.0D0
        GA4(6)=GA4(6)-1.0D0
        A4=A4-X(KA+M+M)
      END IF
      IF (J.EQ.M) THEN
        IF (I.GT.1) THEN
          A3=A3-H
          A4=A4+H
        END IF
        A1=A1-H
        IF (I.LT.M) THEN
          A3=A3+H
          A4=A4+H
        END IF
        A4=A4-H
      END IF
      IF (J.EQ.M-1) THEN
        A4=A4+H
      END IF
      IF (KA.GT.M+M) GA(KA-M-M)=GA(KA-M-M)+0.25D0*PAR*(-A2*GA4(1))
      IF (KA.GT.M+1) GA(KA-M-1)=GA(KA-M-1)+0.25D0*PAR*(+A1*GA3(1)-A2*
     &GA4(2))
      IF (KA.GT.M) GA(KA-M)=GA(KA-M)+0.25D0*PAR*(GA1(1)*A3)
      IF (KA.GT.M-1) GA(KA-M+1)=GA(KA-M+1)+0.25D0*PAR*(+A1*GA3(2)-A2*
     &GA4(3))
      IF (KA.GT.2) GA(KA-2)=GA(KA-2)+0.25D0*PAR*(+A1*GA3(3))
      IF (KA.GT.1) GA(KA-1)=GA(KA-1)+0.25D0*PAR*(-GA2(1)*A4)
      IF (KA.LE.N-1) GA(KA+1)=GA(KA+1)+0.25D0*PAR*(-GA2(2)*A4)
      IF (KA.LE.N-2) GA(KA+2)=GA(KA+2)+0.25D0*PAR*(+A1*GA3(4))
      IF (KA.LE.N-M+1) GA(KA+M-1)=GA(KA+M-1)+0.25D0*PAR*(+A1*GA3(5)-A2*
     &GA4(4))
      IF (KA.LE.N-M) GA(KA+M)=GA(KA+M)+0.25D0*PAR*(GA1(2)*A3)
      IF (KA.LE.N-M-1) GA(KA+M+1)=GA(KA+M+1)+0.25D0*PAR*(+A1*GA3(6)-A2*
     &GA4(5))
      IF (KA.LE.N-M-M) GA(KA+M+M)=GA(KA+M+M)+0.25D0*PAR*(-A2*GA4(6))
      RETURN
  340 H=1.0D0/DBLE(N+1)
      A1=DBLE(KA)*H
      A2=(A1-0.5D0)**2
      IF (A1.GE.0.5D0) THEN
        A3=1.0D6
      ELSE
        A3=-1.0D6
      END IF
      IF (KA.LT.2) THEN
        GA(KA)=2.0D0+H**2*(3.0D0*X(KA)**2*EXP(X(KA))+X(KA)**3*EXP(X(KA))
     &   )
        GA(KA+1)=-1.0D0+H**2*5.0D8*EXP(-1.0D4*A2)*SQRT(ABS(A1-0.5D0))
      ELSE IF (KA.LT.N) THEN
        GA(KA-1)=-1.0D0-H**2*5.0D8*EXP(-1.0D4*A2)*SQRT(ABS(A1-0.5D0))
        GA(KA)=2.0D0+H**2*(3.0D0*X(KA)**2*EXP(X(KA))+X(KA)**3*EXP(X(KA))
     &   )
        GA(KA+1)=-1.0D0+H**2*5.0D8*EXP(-1.0D4*A2)*SQRT(ABS(A1-0.5D0))
      ELSE
        GA(KA-1)=-1.0D0+H**2*5.0D8*EXP(-1.0D4*A2)*SQRT(ABS(A1-0.5D0))
        GA(KA)=2.0D0+H**2*(3.0D0*X(KA)**2*EXP(X(KA))+X(KA)**3*EXP(X(KA))
     &   )
      END IF
      RETURN
  350 H=1.0D0/DBLE(N+1)
      A1=DBLE(KA)*H
      A2=(A1-0.5D0)**2
      IF (KA.LT.2) THEN
        GA(KA+1)=-1.0D0
        GA(KA)=2.0D0+H**2*(3*X(KA)**2+2.0D-4*(2.0D-4*A2-1.0D0))
      ELSE IF (KA.LT.N) THEN
        GA(KA-1)=-1.0D0
        GA(KA)=2.0D0+H**2*(3*X(KA)**2+2.0D-4*(2.0D-4*A2-1.0D0))
        GA(KA+1)=-1.0D0
      ELSE
        GA(KA-1)=-1.0D0
        GA(KA)=2.0D0+H**2*(3*X(KA)**2+2.0D-4*(2.0D-4*A2-1.0D0))
      END IF
      RETURN
      END
! SUBROUTINE TYTIM1                MS DOS                     91/12/01
! PORTABILITY : MS DOS / MS FORTRAN V.5.0
! 91/12/01 SI : ORIGINAL VERSION
!
! PURPOSE :
!  GET TIME IN 100TH OF SEC.
!
      SUBROUTINE TYTIM1 (ITIME)
      INTEGER ITIME
      REAL TIME
      CALL CPU_TIME (TIME)
      ITIME=1.0D2*TIME
      END
! SUBROUTINE TYTIM2                ALL SYSTEMS                91/12/01
! PORTABILITY : ALL SYSTEMS
! 91/12/01 SI : ORIGINAL VERSION
!
! PURPOSE :
!  PRINT TIME ELAPSED.
!
      SUBROUTINE TYTIM2 (ITIME)
      INTEGER ITIME
      INTEGER IHR,IT,IMIN,ISEC
      CALL TYTIM1 (IT)
      IT=IT-ITIME
      IHR=IT/(60*60*100)
      IT=IT-IHR*60*60*100
      IMIN=IT/(60*100)
      IT=IT-IMIN*60*100
      ISEC=IT/100
      IT=IT-ISEC*100
      WRITE (6,10) IHR,IMIN,ISEC,IT
   10 FORMAT (' TIME=',I2,':',I2.2,':',I2.2,'.',I2.2)
      END

	subroutine bumpvec(N, H,ldh,sr1,si1,sr2,si2,v)
	implicit none
	integer ldh, n
	double precision h(ldh,2), sr1, si1, sr2, si2
c
c	returns a vector v that is parallel to the
c	first three entries of the first column of 
c
c		(H-sr1)*(H-sr2) - si1*si2*I.
c
c	This is most useful  is the cases that
c
c		si1 = si2 = 0  	(real shift pair)
c	and
c               sr1 =  sr2 	(complex conjugate shifts)
c		si1 = -si2   
c	
c
c	H:  first two columns of a 2x2 or 3x3 upper Hessenberg
c	N:  order of H, i.e., 2 or 3  (If N = 2, then H(3,:) and 
c		V(3) are unreferenced
c	sr1: 
c	si1:
c	sr2:
c	si2:
c	v:  length 3 output vector
c
	double precision zero
	parameter (zero = 0.0d0)
	double precision cn1, cn2, sn1, sn2, da, db, v(3)

	da = H(1,1) - sr2
	db = H(2,1) 
	call drotg(da,db,cn1,sn1)
	db = si2
	call drotg(da,db,cn2,sn2)

c	if the first column of H-s2*I  is zero, then 
c	the results are undefined.
c	if (da .eq. zero) then
c		v(1) = zero
c		v(2) = zero
c		if (n.GT.2) v(3) = zero
c	else
		v(1) = (H(1,1)-sr1)*cn1 + H(1,2)*sn1
		v(2) = H(2,1)*cn1 + (H(2,2)-sr1)*sn1
		if(N.GT.2) v(3) = H(3,2)*sn1

		v(1) = v(1)*cn2 - si1*sn2
		v(2) = v(2)*cn2
		if (N.GT.2) v(3) = v(3)*cn2
c	endif
	end
      SUBROUTINE CHASESINTRO(A, B, Q, Z,
     $     DESCA, DESCB, DESCQ, DESCZ,
     $     LDA, LDB, LDQ, LDZ,
     $     TMPW, TMPW2,
     $	 LDTMP,
     $     IBULGE, NBULGE, 
     $     IFIRST, ILAST, ILASTM, IFRSTM,
     $	 SAFMIN, ILO, NB, N,      
     $     ILQ, ILZ,
     $     DIST, LDIST, TALPHAR, TALPHAI, TBETA,
     $     WORK, JITER,ATOL,U, V, LDUV,
     $	 TT1, TT2, TT3, TT4, TT5)
      IMPLICIT NONE
*
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDQ, LDZ
      INTEGER            IBULGE, NBULGE, DIST
      INTEGER            IFIRST, ILAST, ILASTM, IFRSTM
      INTEGER            ILO, NB, N, JITER, LDIST
      LOGICAL            ILQ, ILZ
	INTEGER            UCNT, DCNT,LDTMP,LDUV
	DOUBLE PRECISION   TT1, TT2, TT3, TT4, TT5
	DOUBLE PRECISION   SAFMIN,ATOL
*     ..
*     .. Array Arguments ..
*     ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
      DOUBLE PRECISION   Q( LDQ, * ), Z( LDZ, * )
      INTEGER            DESCA(*), DESCB(*)
      INTEGER			   DESCQ(*), DESCZ(*)
      DOUBLE PRECISION   WORK( * )
      DOUBLE PRECISION   TMPW(LDTMP, *), TMPW2(LDTMP, *)
	DOUBLE PRECISION   U(LDUV,*), V(LDUV, *)
	DOUBLE PRECISION   TALPHAI(*), TALPHAR(*), TBETA(*)
*
*     Purpose
*     =======
*
*     CHASESHIFTS chases shifts numstep steps downwards the diagonal.
*     
*     KOLLA UPP NREF vs RT, J används ej i INDX som det är nu...
*     möjligen använda JJ i INDX
*     möjligen använda (NSBULGE+1 + DIST) som NNB i övre delen ??
*
*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )

      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
      INTEGER            ROWWISE, COLUMNWISE, ALL
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2, ALL=3)


*     ..
*     .. Local Arrays ..
*     ..      
      DOUBLE PRECISION   CS((NBULGE+1)*4)
	DOUBLE PRECISION   VV(4)
	
	LOGICAL            NDONE(NBULGE)
	LOGICAL            NDONE2(NBULGE)


*     ..
*     .. Local Scalars ..
*     ..
      INTEGER            J, JJ, IACOL, IAROW, IAROW2, IACOL2
	INTEGER			   OLDIAROW, OLDIACOL, IC, IR
      INTEGER            IB, RT, LT, RET, NNB, JI, SCNT
      INTEGER            ICTXT, MYROW, MYCOL, NPROW, NPCOL 
      DOUBLE PRECISION   D1, D2, D3, D4, D5, C, S, TEMP, TAU
	DOUBLE PRECISION   TAUU, TAUV, V1(3), V2(3)
	DOUBLE PRECISION   T1, T2

	INTEGER            I, INDX, NCOLS, JC, JR, DISPL,NSBULGE
	INTEGER			   LROW1, LCOL1, LROW2, LCOL2,ICHASE
	INTEGER            INDJ1, INDJ2, INDJ3, INDJ4, INDJ5
	INTEGER			   INDJ6, INDJ7, INDJ8, RSCNT, CSCNT
	INTEGER			   F,WSIZE, ICOL, IROW, NCHASE,ORIGJ
	INTEGER			   WSIZE2, DEBUG1, LR,LC
	LOGICAL			   SENTV, INCINDJ1, INCINDJ2, INCINDJ5
	EXTERNAL		   MPI_WTIME
	DOUBLE PRECISION   MPI_WTIME
*     ..
*     .. External Subroutines ..
*     ..
      EXTERNAL           BLACS_BARRIER, PDLACP3

*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT, MOD
*     ..
*     .. External Functions ..
*     ..
      INTEGER            INDXG2L, INDXG2P
      EXTERNAL           INDXG2L, INDXG2P

	IF (IFIRST.GE.ILAST) RETURN

	ICTXT = DESCA(CTXT_)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )  
	NSBULGE = 3
      J = IFIRST
      LT = 1
      RT = NBULGE*3*(NSBULGE + 2 + DIST+NB) + NBULGE*2
	DEBUG1 = 0
      IF (IFIRST+1.EQ.ILAST) THEN
         RETURN
      END IF
		
		
      DO 10 IB = 1, NBULGE
         NDONE(IB) = .TRUE.
         NDONE2(IB) = .TRUE.
 10   CONTINUE



		ICHASE = 0
		NCHASE = 1
		WSIZE = ((NSBULGE+1 + DIST) + 2) * NBULGE + 2

		WSIZE = MIN(WSIZE, ILAST-J+1)
		WSIZE = MIN(WSIZE, NB)

		WSIZE2 = MIN(NB, ILAST-J+1)
		WSIZE2 = WSIZE
*                WSIZE2 = MIN(WSIZE2, WSIZE+20)
          DO 500 IB = IBULGE, 1, -1

		  NCOLS = 3		  		  
		  JJ = J + (IB-1)*(NSBULGE+1 + DIST)
		  INDX = 3*(NSBULGE+DIST+2+MAX(NB,0))*(IB-1)

            
		              
            NNB = MIN((NSBULGE+1 + DIST) + NCOLS +1, ILAST-JJ+1)
            IF (NNB .LE. 0) GOTO 500
            
*           Last block ?
            IF ((JJ + NNB -1) .EQ. ILAST) THEN
               NCOLS = 2
               NDONE(IB) = .FALSE.
            ENDIF
	        
		  IR = (IB-1+ICHASE)*(NSBULGE+1 + DIST)+1
		  IC = (IB-1+ICHASE)*(NSBULGE+1 + DIST)+1

		  CALL QZ2Local( NNB, 
     $			TMPW(IR,IC),
     $			LDTMP,  
     $			TMPW2(IR,IC), 
     $			LDTMP,
     $			WORK(LT+INDX), WORK(RT+INDX), NCOLS, IB,
     $			WSIZE2-IR+1,
     $			ATOL, NB, (IBULGE-IB+1)*2,TALPHAI, TALPHAR, TBETA)
*		  .. Some column update is required before moving the next bulge
		  IF ((IB+ICHASE).GT.1) THEN            
			IF ((JJ-IFRSTM.GT.0).AND.(JJ-(NSBULGE+1+DIST).GT.0)) THEN

			  IC = (IB-1+ICHASE)*(NSBULGE+1+DIST)
			  JI = 1
			  CALL MYDLAREF( 'Cols', TMPW, 
     $			LDTMP, .FALSE., U, LDUV,
     $			.TRUE.,IR,IC+2,
     $			1,NNB-NCOLS-1,
     $            JI,IC, 0, 0,
     $            WORK(RT+INDX), D1, D2, D3, D4, D5 )	
			  IC = (IB-1+ICHASE)*(NSBULGE+1 + DIST)
			  JI = 1
     	  		  CALL MYDLAREF( 'Cols', TMPW2, 
     $			LDTMP, .FALSE., V, LDUV,
     $			.TRUE.,IR,IC+2,
     $			1,NNB-NCOLS-1,
     $            JI,IC, 0, 0,
     $            WORK(RT+INDX), D1, D2, D3, D4, D5 )	
		
              END IF
		  END IF





		  IF ((NCOLS.EQ.2).AND.(NDONE(1))) THEN

              JJ = ILAST - 1
               
			IR = WSIZE2 -2  
			IC = WSIZE2 -2

			TEMP = TMPW( IR+1, IC )
			CALL DLARTG( TEMP,TMPW(IR+2,IC),C,S,TMPW(IR+1,IC ) )
			TMPW(IR+2,IC) = ZERO
			CS(1+(IB-1)*4) = C
			CS(2+(IB-1)*4) = S	
			
*				Update rows of A				
			IR = WSIZE2 -1  
			IC = WSIZE2 -1
			CALL MYROT( 'R', TMPW, LDTMP, 
     $			IR, IC, 1,1,
     $			IC, IC+1, CS(1+(IB-1)*4))

*				Update rows of B
			CALL MYROT( 'R', TMPW2, LDTMP, 
     $			IR, IC, 1,1,
     $			IC, IC+1, CS(1+(IB-1)*4))

			IR = WSIZE2 -2  
			IC = WSIZE2 -2
				
							
			TEMP = TMPW2(IR+2,IC+2)
			CALL DLARTG(TEMP,TMPW2(IR+2,IC+1),C,S,TMPW2(IR+2,IC+2))
			TMPW2(IR+2,IC+1) = ZERO
			CS(3+(IB-1)*4) = C
			CS(4+(IB-1)*4) = S

*				Update cols of A				
			IR = WSIZE2 -2  
			IC = WSIZE2 -1
			CALL MYROT( 'C', TMPW, LDTMP, 
     $			IR, IC, 1,1,
     $			1, WSIZE2, CS(3+(IB-1)*4))

*				Update col of B				
			IR = WSIZE2 -2  
			IC = WSIZE2 -1
			CALL MYROT( 'C', TMPW2, LDTMP, 
     $			IR, IC, 1,1,
     $			1, WSIZE2-1, CS(3+(IB-1)*4))
	



		    

               
            END IF
 500      CONTINUE
         

          DO 600 IB = IBULGE, 1, -1

		  NCOLS = 3		  		  
		  JJ = J + (IB-1)*(NSBULGE+1+DIST)
		  INDX = 3*(NSBULGE+DIST+2+MAX(NB,0))*(IB-1)
		  
            IF (.NOT. NDONE2(IB)) THEN 
              GOTO 600
            END IF	
         
            NNB = MIN((NSBULGE+1 + DIST) + NCOLS +1, ILAST-JJ+1)
            IF (NNB .LE. 0) GOTO 600
            
*           Last block ?
            IF ((JJ + NNB -1) .EQ. ILAST) THEN
               NCOLS = 2
               NDONE2(IB) = .FALSE.
            ENDIF


							  
*		  Apply rowrotations to U and V 
		  DO 610 IR = 1, NNB-NCOLS-1
			VV(1) = ONE
			VV(2) = WORK(LT+INDX + (IR-1)*3)
			VV(3) = WORK(LT+INDX + (IR-1)*3+1)				
			TAU = WORK(LT+INDX + (IR-1)*3+2)

              CALL DLARFX( 'Right', LDUV, 3, VV, TAU,
     $            U(1,IR+ (IB-1+ICHASE)*(NSBULGE+1+DIST)+1),
     $			LDUV, WORK)
			VV(1) = ONE
			VV(2) = WORK(RT+INDX + (IR-1)*3)
			VV(3) = WORK(RT+INDX + (IR-1)*3+1) 				
			TAU = WORK(RT+INDX + (IR-1)*3+2)

              CALL DLARFX( 'Right', LDUV, 3, VV, TAU,
     $            V(1,IR+ (IB-1+ICHASE)*(NSBULGE+1+DIST)+1), 
     $			LDUV, WORK)
610		  CONTINUE
		  
*		  Clean up at border
	      IF ((NCOLS.EQ.2).AND.(NDONE2(1))) THEN
			IR = WSIZE2 -2  
			IC = WSIZE2 -1
			CALL MYROT( 'X', U, LDUV, 
     $				IR, IC, 1,1,
     $				1, LDUV, CS(1+(IB-1)*4))

			
			IR = WSIZE2 -2  
			IC = WSIZE2 -1
			CALL MYROT( 'C', V, LDUV, 
     $				IR, IC, 1,1,
     $				1, LDUV, CS(3+(IB-1)*4))
		  END IF

600		CONTINUE      
      RETURN
*
*     End of CHASESHIFTS
*
 800  FORMAT(40f7.2)
 850  FORMAT(40F7.2)
      END
      SUBROUTINE CHASESHIFTS2(A, B, Q, Z,
     $     DESCA, DESCB, DESCQ, DESCZ,
     $     LDA, LDB, LDQ, LDZ,
     $     TMPW, TMPW2,
     $	 LDTMP,
     $     IBULGE, NBULGE, 
     $     IFIRST, ILAST, ILASTM, IFRSTM,
     $     SAFMIN, ILO, NB, N,
     $     ILQ, ILZ, 
     $     DIST, LDIST, TALPHAR, TALPHAI, TBETA, 
     $     WORK, JITER,ATOL, U,V,LDUV,
     $	 TT1, TT2, TT3, TT4, TT5)
      IMPLICIT NONE

*
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDQ, LDZ
      INTEGER            IBULGE, NBULGE, DIST
      INTEGER            IFIRST, ILAST, ILASTM, IFRSTM
      INTEGER            ILO, NB, N, JITER, LDIST
      LOGICAL            ILQ, ILZ
	INTEGER            UCNT, DCNT, LDUV, LDTMP
	DOUBLE PRECISION   ATOL, SAFMIN
	DOUBLE PRECISION   TT1, TT2, TT3, TT4, TT5
*     ..
*     .. Array Arguments ..
*     ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
      DOUBLE PRECISION   Q( LDQ, * ), Z( LDZ, * )
      INTEGER            DESCA(*), DESCB(*)
      INTEGER			   DESCQ(*), DESCZ(*)
      DOUBLE PRECISION   WORK( * )
      DOUBLE PRECISION   U( * )
      DOUBLE PRECISION   V( * )

      DOUBLE PRECISION   TMPW(LDTMP, *), TMPW2(LDTMP, *)
	DOUBLE PRECISION   TALPHAI(*), TALPHAR(*), TBETA(*)


*
*     Purpose
*     =======
*
*     CHASESHIFTS chases shifts numstep steps downwards the diagonal.
*     
*     KOLLA UPP NREF vs RT, J används ej i INDX som det är nu...
*     möjligen använda JJ i INDX
*     möjligen använda DIST som NNB i övre delen ??
*
*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )

      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
      INTEGER            ROWWISE, COLUMNWISE, ALL
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2, ALL=3)


*     ..
*     .. Local Arrays ..
*     ..      
      DOUBLE PRECISION   CS((NBULGE+1)*4)
	LOGICAL            NDONE(NBULGE)
	LOGICAL            NDONE2(NBULGE)
	LOGICAL            NDONE21(NBULGE)
	LOGICAL            NDONE3(NBULGE)
	LOGICAL            NDONE4(NBULGE)
	LOGICAL            NDONE41(NBULGE)
	LOGICAL            NDONE5(NBULGE)


      INTEGER			   SENDARR(64,64)
	INTEGER			   RECVARR(64,64)

*     ..
*     .. Local Scalars ..
*     ..
      INTEGER            J, JJ, IACOL, IAROW, IAROW2, IACOL2
      INTEGER            IB, RT, LT, NREF, RET, NBL, NNB
      INTEGER            ICTXT, MYROW, MYCOL, NPROW, NPCOL 
      DOUBLE PRECISION   D1, D2, D3, D4, D5, C, S, TEMP
	DOUBLE PRECISION   AA(20,20), BB(20,20)
      INTEGER            DEBUG1, DEBUG2, DEBUG3, DEBUG4
	INTEGER            I, INDX, NCOLS, JC, JR, DISPL
	INTEGER			   LROW1, LCOL1, LROW2, LCOL2
	INTEGER			   NSBULGE, WSIZE, IR,IC,OLDIAROW, OLDIACOL
	INTEGER            INDJ1, INDJ2, INDJ3, INDJ4, INDJ5
	INTEGER			   INDJ6, INDJ7, INDJ8
*     ..
*     .. External Subroutines ..
*     ..
      EXTERNAL           BLACS_BARRIER, PDLACP3

*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT, MOD
*     ..
*     .. External Functions ..
*     ..
      INTEGER            INDXG2L, INDXG2P
      EXTERNAL           INDXG2L, INDXG2P

	ICTXT = DESCA(CTXT_)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )  

	NSBULGE = 3

      J = IFIRST
      LT = 1
      RT = NBULGE*3*(NSBULGE + 2 + DIST+NB) + NBULGE*2

      NREF = 1
C      NDONE = .TRUE.
      DEBUG1 = 0
      DEBUG2 = 0
      DEBUG3 = 0
	DEBUG4 = 0

	WSIZE = ((NSBULGE+1 + DIST) + 2) * IBULGE + 2
        WSIZE = NB - 5
	WSIZE = MIN(WSIZE, ILAST-J+1)
	WSIZE = MIN(WSIZE, NB)
      
      IF (IFIRST+1.GE.ILAST) THEN
         RETURN
      END IF

      DO 10 IB = 1, NBULGE
         NDONE(IB) = .TRUE.
         NDONE2(IB) = .TRUE.
	   NDONE4(IB) = .TRUE.
	   NDONE5(IB) = .TRUE.
 10   CONTINUE
      
         
 1000    CONTINUE


         DO 2000 IB = IBULGE, 1, -1
            
            IF (.NOT. NDONE(IB)) THEN                     									  
			GOTO 2000
            END IF	
*		  13, 9, 5, 1             
		  JJ = J + (IB-1)*(NSBULGE+1 + DIST)
		  INDX = 3*(NSBULGE+DIST+2+NB)*(IB-1)
		  
            
            IF(.NOT. NDONE(IB).OR. JJ.GT.(ILAST-2)) THEN
	         NDONE(IB) = .FALSE.
               GOTO 2000
            END IF
		  
		  NCOLS = 3
		              
            NNB = MIN(NSBULGE+DIST+NCOLS+2+MAX(NB-WSIZE-1,0),ILAST-JJ+1)
            IF (NNB .LE. 0) GOTO 2000
            
*           Last block ?
            IF ((JJ + NNB -1) .EQ. ILAST) THEN
               NCOLS = 2 
               NDONE(IB) = .FALSE.
            ENDIF
	         



            CALL PQZ2(NNB, NB, A, LDA, DESCA, JJ, 
     $           B, LDB, DESCB, WORK(LT+INDX), 
     $           WORK(RT+INDX), NCOLS,
     $           TALPHAI, TALPHAR, TBETA, SAFMIN, (IBULGE-IB+1)*2, 
     $		   TMPW, TMPW2,LDTMP,JITER,ATOL)
		 
		  


*		  .. Minor column update is required before moving the next bulge
		  IF (IB.GT.1) THEN            
			IF ((JJ-IFRSTM.GT.0).AND.(J.GT.0)) THEN
			IAROW = INDXG2P(JJ, NB, 0, 0, NPROW)
			IAROW2 = INDXG2P(J, NB, 0, 0, NPROW)

			IACOL = INDXG2P(JJ+1, NB, 0, 0, NPCOL)
			IACOL2 = INDXG2P(JJ+(NNB-(NCOLS+1)), NB, 0, 0, NPCOL)

			IF (IAROW.NE.IAROW2) THEN
				IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
					CALL Send(3*(NNB-(NCOLS+1)), 1, UP, 
     $					WORK(RT+INDX), 3*(NNB-(NCOLS+1)), 
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				END IF
				IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.IACOL)) THEN
					CALL Recv(3*(NNB-(NCOLS+1)), 1, DOWN, 
     $					WORK(RT+INDX), 3*(NNB-(NCOLS+1)),
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				END IF	
			END IF
			IF ((IACOL.NE.IACOL2).AND.(IAROW.NE.IAROW2)) THEN
				IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL2)) THEN
					CALL Send(3*(NNB-(NCOLS+1)), 1, UP, 
     $					WORK(RT+INDX), 3*(NNB-(NCOLS+1)),
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				END IF
				IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.IACOL2)) THEN
					CALL Recv(3*(NNB-(NCOLS+1)), 1, DOWN, 
     $					WORK(RT+INDX), 3*(NNB-(NCOLS+1)),
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				END IF	
			END IF

			CALL PDLAREF('Cols', A, LDA,
     $              N, NB, DESCA, 
     $              .TRUE.,
     $			  0, JJ+1, 
     $			  1, NNB-(NCOLS+1), 
     $              J, JJ-1,
     $              WORK(RT+INDX), D1, D2, D3, D4, D5,
     $              TMPW, RET,LDTMP)
			CALL PDLAREF('Cols', B, LDB,
     $              N, NB, DESCB, 
     $              .TRUE., 
     $			  0, JJ+1, 
     $			  1, NNB-(NCOLS+1), 
     $              J, JJ-1,
     $              WORK(RT+INDX), D1, D2, D3, D4, D5,
     $              TMPW, RET,LDTMP)
              END IF
		  END IF
		 


*		  .. Move the current bulge over the edge 
*		  .. only if the first bulge is not done
*		  .. If the first bulge is done then logic will take of this futher down
		  IF ((NCOLS.EQ.2).AND.(NDONE(1))) THEN
			
              JJ = ILAST - 1
               

			IAROW = INDXG2P(JJ, NB, 0, 0, NPROW)
			IACOL = INDXG2P(JJ, NB, 0, 0, NPCOL)
			IACOL2 = INDXG2P(JJ+1, NB, 0, 0, NPCOL)
			IAROW2 = INDXG2P(JJ+1, NB, 0, 0, NPROW)

			CALL PDLACP3(3, JJ-1, A, DESCA, TMPW,LDTMP,IAROW,IACOL, 0)
			CALL PDLACP3(2, JJ, B, DESCB, TMPW2, LDTMP, IAROW,IACOL,0)
			IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
				TEMP = TMPW( 2, 1 )
				CALL DLARTG( TEMP, TMPW( 3, 1 ), C, S, TMPW( 2, 1 ) )
				TMPW( 3, 1 ) = ZERO
			    CS(1+(IB-1)*4) = C
		 	    CS(2+(IB-1)*4) = S		
*				Update rows of A
                  TEMP = C*TMPW(2,2) + S*TMPW(3,2)
                  TMPW(3,2)  = -S*TMPW(2,2) + C*TMPW(3,2)
                  TMPW(2,2) = TEMP

                  TEMP = C*TMPW(2,3) + S*TMPW(3,3)
                  TMPW(3,3)  = -S*TMPW(2,3) + C*TMPW(3,3)
                  TMPW(2,3) = TEMP
*				Update rows of B
                  TEMP = C*TMPW2(1,1) + S*TMPW2(2,1)
                  TMPW2(2,1)  = -S*TMPW2(1,1) + C*TMPW2(2,1)
                  TMPW2(1,1) = TEMP

                  TEMP = C*TMPW2(1,2) + S*TMPW2(2,2)
                  TMPW2(2,2)  = -S*TMPW2(1,2) + C*TMPW2(2,2)
                  TMPW2(1,2) = TEMP
							
				TEMP = TMPW2( 2, 2 )
				CALL DLARTG( TEMP, TMPW2( 2, 1 ), C, S, TMPW2( 2, 2 ))
				TMPW2( 2, 1 ) = ZERO
			    CS(3+(IB-1)*4) = C
			    CS(4+(IB-1)*4) = S

*				Update cols of A
	            TEMP = C*TMPW(3,3) + S*TMPW(3,2)
		        TMPW(3,2) = -S*TMPW(3,3) + 
     $			     C*TMPW(3,2)
				TMPW(3,3) = TEMP

	            TEMP = C*TMPW(2,3) + S*TMPW(2,2)
		        TMPW(2,2) = -S*TMPW(2,3) + C*TMPW(2,2)
				TMPW(2,3) = TEMP

*				Update col of B
	            TEMP = C*TMPW2(1,2) + S*TMPW2(1,1)
		        TMPW2(1,1 ) = -S*TMPW2(1,2) + C*TMPW2(1,1)
				TMPW2(1,2) = TEMP
			END IF

			CALL PDLACP3(3, JJ-1, A,DESCA,TMPW,LDTMP, IAROW, IACOL,-1)
			CALL PDLACP3(2, JJ, B, DESCB, TMPW2, LDTMP,IAROW,IACOL,-1)
			

*             ..
*             .. Bcast Rotations from the left.and right
*             ..
		    IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
			   Call BSend((NBULGE+1)*4, 1,ROWWISE, CS, 
     $				(NBULGE+1)*4, ICTXT)
			   Call BSend((NBULGE+1)*4, 1, COLUMNWISE, CS, 
     $				(NBULGE+1)*4, ICTXT)
	        ELSE IF (MYROW.EQ.IAROW) THEN
			  Call BRecv((NBULGE+1)*4, 1, ROWWISE, CS, 
     $			(NBULGE+1)*4, 
     $			IAROW, IACOL, ICTXT)
		    ELSE IF (MYCOL.EQ.IACOL) THEN
			  Call BRecv((NBULGE+1)*4, 1, COLUMNWISE, CS, (NBULGE+1)*4, 
     $			IAROW, IACOL, ICTXT)
		    END IF
		    IF (IAROW.NE.IAROW2) THEN
		      DO 23 IC = 1, NPCOL 
		        IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.(IC-1))) THEN
					CALL Send((NBULGE+1)*4, 1, DOWN, CS, (NBULGE+1)*4,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
			     ELSE IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.(IC-1)))THEN
					CALL Recv((NBULGE+1)*4, 1, UP, CS, (NBULGE+1)*4,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
			     END IF
 23		      CONTINUE
	        END IF

		    IF (IACOL.NE.IACOL2) THEN
		      DO 24 IR = 1, NPROW 
		        IF ((MYCOL.EQ.IACOL).AND.(MYROW.EQ.(IR-1))) THEN
					CALL Send((NBULGE+1)*4, 1, RIGHT, CS, (NBULGE+1)*4,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
			    ELSE IF ((MYCOL.EQ.IACOL2).AND.(MYROW.EQ.(IR-1))) THEN
				  	CALL Recv((NBULGE+1)*4, 1, LEFT, CS, (NBULGE+1)*4,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
			    END IF 
 24		      CONTINUE
	        END IF

			CALL PDMYROT( 'C', A, DESCA, LDA, NB, JJ, JJ, 
     $			1, 1, J, ILAST-2,
     $			CS(3+(IB-1)*4), TMPW,LDTMP)  

			CALL PDMYROT( 'C', B, DESCB, LDB, NB, JJ, JJ, 
     $			1, 1, J, ILAST-2,
     $			CS(3+(IB-1)*4), TMPW,LDTMP) 



			
		  END IF
            
 2000    CONTINUE




*	   .. Accumulate transformations
	   IF ((NPROW.EQ.1).AND.(NPCOL.EQ.1)) THEN		
			GOTO 2090
	   END IF
	   

	   INDJ1 = 0
	   INDJ2 = 0

         DO 2040 IAROW=1, NPROW
		 DO 2040 IACOL=1, NPCOL
			SENDARR(IAROW, IACOL) = 0
			RECVARR(IAROW, IACOL) = 0
 2040	   CONTINUE

*	   GOTO 2055

         DO 2050 IB = 1, IBULGE

            IF (.NOT. NDONE2(IB)) THEN              
              GOTO 2050
            END IF	
            
		  JJ = J + (IB-1)*(NSBULGE+1 + DIST)
		  INDX = 3*(NSBULGE+DIST+2+NB)*(IB-1)
		              
            IF(.NOT. NDONE2(IB).OR.JJ.GT.(ILAST-2)) THEN
			 NDONE2(IB) = .FALSE.
               GOTO 2050
            END IF
		  
		  NCOLS = 3            
            NNB = MIN(NSBULGE+DIST+NCOLS+2+MAX(NB-WSIZE-1,0),ILAST-JJ+1)

            IF (NNB .LE. 0) GOTO 2050
            
*           Last block ?
            IF ((JJ + NNB -1) .EQ. ILAST) THEN
               NCOLS = 2
               NDONE2(IB) = .FALSE.
            ENDIF 
	                         
		  IF ((NNB-(NCOLS+1)).LE.0) THEN 
			GOTO 2050
		  END IF

            IAROW = INDXG2P(JJ+1, NB, 0, 0, NPROW)
            IACOL = INDXG2P(JJ+1, NB, 0, 0, NPCOL)

            IAROW2 = INDXG2P(JJ+(NNB-NCOLS)+1, NB, 0, 0, NPROW)
            IACOL2 = INDXG2P(JJ+(NNB-NCOLS)+1, NB, 0, 0, NPCOL)



		  UCNT = NB - MOD(JJ+1-1, NB)
		  UCNT = MIN(UCNT, NNB-(NCOLS+1))
		  DCNT = MOD(JJ+(NNB-NCOLS)+1, NB)
		  DCNT = MIN(DCNT, NNB-(NCOLS+1))

		  IF ((IAROW.EQ.IAROW2).AND.(IACOL.EQ.IACOL2)) THEN			 
			 SENDARR(IAROW+1,IACOL+1)=SENDARR(IAROW+1,IACOL+1)+1
			 IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN				
			    INDJ1 = INDJ1 + 1
				INDJ2 = INDJ2 + 1
				CALL DCOPY( 3*UCNT, WORK(LT+INDX), 1,
     $			   TMPW(1,INDJ1), 1)				
			    CALL DCOPY( 3*UCNT, WORK(RT+INDX), 1,
     $			   TMPW2(1,INDJ2), 1)
			 END IF	 
		  ELSE		 
			 SENDARR(IAROW+1,IACOL+1)=SENDARR(IAROW+1,IACOL+1)+1
			 SENDARR(IAROW2+1,IACOL2+1)=SENDARR(IAROW2+1,IACOL2+1)+1

			 IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN				
			    INDJ1 = INDJ1 + 1
				INDJ2 = INDJ2 + 1
			    CALL DCOPY( 3*UCNT, WORK(LT+INDX), 1,
     $			   TMPW(1,INDJ1), 1)				
			    CALL DCOPY( 3*UCNT, WORK(RT+INDX), 1,
     $			   TMPW2(1,INDJ2), 1)
			 END IF	 

			 IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.IACOL2).AND.
     $			((IACOL.NE.IACOL2).OR.(ILQ.AND.(IAROW.NE.IAROW2)))) THEN
			    INDJ1 = INDJ1 + 1 
				CALL DCOPY( 3*DCNT, 
     $				WORK(LT+INDX+3*((NNB-(NCOLS+1))-DCNT)), 1,
     $			    TMPW(1,INDJ1), 1)
			 END IF
			 IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.IACOL2).AND.
     $			(IAROW.NE.IAROW2)) THEN
				INDJ2 = INDJ2 + 1
			    CALL DCOPY( 3*DCNT,
     $				WORK(RT+INDX+3*((NNB-(NCOLS+1))-DCNT)), 1,
     $			    TMPW2(1,INDJ2), 1)
			 END IF	 
	      END IF
 2050    CONTINUE
 2055	   CONTINUE


*	   .. Distribute accumulated transformations

*	   GOTO 2065	   

*	   Previous start indx
	   INDJ3 = 1 
*	   NNB = MIN(DIST + 4, ILAST-J+1 )
	   DO 2060 JJ = J+1, ILAST, 1
*		  .. Distribute the accumulated transformations
		  IAROW = INDXG2P(JJ, NB, 0, 0, NPROW)
            IACOL = INDXG2P(JJ, NB, 0, 0, NPCOL)
		  
		  INDJ4 = SENDARR(IAROW+1, IACOL+1)
		  IF (INDJ4.EQ.0) THEN
			 GOTO 2060
	      END IF
		  SENDARR(IAROW+1,IACOL+1) = 0
		  
		  IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN			
			 IF (NPCOL.GT.1) THEN
			     Call BSend(LDTMP, INDJ4, 
     $				ROWWISE, TMPW(1, 1), LDTMP, ICTXT)
               END IF
	         IF (NPROW.GT.1) THEN
				Call BSend(LDTMP, INDJ4, 
     $				COLUMNWISE, TMPW2(1, 1), LDTMP, ICTXT)
			 END IF
			 IF (ILQ.AND.(NPROW.GT.1)) THEN
			     Call BSend(LDTMP, INDJ4, 
     $				COLUMNWISE, TMPW(1, 1), LDTMP, ICTXT)
	         END IF
		  ELSE IF ((MYROW.EQ.IAROW).AND.(NPCOL.GT.1)) THEN
			 Call BRecv(LDTMP, INDJ4, ROWWISE, 
     $				TMPW(1, INDJ3+LDTMP/2), LDTMP, IAROW, IACOL, ICTXT)
			 IF (RECVARR(IAROW+1, IACOL+1).EQ.0) THEN
				RECVARR(IAROW+1, IACOL+1) = INDJ3
	         END IF
		  ELSE IF ((MYCOL.EQ.IACOL).AND.(NPROW.GT.1)) THEN	
			 Call BRecv(LDTMP, INDJ4, COLUMNWISE, 
     $				TMPW2(1,INDJ3+LDTMP/2), LDTMP, IAROW, IACOL, ICTXT)
			 IF (RECVARR(IAROW+1, IACOL+1).EQ.0) THEN
				RECVARR(IAROW+1, IACOL+1) = INDJ3
	         END IF

			 IF (ILQ) THEN
				Call BRecv(LDTMP, INDJ4, COLUMNWISE, 
     $				TMPW(1, INDJ3+LDTMP/2), LDTMP, IAROW, IACOL, ICTXT)
	         END IF
		  END IF

		  INDJ3 = INDJ3 + INDJ4
 2060    CONTINUE

 2065    CONTINUE

*	   .. Restore accumulated transformations	  
	   INDJ1 = 0
	   INDJ2 = 0
	   INDJ5 = 0
	   INDJ6 = 0
	   INDJ7 = 0
	   INDJ8 = 0
*	   GOTO 2090	   
         DO 2070 IB = 1, IBULGE		
            

            IF (.NOT. NDONE4(IB)) THEN              
              GOTO 2070
            END IF	
            
		  JJ = J + (IB-1)*(NSBULGE+1 + DIST)
		  INDX = 3*(NSBULGE+DIST+2+NB)*(IB-1)
		  		            
            IF(.NOT. NDONE4(IB).OR.JJ.GT.(ILAST-2)) THEN
	         NDONE4(IB) = .FALSE.
               GOTO 2070
            END IF
		  
		  NCOLS = 3
            NNB = MIN(NSBULGE+DIST+NCOLS+2+MAX(NB-WSIZE-1,0),ILAST-JJ+1)

            IF (NNB .LE. 0) GOTO 2070            
*           .. Last block ?
            IF ((JJ + NNB -1) .EQ. ILAST) THEN
               NCOLS = 2
               NDONE4(IB) = .FALSE.
            ENDIF
	                         
		  IF ((NNB-(NCOLS+1)).LE.0) THEN 
               GOTO 2070
		  END IF

            IAROW = INDXG2P(JJ+1, NB, 0, 0, NPROW)
            IACOL = INDXG2P(JJ+1, NB, 0, 0, NPCOL)

            IAROW2 = INDXG2P(JJ+(NNB-NCOLS)+1, NB, 0, 0, NPROW)
            IACOL2 = INDXG2P(JJ+(NNB-NCOLS)+1, NB, 0, 0, NPCOL)

		  IF (IB.EQ.1) THEN
			OLDIAROW = IAROW
			OLDIACOL = IACOL
		  ELSE
			IF ((IAROW.NE.OLDIAROW).OR.(IACOL.NE.OLDIACOL)) THEN
				INDJ1 = INDJ2
				INDJ2 = 0
				OLDIAROW = IAROW
				OLDIACOL = IACOL
*				RECVARR(IAROW+1, IACOL+1) = RECVARR(OLDIAROW+1, OLDIACOL+1)
			END IF
		  END IF

		  INDJ3 = RECVARR(IAROW+1, IACOL+1)-1
		  INDJ4 = RECVARR(IAROW2+1, IACOL2+1)-1
	  
		  UCNT = NB - MOD(JJ+1-1, NB)
		  UCNT = MIN(UCNT, NNB-(NCOLS+1))
		  DCNT = MOD(JJ+(NNB-NCOLS)+1, NB)
		  DCNT = MIN(DCNT, NNB-(NCOLS+1))

		  IF ((IAROW.EQ.IAROW2).AND.(IACOL.EQ.IACOL2)) THEN
		     INDJ1 = INDJ1 + 1
			 IF ((MYROW.EQ.IAROW).AND.(MYCOL.NE.IACOL)) THEN     			 				
				CALL DCOPY( 3*UCNT, TMPW(1,INDJ3+INDJ1+LDTMP/2), 1,
     $			   WORK(LT+INDX), 1)		
     			 END IF

		     IF ((MYROW.NE.IAROW).AND.(MYCOL.EQ.IACOL)) THEN				
				CALL DCOPY( 3*UCNT, TMPW2(1,INDJ3+INDJ1+LDTMP/2), 1,
     $			   WORK(RT+INDX), 1)
			 END IF	 

			 IF ((MYROW.NE.IAROW).AND.(MYCOL.EQ.IACOL).AND.ILQ) THEN
				CALL DCOPY( 3*UCNT, TMPW(1,INDJ3+INDJ1+LDTMP/2), 1 ,
     $			   WORK(LT+INDX), 1)		
			 END IF




		  ELSE
			 INDJ1 = INDJ1 + 1
			 INDJ2 = INDJ2 + 1
			 IF ((MYROW.EQ.IAROW).AND.(MYCOL.NE.IACOL)) THEN
				CALL DCOPY( 3*UCNT, TMPW(1,INDJ3+INDJ1+LDTMP/2), 1 ,
     $			   WORK(LT+INDX), 1)		
			 END IF
			 			 
			 IF ((MYROW.NE.IAROW).AND.(MYCOL.EQ.IACOL)) THEN				
				CALL DCOPY( 3*UCNT, TMPW2(1,INDJ3+INDJ1+LDTMP/2), 1 ,
     $			   WORK(RT+INDX), 1)
			 END IF

			 IF ((MYROW.NE.IAROW).AND.(MYCOL.EQ.IACOL).AND.ILQ) THEN
				CALL DCOPY( 3*UCNT, TMPW(1,INDJ3+INDJ1+LDTMP/2), 1 ,
     $			   WORK(LT+INDX), 1)		
			 END IF

		     IF ((MYROW.EQ.IAROW2).AND.(MYCOL.NE.IACOL2).AND.
     $			((IAROW.NE.IAROW2).OR.(NPROW.EQ.1))) THEN			    
				CALL DCOPY( 3*DCNT, TMPW(1,INDJ4+INDJ2+LDTMP/2), 1 ,
     $				WORK(LT+INDX+3*((NNB-(NCOLS+1))-DCNT)), 1)
			 END IF	 
	 			 								 		
			 IF ((MYROW.NE.IAROW2).AND.(MYCOL.EQ.IACOL2).AND.
     $			((IACOL.NE.IACOL2).OR.(NPCOL.EQ.1))) THEN
				CALL DCOPY( 3*DCNT, TMPW2(1,INDJ4+INDJ2+LDTMP/2), 1,
     $				WORK(RT+INDX+3*((NNB-(NCOLS+1))-DCNT)), 1)
			 END IF	 

			 IF ((MYROW.NE.IAROW2).AND.(MYCOL.EQ.IACOL2).AND.
     $			((IACOL.NE.IACOL2).OR.(NPCOL.EQ.1)).AND.ILQ) THEN
				CALL DCOPY( 3*DCNT, TMPW(1,INDJ4+INDJ2+LDTMP/2), 1 ,
     $				WORK(LT+INDX+3*((NNB-(NCOLS+1))-DCNT)), 1)
			 END IF

	      END IF 
 2070    CONTINUE

*	   .. Apply row and column transformations

 2090	   CONTINUE
         DO 2100 IB = IBULGE, 1, -1

		  		  
            IF (.NOT. NDONE5(IB)) THEN 
              GOTO 2100
            END IF	
            
		  JJ = J + (IB-1)*(NSBULGE+1 + DIST)

		  INDX = 3*(NSBULGE+DIST+2+NB)*(IB-1)		  
            
            IF(.NOT. NDONE5(IB) .OR. JJ.GT.(ILAST-2)) THEN
			 NDONE5(IB) = .FALSE.
               GOTO 2100
            END IF	

		  NCOLS = 3            
            NNB = MIN(NSBULGE+DIST+NCOLS+2+MAX(NB-WSIZE-1,0),ILAST-JJ+1)

            IF (NNB .LE. 0) GOTO 2100
            
*           Last block ?
            IF ((JJ + NNB -1) .EQ. ILAST) THEN
               NCOLS = 2
               NDONE5(IB) = .FALSE.
            ENDIF



            CALL PDLAREF('Rows', A, LDA,
     $           N, NB, DESCA, 
     $          .TRUE., JJ+1, JJ+1, 1, NNB-(NCOLS+1),
     $           JJ+NNB, ILASTM,
     $           WORK(LT+INDX), D1, D2, D3, D4, D5,
     $           TMPW, RET,LDTMP)
            CALL PDLAREF('Rows', B, LDB,
     $           N, NB, DESCB, 
     $           .TRUE., JJ+1, JJ+1, 1, NNB-(NCOLS+1),
     $           JJ+NNB, ILASTM,
     $           WORK(LT+INDX), D1, D2, D3, D4, D5,
     $           TMPW, RET,LDTMP)

  		  IF (IB.GT.1) THEN
			IF (JJ-IFRSTM .GT. 0) THEN
			CALL PDLAREF('Cols', A, LDA,
     $              N, NB, DESCA, 
     $              .TRUE., IFRSTM, JJ+1, 1, NNB-(NCOLS+1), 
     $              IFRSTM, JJ-(NSBULGE+1+DIST)*(IB-1)-1,
     $              WORK(RT+INDX), D1, D2, D3, D4, D5,
     $              TMPW, RET,LDTMP)
			CALL PDLAREF('Cols', B, LDB,
     $              N, NB, DESCB, 
     $              .TRUE., IFRSTM, JJ+1, 1, NNB-(NCOLS+1), 
     $              IFRSTM, JJ-(NSBULGE+1+DIST)*(IB-1)-1,
     $              WORK(RT+INDX), D1, D2, D3, D4, D5,
     $              TMPW, RET,LDTMP)
			END IF
		  ELSE

			CALL PDLAREF('Cols', A, LDA,
     $              N, NB, DESCA, 
     $              .TRUE., IFRSTM, JJ+1, 1, NNB-(NCOLS+1), 
     $              IFRSTM, JJ-1,
     $              WORK(RT+INDX), D1, D2, D3, D4, D5,
     $              TMPW, RET,LDTMP)
			CALL PDLAREF('Cols', B, LDB,
     $              N, NB, DESCB, 
     $              .TRUE., IFRSTM, JJ+1, 1, NNB-(NCOLS+1), 
     $              IFRSTM, JJ-1,
     $              WORK(RT+INDX), D1, D2, D3, D4, D5,
     $              TMPW, RET,LDTMP)


		  END IF



	
            IF (ILQ) THEN
               CALL PDLAREF('Cols', Q, LDQ,
     $              N, NB, DESCQ, 
     $              .TRUE., 1, JJ+1, 1, NNB-(NCOLS+1), 
     $              1, N, WORK(LT+INDX), 
     $              D1, D2, D3, D4, D5,
     $              TMPW, RET,LDTMP)
            ENDIF

            IF (ILZ) THEN
               CALL PDLAREF('Cols', Z, LDZ,
     $              N, NB, DESCZ, 
     $              .TRUE., 1, JJ+1, 1, NNB-(NCOLS+1), 
     $              1, N, WORK(RT+INDX), 
     $              D1, D2, D3, D4, D5,
     $              TMPW, RET,LDTMP)
            ENDIF
 
		  IF ((NCOLS.EQ.2).AND.(NDONE5(1))) THEN
			
              JJ = ILAST - 1
             
			CALL PDMYROT( 'R', A, DESCA, LDA, NB, JJ, JJ, 
     $			1, 1, JJ+2, ILASTM,
     $	        CS(1+(IB-1)*4), TMPW,LDTMP)  


			CALL PDMYROT( 'R', B, DESCB, LDB, NB, JJ, JJ, 
     $			1, 1, JJ+2, ILASTM,
     $			CS(1+(IB-1)*4), TMPW,LDTMP)

			IF( ILQ ) THEN
				CALL PDMYROT( 'X', Q, DESCQ, LDQ, NB, JJ, JJ, 
     $				1, 1, 1, N,
     $				CS(1+(IB-1)*4), TMPW,LDTMP)
     			END IF

			CALL PDMYROT( 'C', A, DESCA, LDA, NB, JJ, JJ, 
     $			1, 1, IFRSTM, J - 1,
     $			CS(3+(IB-1)*4), TMPW,LDTMP)  			
			CALL PDMYROT( 'C', B, DESCB, LDB, NB, JJ, JJ, 
     $			1, 1, IFRSTM, J - 1,
     $			CS(3+(IB-1)*4), TMPW,LDTMP) 




			IF( ILZ ) THEN
				CALL PDMYROT( 'C', Z, DESCZ, LDZ, NB, JJ, JJ, 
     $			1, 1, 1, N,
     $			CS(3+(IB-1)*4), TMPW,LDTMP)
			END IF
			
		  END IF
 2100    CONTINUE







		
         
	    J = J + (NSBULGE+1+DIST+MAX(NB-WSIZE-1,0)) 



		IF (LDIST.GT.0.AND.J.LE.(ILAST)) THEN
			RETURN
		END IF

*	    .. Start with next blockmovement at j+(NSBULGE+1 + DIST) if the last bulge is not done

          IF ((J.LE.ILAST))  GOTO 1000

*	    .. Start with next blockmovement at j+dist if the last bulge is not done
*          IF (NDONE5(1))  GOTO 1000

*         J = J + NB - 4 
*         IF (NDONE(1))  GOTO 1000

*        .. Last elements: Use Givens rotations
*        ..
*        ..Rotations from the left
*        ..

         JJ = ILAST - 1

         IAROW = INDXG2P(JJ, NB, 0, 0, NPROW)
         IACOL = INDXG2P(JJ-1, NB, 0, 0, NPCOL)
	   IAROW2 = INDXG2P(JJ+1, NB, 0, 0, NPROW)
	   IF (IAROW2.EQ.IAROW) THEN
		  IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
			 LROW1 = INDXG2L(JJ, NB, 0, 0, NPROW)
			 LROW2 = INDXG2L(JJ+1, NB, 0, 0, NPROW)
			 LCOL1 = INDXG2L(JJ-1, NB, 0, 0, NPCOL)
			 TMPW(2,1) = A(LROW1, LCOL1)
			 TMPW(3,1) = A(LROW2, LCOL1)		
		  END IF
	   ELSE
		  CALL PDLACP3(3, JJ-1, A, DESCA, TMPW, NB, IAROW, IACOL, 0)
	   END IF

	   IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
			TEMP = TMPW( 2, 1 )
			CALL DLARTG( TEMP, TMPW( 3, 1 ), C, S, TMPW( 2, 1 ) )
			TMPW( 3, 1 ) = ZERO
	        CS(1) = C
              CS(2) = S	
	   END IF

	   IF (IAROW2.EQ.IAROW) THEN
		  IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
			 A(LROW1, LCOL1) = TMPW(2,1)
			 A(LROW2, LCOL1) = TMPW(3,1)
		  END IF
	   ELSE
		  CALL PDLACP3(3, JJ-1, A, DESCA, TMPW, NB, IAROW, IACOL, -1)
	   END IF


	   IF (IAROW.NE.IAROW2) THEN
		  IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
				CALL Send(2, 1, DOWN, 
     $				CS, 2, MYROW, MYCOL, NPCOL, NPROW, ICTXT)
		  END IF
		  IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.IACOL)) THEN
				CALL Recv(2, 1, UP, 
     $				CS, 2, MYROW, MYCOL, NPCOL, NPROW, ICTXT)
		  END IF
	   END IF

	   IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
			Call BSend(2, 1, ROWWISE, CS, 2, ICTXT)
	   ELSE IF (MYROW.EQ.IAROW) THEN
		    Call BRecv(2, 1, ROWWISE, CS, 2, 
     $			IAROW, IACOL, ICTXT)
	   END IF

	   IF (IAROW.NE.IAROW2) THEN
		  IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.IACOL)) THEN
			Call BSend(2, 1, ROWWISE, CS, 2, ICTXT)
		  ELSE IF (MYROW.EQ.IAROW2) THEN
		    Call BRecv(2, 1, ROWWISE, CS, 2, 
     $			IAROW2, IACOL, ICTXT)
	      END IF
	   END IF

         
*     
         IAROW = INDXG2P(JJ, NB, 0, 0, NPROW)
         IACOL = INDXG2P(JJ, NB, 0, 0, NPCOL)
	   IACOL2 = INDXG2P(JJ+1, NB, 0, 0, NPCOL)
	   
	   IF (ILQ) THEN
		  IF ((MYCOL.EQ.IACOL).AND.(MYROW.EQ.IAROW)) THEN
			  Call BSend(2, 1, COLUMNWISE, CS, 2, ICTXT)
		  ELSE IF (MYCOL.EQ.IACOL) THEN
		    Call BRecv(2, 1, COLUMNWISE, CS, 2, 
     $			IAROW, IACOL, ICTXT)
		  END IF
		  IF (IACOL.NE.IACOL2) THEN
			 IF ((MYCOL.EQ.IACOL2).AND.(MYROW.EQ.IAROW)) THEN
				Call BSend(2, 1, COLUMNWISE, CS, 2, ICTXT)
			 ELSE IF (MYCOL.EQ.IACOL2) THEN
				Call BRecv(2, 1, COLUMNWISE, CS, 2, 
     $				IAROW, IACOL2, ICTXT)
			 END IF
		  END IF

	   END IF
	   

         CALL PDMYROT( 'R', A, DESCA, LDA, NB, JJ, JJ, 
     $        1, 1, JJ, ILASTM,
     $        CS, TMPW,LDTMP)  


         CALL PDMYROT( 'R', B, DESCB, LDB, NB, JJ, JJ, 
     $        1, 1, JJ, ILASTM,
     $        CS, TMPW,LDTMP)

         IF( ILQ ) THEN
            CALL PDMYROT( 'X', Q, DESCQ, LDQ, NB, JJ, JJ, 
     $           1, 1, 1, N,
     $           CS, TMPW,LDTMP)
         END IF
*     0706099746



*        ..
*        .. Rotations from the right.
*        ..
         IAROW = INDXG2P(JJ+1, NB, 0, 0, NPROW)
         IACOL = INDXG2P(JJ, NB, 0, 0, NPCOL)
	   IACOL2 = INDXG2P(JJ+1, NB, 0, 0, NPCOL)
         IF (IACOL.EQ.IACOL2) THEN
		 IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN	   
			 LROW1 = INDXG2L(JJ+1, NB, 0, 0, NPROW)
			 LCOL1 = INDXG2L(JJ, NB, 0, 0, NPCOL)
			 LCOL2 = INDXG2L(JJ+1, NB, 0, 0, NPCOL)			 			
			 TMPW(2,1) = B(LROW1, LCOL1)
			 TMPW(2,2) = B(LROW1, LCOL2)		 
		 END IF	
	   ELSE
		  CALL PDLACP3(2, JJ, B, DESCB, TMPW, NB, IAROW, IACOL, 0)
	   END IF

	   IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN	   
			TEMP = TMPW( 2, 2 )
			CALL DLARTG( TEMP, TMPW( 2, 1 ), C, S, TMPW( 2, 2 ))
			TMPW( 2, 1 ) = ZERO
			CS(1) = C
			CS(2) = S
	   END IF

         IF (IACOL.EQ.IACOL2) THEN
		 IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN	   		
			 B(LROW1, LCOL1) = TMPW(2,1)
			 B(LROW1, LCOL2) = TMPW(2,2)
		 END IF	
	   ELSE
		  CALL PDLACP3(2, JJ, B, DESCB, TMPW, NB, IAROW, IACOL, -1)
	   END IF

	   IF (IACOL.NE.IACOL2) THEN
		  IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
				CALL Send(2, 1, RIGHT, 
     $				CS, 2, MYROW, MYCOL, NPCOL, NPROW, ICTXT)
		  END IF
		  IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL2)) THEN
				CALL Recv(2, 1, LEFT, 
     $				CS, 2, MYROW, MYCOL, NPCOL, NPROW, ICTXT)
		  END IF

	   END IF


	   IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
			Call BSend(2, 1, COLUMNWISE, CS, 2, ICTXT)
	   ELSE IF (MYCOL.EQ.IACOL) THEN
		    Call BRecv(2, 1, COLUMNWISE, CS, 2, 
     $			IAROW, IACOL, ICTXT)
	   END IF

	   IF (IACOL.NE.IACOL2) THEN
		  IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL2)) THEN
			Call BSend(2, 1, COLUMNWISE, CS, 2, ICTXT)
		  ELSE IF (MYCOL.EQ.IACOL2) THEN
		    Call BRecv(2, 1, COLUMNWISE, CS, 2, 
     $			IAROW, IACOL2, ICTXT)
	      END IF
	   END IF

         



         CALL PDMYROT( 'C', A, DESCA, LDA, NB, JJ, JJ, 
     $        1, 1, IFRSTM, ILAST,
     $        CS, TMPW,LDTMP)  


         CALL PDMYROT( 'C', B, DESCB, LDB, NB, JJ, JJ, 
     $        1, 1, IFRSTM, ILAST-1,
     $        CS, TMPW,LDTMP) 

         IF( ILZ ) THEN
            CALL PDMYROT( 'C', Z, DESCZ, LDZ, NB, JJ, JJ, 
     $           1, 1, 1, N,
     $           CS, TMPW,LDTMP)
         END IF



 3000 CONTINUE

      RETURN
*
*     End of CHASESHIFTS
*
 800  FORMAT(20f7.2)
      END
      SUBROUTINE CHASE_DISTUV_EXCH(
     $     IBULGE, 
     $     IFIRST, ILAST, 
     $	 NB, N,      
     $     DIST, 
     $     U, V, LDUV,
     $	 MOVEDIST,ILQ,DESCA,J)
      IMPLICIT NONE

*
*
*     .. Scalar Arguments ..
      INTEGER            IBULGE, DIST
      INTEGER            IFIRST, ILAST
      INTEGER            NB, N, LDIST,MOVEDIST
	INTEGER            LDTMP,LDUV
      LOGICAL            ILQ
*     ..
*     .. Array Arguments ..
*     ..
      INTEGER            DESCA(*)
	DOUBLE PRECISION   U(LDUV,*), V(LDUV, *)

      INTEGER            ICTXT, MYROW, MYCOL, NPROW, NPCOL 

      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
      INTEGER            ROWWISE, COLUMNWISE, ALL
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2, ALL=3)

	
	INTEGER			   NSBULGE, WSIZE, WSIZE2
      INTEGER            J, JJ, IACOL, IAROW, IAROW2, IACOL2
	INTEGER			   IC, IR

*     ..
*     .. External Functions ..
*     ..
      INTEGER            INDXG2L, INDXG2P
      EXTERNAL           INDXG2L, INDXG2P



	NSBULGE = 3


	IF (IFIRST.GE.ILAST) RETURN
	IF (J.GE.ILAST) RETURN

	ICTXT = DESCA(CTXT_)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )  
      IF (IFIRST+1.EQ.ILAST) THEN
         RETURN
      END IF
      IF (J+1.EQ.ILAST) THEN
         RETURN
      END IF

	    IAROW = INDXG2P(J, NB, 0, 0, NPROW)
	    IACOL = INDXG2P(J, NB, 0, 0, NPCOL)
		WSIZE = ((NSBULGE+1 + DIST) + 2) * IBULGE + 2
		WSIZE = MIN(WSIZE, ILAST-J+1)
		WSIZE = MIN(WSIZE, NB)
		WSIZE2 = MIN(NB, ILAST-J+1)
		WSIZE2 = MIN(WSIZE2, WSIZE+MOVEDIST)
	    IACOL2 = INDXG2P(J+WSIZE2-1, NB, 0, 0, NPCOL)
		IAROW2 = INDXG2P(J+WSIZE2-1, NB, 0, 0, NPROW)



		IF (IAROW.NE.IAROW2) THEN
		   DO 3501 IC = 1, NPCOL 
			 IF ((IC-1).EQ.IACOL) GOTO 3501
		     IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.(IC-1))) THEN
				CALL Send(WSIZE2,WSIZE2, DOWN, V, LDUV,
     $				MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				IF (ILQ) THEN
					CALL Send(WSIZE2,WSIZE2, DOWN, U, LDUV,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				END IF
			 ELSE IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.(IC-1))) THEN
				CALL Recv(WSIZE2,WSIZE2, UP, V, LDUV,
     $				MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				IF (ILQ) THEN
					CALL Recv(WSIZE2,WSIZE2, UP, U, LDUV,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				END IF
			 END IF
 3501		   CONTINUE
	    END IF

		IF (IACOL.NE.IACOL2) THEN
		   DO 3601 IR = 1, NPROW 
			 IF ((IR-1).EQ.IAROW) GOTO 3601
		     IF ((MYCOL.EQ.IACOL).AND.(MYROW.EQ.(IR-1))) THEN
				CALL Send(WSIZE2,WSIZE2, RIGHT, U, LDUV,
     $				MYROW, MYCOL, NPCOL, NPROW, ICTXT)
			 ELSE IF ((MYCOL.EQ.IACOL2).AND.(MYROW.EQ.(IR-1))) THEN
				CALL Recv(WSIZE2,WSIZE2, LEFT, U, LDUV,
     $				MYROW, MYCOL, NPCOL, NPROW, ICTXT)

			 END IF
 3601		   CONTINUE
	    END IF 



	END



      SUBROUTINE CHASE_DISTUV_ROW(
     $     IBULGE, 
     $     IFIRST, ILAST, 
     $	 NB, N,      
     $     DIST, 
     $     U, V, LDUV,
     $	 MOVEDIST,ILQ,DESCA,J)
      IMPLICIT NONE

*
*
*     .. Scalar Arguments ..
      INTEGER            IBULGE, DIST
      INTEGER            IFIRST, ILAST
      INTEGER            NB, N, LDIST,MOVEDIST
	INTEGER            LDTMP,LDUV
      LOGICAL            ILQ
*     ..
*     .. Array Arguments ..
*     ..
      INTEGER            DESCA(*)
	DOUBLE PRECISION   U(LDUV,*), V(LDUV, *)

      INTEGER            ICTXT, MYROW, MYCOL, NPROW, NPCOL 

      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
      INTEGER            ROWWISE, COLUMNWISE, ALL
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2, ALL=3)

	
	INTEGER			   NSBULGE, WSIZE, WSIZE2
      INTEGER            J, JJ, IACOL, IAROW, IAROW2, IACOL2
	INTEGER			   IC, IR

*     ..
*     .. External Functions ..
*     ..
      INTEGER            INDXG2L, INDXG2P
      EXTERNAL           INDXG2L, INDXG2P



	NSBULGE = 3

	IF (IFIRST.GE.ILAST) RETURN
	IF (J.GE.ILAST) RETURN
	ICTXT = DESCA(CTXT_)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )  
      IF (IFIRST+1.EQ.ILAST) THEN
         RETURN
      END IF
      IF (J+1.EQ.ILAST) THEN
         RETURN
      END IF

	    IAROW = INDXG2P(J, NB, 0, 0, NPROW)
	    IACOL = INDXG2P(J, NB, 0, 0, NPCOL)
		WSIZE = ((NSBULGE+1 + DIST) + 2) * IBULGE + 2
		WSIZE = MIN(WSIZE, ILAST-J+1)
		WSIZE = MIN(WSIZE, NB)
		WSIZE2 = MIN(NB, ILAST-J+1)
		WSIZE2 = MIN(WSIZE2, WSIZE+MOVEDIST)
	    IACOL2 = INDXG2P(J+WSIZE2-1, NB, 0, 0, NPCOL)
		IAROW2 = INDXG2P(J+WSIZE2-1, NB, 0, 0, NPROW)





C		Distribute U (=Q) and V(=Z)
C	    U both directions(needed for Q accumulations), V only columnwize


		IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
			
			IF (NPCOL.GT.1) THEN
				Call BSend(WSIZE2, WSIZE2, ROWWISE, U,
     $				LDUV, ICTXT)
			END IF

		ELSE IF (MYROW.EQ.IAROW.AND.NPCOL.GT.1) THEN
			Call BRecv(WSIZE2,WSIZE2, ROWWISE, U,
     $				LDUV, IAROW,IACOL, ICTXT)
			
		END IF

		IF (IAROW.NE.IAROW2) THEN
		   DO 3501 IC = 1, NPCOL 
			 IF ((IC-1).EQ.IACOL) GOTO 3501
		     IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.(IC-1))) THEN
				CALL Send(WSIZE2,WSIZE2, DOWN, V, LDUV,
     $				MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				IF (ILQ) THEN
					CALL Send(WSIZE2,WSIZE2, DOWN, U, LDUV,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				END IF
			 ELSE IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.(IC-1))) THEN
				CALL Recv(WSIZE2,WSIZE2, UP, V, LDUV,
     $				MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				IF (ILQ) THEN
					CALL Recv(WSIZE2,WSIZE2, UP, U, LDUV,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				END IF
			 END IF
 3501		   CONTINUE
	    END IF

		IF (IACOL.NE.IACOL2) THEN
		   DO 3601 IR = 1, NPROW 
			 IF ((IR-1).EQ.IAROW) GOTO 3601
		     IF ((MYCOL.EQ.IACOL).AND.(MYROW.EQ.(IR-1))) THEN
				CALL Send(WSIZE2,WSIZE2, RIGHT, U, LDUV,
     $				MYROW, MYCOL, NPCOL, NPROW, ICTXT)
			 ELSE IF ((MYCOL.EQ.IACOL2).AND.(MYROW.EQ.(IR-1))) THEN
				CALL Recv(WSIZE2,WSIZE2, LEFT, U, LDUV,
     $				MYROW, MYCOL, NPCOL, NPROW, ICTXT)

			 END IF
 3601		   CONTINUE
	    END IF 



	END



      SUBROUTINE CHASE_DISTUV_COL(
     $     IBULGE, 
     $     IFIRST, ILAST, 
     $	 NB, N,      
     $     DIST, 
     $     U, V, LDUV,
     $	 MOVEDIST,ILQ,DESCA,J)
      IMPLICIT NONE

*
*
*     .. Scalar Arguments ..
      INTEGER            IBULGE, DIST
      INTEGER            IFIRST, ILAST
      INTEGER            NB, N, LDIST,MOVEDIST
	INTEGER            LDTMP,LDUV
      LOGICAL            ILQ
*     ..
*     .. Array Arguments ..
*     ..
      INTEGER            DESCA(*)
	DOUBLE PRECISION   U(LDUV,*), V(LDUV, *)

      INTEGER            ICTXT, MYROW, MYCOL, NPROW, NPCOL 

      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
      INTEGER            ROWWISE, COLUMNWISE, ALL
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2, ALL=3)

	
	INTEGER			   NSBULGE, WSIZE, WSIZE2
      INTEGER            J, JJ, IACOL, IAROW, IAROW2, IACOL2
	INTEGER			   IC, IR

*     ..
*     .. External Functions ..
*     ..
      INTEGER            INDXG2L, INDXG2P
      EXTERNAL           INDXG2L, INDXG2P



	NSBULGE = 3

	IF (IFIRST.GE.ILAST) RETURN
	IF (J.GE.ILAST) RETURN

	ICTXT = DESCA(CTXT_)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )  
      IF (IFIRST+1.EQ.ILAST) THEN
         RETURN
      END IF
      IF (J+1.EQ.ILAST) THEN
         RETURN
      END IF

	    IAROW = INDXG2P(J, NB, 0, 0, NPROW)
	    IACOL = INDXG2P(J, NB, 0, 0, NPCOL)
		WSIZE = ((NSBULGE+1 + DIST) + 2) * IBULGE + 2
		WSIZE = MIN(WSIZE, ILAST-J+1)
		WSIZE = MIN(WSIZE, NB)
		WSIZE2 = MIN(NB, ILAST-J+1)
		WSIZE2 = MIN(WSIZE2, WSIZE+MOVEDIST)
	    IACOL2 = INDXG2P(J+WSIZE2-1, NB, 0, 0, NPCOL)
		IAROW2 = INDXG2P(J+WSIZE2-1, NB, 0, 0, NPROW)





C		Distribute U (=Q) and V(=Z)
C	    U both directions(needed for Q accumulations), V only columnwize


		IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
			
			IF (ILQ.AND.NPROW.GT.1) THEN
     				Call BSend(WSIZE2, WSIZE2, COLUMNWISE, U,
     $				LDUV, ICTXT)
			END IF

			IF (NPROW.GT.1) THEN
     				Call BSend(WSIZE2,WSIZE2, COLUMNWISE, V,
     $				LDUV, ICTXT)
			END IF
			
		ELSE IF (MYCOL.EQ.IACOL.AND.NPROW.GT.1) THEN
			IF (ILQ) THEN
				Call BRecv(WSIZE2,WSIZE2, COLUMNWISE, U,
     $				LDUV, IAROW,IACOL, ICTXT)
			END IF
			Call BRecv(WSIZE2,WSIZE2, COLUMNWISE, V,
     $				LDUV, IAROW,IACOL, ICTXT)
		END IF



	END



      SUBROUTINE CHASESHIFTS(A, B, Q, Z,
     $     DESCA, DESCB, DESCQ, DESCZ,
     $     LDA, LDB, LDQ, LDZ,
     $     TMPW, TMPW2,
     $	 LDTMP,
     $     IBULGE, NBULGE, 
     $     IFIRST, ILAST, ILASTM, IFRSTM,
     $	 SAFMIN, ILO, NB, N,      
     $     ILQ, ILZ,
     $     DIST, LDIST, TALPHAR, TALPHAI, TBETA,
     $     WORK, JITER,ATOL,U, V, LDUV,
     $	 TT1, TT2, TT3, TT4, TT5,MOVEDIST,J)
      IMPLICIT NONE
*
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDQ, LDZ
      INTEGER            IBULGE, NBULGE, DIST
      INTEGER            IFIRST, ILAST, ILASTM, IFRSTM
      INTEGER            ILO, NB, N, JITER, LDIST,MOVEDIST
      LOGICAL            ILQ, ILZ
	INTEGER            UCNT, DCNT,LDTMP,LDUV
	DOUBLE PRECISION   TT1, TT2, TT3, TT4, TT5
	DOUBLE PRECISION   SAFMIN,ATOL
*     ..
*     .. Array Arguments ..
*     ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
      DOUBLE PRECISION   Q( LDQ, * ), Z( LDZ, * )
      INTEGER            DESCA(*), DESCB(*)
      INTEGER			   DESCQ(*), DESCZ(*)
      DOUBLE PRECISION   WORK( * )
      DOUBLE PRECISION   TMPW(LDTMP, *), TMPW2(LDTMP, *)
	DOUBLE PRECISION   U(LDUV,*), V(LDUV, *)
	DOUBLE PRECISION   TALPHAI(*), TALPHAR(*), TBETA(*)
*
*     Purpose
*     =======
*
*     CHASESHIFTS chases shifts numstep steps downwards the diagonal.
*     
*
*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )

      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
      INTEGER            ROWWISE, COLUMNWISE, ALL
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2, ALL=3)


*     ..
*     .. Local Arrays ..
*     ..      
      DOUBLE PRECISION   CS((NBULGE+1)*4)
	DOUBLE PRECISION   VV(4)
	
	LOGICAL            NDONE(NBULGE)
	LOGICAL            NDONE2(NBULGE)


*     ..
*     .. Local Scalars ..
*     ..
      INTEGER            J, JJ, IACOL, IAROW, IAROW2, IACOL2
	INTEGER			   OLDIAROW, OLDIACOL, IC, IR
      INTEGER            IB, RT, LT, RET, NNB, JI, SCNT
      INTEGER            ICTXT, MYROW, MYCOL, NPROW, NPCOL 
      DOUBLE PRECISION   D1, D2, D3, D4, D5, C, S, TEMP, TAU
	DOUBLE PRECISION   TAUU, TAUV, V1(3), V2(3)
	DOUBLE PRECISION   T1, T2

	INTEGER            I, INDX, NCOLS, JC, JR, DISPL,NSBULGE
	INTEGER			   LROW1, LCOL1, LROW2, LCOL2,ICHASE
	INTEGER            INDJ1, INDJ2, INDJ3, INDJ4, INDJ5
	INTEGER			   INDJ6, INDJ7, INDJ8, RSCNT, CSCNT
	INTEGER			   F,WSIZE, ICOL, IROW, NCHASE,ORIGJ
	INTEGER			   WSIZE2, DEBUG1, LR,LC
	LOGICAL			   SENTV, INCINDJ1, INCINDJ2, INCINDJ5
	EXTERNAL		   MPI_WTIME
	DOUBLE PRECISION   MPI_WTIME
*     ..
*     .. External Subroutines ..
*     ..
      EXTERNAL           BLACS_BARRIER, PDLACP3

*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT, MOD
*     ..
*     .. External Functions ..
*     ..
      INTEGER            INDXG2L, INDXG2P
      EXTERNAL           INDXG2L, INDXG2P

	IF (IFIRST.GE.ILAST) RETURN
	IF (J.GE.ILAST) RETURN

	ICTXT = DESCA(CTXT_)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )  
	NSBULGE = 3
      LT = 1
      RT = NBULGE*3*(NSBULGE + 2 + DIST+NB) + NBULGE*2
	DEBUG1 = 0
      IF (IFIRST+1.EQ.ILAST) THEN
         RETURN
      END IF
      IF (J+1.EQ.ILAST) THEN
         RETURN
      END IF		
		
      DO 10 IB = 1, NBULGE
         NDONE(IB) = .TRUE.
         NDONE2(IB) = .TRUE.
 10   CONTINUE






	    DISPL = 0

 1000     CONTINUE

		ORIGJ = J
	    T1 = MPI_WTIME()




	    IAROW = INDXG2P(J, NB, 0, 0, NPROW)
	    IACOL = INDXG2P(J, NB, 0, 0, NPCOL)
		WSIZE = ((NSBULGE+1 + DIST) + 2) * IBULGE + 2
		WSIZE = MIN(WSIZE, ILAST-J+1)
		WSIZE = MIN(WSIZE, NB)
		WSIZE2 = MIN(NB, ILAST-J+1)
		WSIZE2 = MIN(WSIZE2, WSIZE+MOVEDIST)
	    IACOL2 = INDXG2P(J+WSIZE2-1, NB, 0, 0, NPCOL)
		IAROW2 = INDXG2P(J+WSIZE2-1, NB, 0, 0, NPROW)

*		IF (J.NE.IFIRST) THEN

		IF (WSIZE2.GT.MOD(J+WSIZE2-1,NB)) THEN
			IF (DEBUG1.EQ.1.AND.
     $			MYROW.EQ.IAROW.AND.MYCOL.EQ.IACOL) THEN
				WRITE(*,*)'Not local'
			END IF
			CALL PDLACP3(WSIZE2,J, A, DESCA, 
     $			TMPW, LDTMP, IAROW,IACOL, 0)
			CALL PDLACP3(WSIZE2, J, B, DESCB, 
     $			TMPW2, LDTMP, IAROW, IACOL, 0)
		
		ELSE IF (MYROW.EQ.IAROW.AND.MYCOL.EQ.IACOL) THEN
		 LR = INDXG2L(J,NB, 0,0,NPROW)
		 LC = INDXG2L(J,NB, 0,0,NPCOL)
		 IF (DEBUG1.EQ.1) THEN
			WRITE(*,*)'Local'
		 END IF
		 CALL DLACPY('All',WSIZE2,WSIZE2,A(LR,LC),
     $		LDA,TMPW,LDTMP)
		 CALL DLACPY('All',WSIZE2,WSIZE2,B(LR,LC),
     $		LDB,TMPW2,LDTMP)
		
		END IF
*		END IF

		CALL DLASET( 'All', LDUV, LDUV, ZERO, ONE, U,
     $                      LDUV )
		CALL DLASET( 'All', LDUV, LDUV, ZERO, ONE, V,
     $                      LDUV )
         ICHASE = 0
	   NCHASE = 1

         T2 = MPI_WTIME()
	   TT1 =TT1 + (T2 - T1)	
		




1100		CONTINUE


     

	   

1200	   CONTINUE     
         T1 = MPI_WTIME()

	   DO 2000 IB = IBULGE, 1, -1

		  NCOLS = 3		  		  
		  JJ = J + (IB-1)*(NSBULGE+1 + DIST)
		  INDX = 3*(NSBULGE+DIST+2+MAX(NB,0))*(IB-1)

            IF (.NOT. NDONE(IB).OR.JJ.GT.(ILAST-2)) THEN  
		    NDONE(IB) = .FALSE.                   									  
			GOTO 2000
            END IF	
            
		              
		  NNB = MIN(NSBULGE+DIST+NCOLS+2+MOVEDIST,ILAST-JJ+1)
            IF (NNB .LE. 0) GOTO 2000
            
*           Last block ?
            IF ((JJ + NNB -1) .EQ. ILAST) THEN
               NCOLS = 2
               NDONE(IB) = .FALSE.
            ENDIF
	        
	      IF (MYROW.EQ.IAROW.AND.MYCOL.EQ.IACOL) THEN

		  IR = (IB-1+ICHASE)*(NSBULGE+1 + DIST)+1
		  IC = (IB-1+ICHASE)*(NSBULGE+1 + DIST)+1

		  CALL QZ2Local( NNB, 
     $			TMPW(IR,IC),
     $			LDTMP,  
     $			TMPW2(IR,IC), 
     $			LDTMP,
     $			WORK(LT+INDX), WORK(RT+INDX), NCOLS, IB,
     $			WSIZE2-IR+1,
     $			ATOL, NB, (IBULGE-IB+1)*2,TALPHAI, TALPHAR, TBETA)
        
		 
*		  .. Some column update is required before moving the next bulge
		  IF ((IB+ICHASE).GT.1) THEN            
			IF ((JJ-IFRSTM.GT.0).AND.(JJ-(NSBULGE+1+DIST).GT.0)) THEN

*			13, 9, 5, 1
              IR = (IB-1+ICHASE)*(NSBULGE+1+DIST)+1
*			12, 8, 4, 0
			IC = (IB-1+ICHASE)*(NSBULGE+1+DIST)
			JI = 1
			CALL MYDLAREF( 'Cols', TMPW, 
     $			LDTMP, .FALSE., U, LDUV,
     $			.TRUE.,
     $			IR,IC+2,
     $			1,NNB-NCOLS-1,
     $            JI,IC, 
     $			0, 0,
     $            WORK(RT+INDX), D1, D2, D3, D4, D5 )	
			IC = (IB-1+ICHASE)*(NSBULGE+1 + DIST)
			JI = 1
     	  		CALL MYDLAREF( 'Cols', TMPW2, 
     $			LDTMP, .FALSE., V, LDUV,
     $			.TRUE.,
     $			IR,IC+2,
     $			1,NNB-NCOLS-1,
     $            JI,IC, 
     $			0, 0,
     $            WORK(RT+INDX), D1, D2, D3, D4, D5 )	
		
              END IF
		  END IF
           
		  
		  	  
*		  .. Move the current bulge over the edge 
*		  .. only if the last bulge is not done
*		  .. If the last bulge is done then logic will take of this futher down
		  IF ((NCOLS.EQ.2).AND.(NDONE(1))) THEN

              JJ = ILAST - 1
               
			IR = WSIZE2 -2  
			IC = WSIZE2 -2

			TEMP = TMPW( IR+1, IC )
			CALL DLARTG( TEMP,TMPW(IR+2,IC),C,S,TMPW(IR+1,IC ) )
			TMPW(IR+2,IC) = ZERO
			CS(1+(IB-1)*4) = C
			CS(2+(IB-1)*4) = S	
			
*				Update rows of A				
			IR = WSIZE2 -1  
			IC = WSIZE2 -1
			CALL MYROT( 'R', TMPW, LDTMP, 
     $				IR, IC, 1,1,
     $                IC, IC+1, CS(1+(IB-1)*4))

*				Update rows of B
			CALL MYROT( 'R', TMPW2, LDTMP, 
     $				IR, IC, 1,1,
     $                IC, IC+1, CS(1+(IB-1)*4))

			IR = WSIZE2 -2  
			IC = WSIZE2 -2
				
							
			TEMP = TMPW2(IR+2,IC+2)
			CALL DLARTG(TEMP,TMPW2(IR+2,IC+1),C,S,TMPW2(IR+2,IC+2))
			TMPW2(IR+2,IC+1) = ZERO
			CS(3+(IB-1)*4) = C
			CS(4+(IB-1)*4) = S

*			Update cols of A				
			IR = WSIZE2 -2  
			IC = WSIZE2 -1
			CALL MYROT( 'C', TMPW, LDTMP, 
     $				IR, IC, 1,1,
     $                1, WSIZE2, CS(3+(IB-1)*4))

*			Update col of B				
			IR = WSIZE2 -2  
			IC = WSIZE2 -1
			CALL MYROT( 'C', TMPW2, LDTMP, 
     $				IR, IC, 1,1,
     $                1, WSIZE2-1, CS(3+(IB-1)*4))
			END IF

			
			
		 
			END IF
 2000    CONTINUE


*	   .. Apply row and column transformations to U and V

 2090	   CONTINUE
                


         DO 2100 IB = IBULGE, 1, -1

		  NCOLS = 3		  		  
		  JJ = J + (IB-1)*(NSBULGE+1+DIST)
		  INDX = 3*(NSBULGE+DIST+2+MAX(NB,0))*(IB-1)
		  
            IF (.NOT. NDONE2(IB).OR.JJ.GT.(ILAST-2)) THEN 
			NDONE2(IB) = .FALSE.
              GOTO 2100
            END IF	
         
		  NNB = MIN(NSBULGE+DIST+NCOLS+2+MOVEDIST,ILAST-JJ+1)

            IF (NNB .LE. 0) GOTO 2100
            
*           Last block ?
            IF ((JJ + NNB -1) .EQ. ILAST) THEN
               NCOLS = 2
               NDONE2(IB) = .FALSE.
            ENDIF

		  	
            IF (MYROW.EQ.IAROW.AND.MYCOL.EQ.IACOL) THEN

*		  Apply rowrotations to U and V 
		  DO 2221 IR = 1, NNB-NCOLS-1
				VV(1) = ONE
				VV(2) = WORK(LT+INDX + (IR-1)*3)
				VV(3) = WORK(LT+INDX + (IR-1)*3+1)				
				TAU = WORK(LT+INDX + (IR-1)*3+2)
 
                  CALL DLARFX( 'Right', LDUV, 3, VV, TAU,
     $                   U(1,IR+ (IB-1+ICHASE)*(NSBULGE+1+DIST)+1),
     $					LDUV, WORK)
				VV(1) = ONE
				VV(2) = WORK(RT+INDX + (IR-1)*3)
				VV(3) = WORK(RT+INDX + (IR-1)*3+1) 				
				TAU = WORK(RT+INDX + (IR-1)*3+2)

                  CALL DLARFX( 'Right', LDUV, 3, VV, TAU,
     $                   V(1,IR+ (IB-1+ICHASE)*(NSBULGE+1+DIST)+1), 
     $					LDUV, WORK)
 2221		  CONTINUE



*		  Clean up at border		
		  IF ((NCOLS.EQ.2).AND.(NDONE2(1))) THEN
			IR = WSIZE2 -2  
			IC = WSIZE2 -1
			CALL MYROT( 'X', U, LDUV, 
     $				IR, IC, 1,1,
     $				1, LDUV, CS(1+(IB-1)*4))

			
			IR = WSIZE2 -2  
			IC = WSIZE2 -1
			CALL MYROT( 'C', V, LDUV, 
     $				IR, IC, 1,1,
     $				1, LDUV, CS(3+(IB-1)*4))
		  END IF
		  END IF
 2100    CONTINUE
 2101		CONTINUE
		T2 = MPI_WTIME()
	    TT2 = TT2 + (T2 - T1)
		
		ICHASE= ICHASE + 1
		IF (ICHASE.LT.NCHASE) THEN
			J = J + (NSBULGE+1 + DIST+MAX(NB-WSIZE-1,0))
			GOTO 1200
		END IF

		T1 = MPI_WTIME()


C	Restore A and B
		IF (WSIZE2.GT.MOD(ORIGJ+WSIZE2-1,NB)) THEN

		CALL PDLACP3(WSIZE2, ORIGJ, A, DESCA, 
     $		TMPW, LDTMP, IAROW, IACOL, -1)
		CALL PDLACP3(WSIZE2, ORIGJ, B, DESCB, 
     $		TMPW2, LDTMP, IAROW, IACOL, -1)

		ELSE IF (MYROW.EQ.IAROW.AND.MYCOL.EQ.IACOL) THEN

		 LR = INDXG2L(ORIGJ,NB, 0,0,NPROW)
		 LC = INDXG2L(ORIGJ,NB, 0,0,NPCOL)
		 
		 CALL DLACPY('All',WSIZE2,WSIZE2,TMPW,LDTMP,
     $		A(LR,LC),LDA)
		 CALL DLACPY('All',WSIZE2,WSIZE2,TMPW2,LDTMP,
     $		B(LR,LC),LDB)
		
		END IF



 2200		CONTINUE

 2201		CONTINUE




		T2 = MPI_WTIME()
	    TT3 = TT3 + (T2 - T1)


		

 3000 CONTINUE


      RETURN
*
*     End of CHASESHIFTS
*
 800  FORMAT(40f7.2)
 850  FORMAT(40F7.2)
      END


      SUBROUTINE CHASES_UPDATE(A, B, Q, Z,
     $     DESCA, DESCB, DESCQ, DESCZ,
     $     LDA, LDB, LDQ, LDZ,
     $     TMPW, TMPW2,
     $	 LDTMP,
     $     IBULGE, NBULGE, 
     $     IFIRST, ILAST, ILASTM, IFRSTM,
     $	 SAFMIN, ILO, NB, N,      
     $     ILQ, ILZ,
     $     DIST, LDIST, TALPHAR, TALPHAI, TBETA,
     $     WORK, JITER,ATOL,U, V, LDUV,
     $	 TT1, TT2, TT3, TT4, TT5,MOVEDIST,J)
      IMPLICIT NONE
*
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDQ, LDZ
      INTEGER            IBULGE, NBULGE, DIST
      INTEGER            IFIRST, ILAST, ILASTM, IFRSTM
      INTEGER            ILO, NB, N, JITER, LDIST,MOVEDIST
      LOGICAL            ILQ, ILZ
	INTEGER            UCNT, DCNT,LDTMP,LDUV
	DOUBLE PRECISION   TT1, TT2, TT3, TT4, TT5
	DOUBLE PRECISION   SAFMIN,ATOL
*     ..
*     .. Array Arguments ..
*     ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
      DOUBLE PRECISION   Q( LDQ, * ), Z( LDZ, * )
      INTEGER            DESCA(*), DESCB(*)
      INTEGER			   DESCQ(*), DESCZ(*)
      DOUBLE PRECISION   WORK( * )
      DOUBLE PRECISION   TMPW(LDTMP, *), TMPW2(LDTMP, *)
	DOUBLE PRECISION   U(LDUV,*), V(LDUV, *)
	DOUBLE PRECISION   TALPHAI(*), TALPHAR(*), TBETA(*)
*
*     Purpose
*     =======
*
*     CHASESHIFTS chases shifts numstep steps downwards the diagonal.
*     
*
*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )

      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
      INTEGER            ROWWISE, COLUMNWISE, ALL
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2, ALL=3)


*     ..
*     .. Local Arrays ..
*     ..      
      DOUBLE PRECISION   CS((NBULGE+1)*4)
	DOUBLE PRECISION   VV(4)
	
	LOGICAL            NDONE(NBULGE)
	LOGICAL            NDONE2(NBULGE)


*     ..
*     .. Local Scalars ..
*     ..
      INTEGER            J, JJ, IACOL, IAROW, IAROW2, IACOL2
	INTEGER			   OLDIAROW, OLDIACOL, IC, IR
      INTEGER            IB, RT, LT, RET, NNB, JI, SCNT
      INTEGER            ICTXT, MYROW, MYCOL, NPROW, NPCOL 
      DOUBLE PRECISION   D1, D2, D3, D4, D5, C, S, TEMP, TAU
	DOUBLE PRECISION   TAUU, TAUV, V1(3), V2(3)
	DOUBLE PRECISION   T1, T2

	INTEGER            I, INDX, NCOLS, JC, JR, DISPL,NSBULGE
	INTEGER			   LROW1, LCOL1, LROW2, LCOL2,ICHASE
	INTEGER            INDJ1, INDJ2, INDJ3, INDJ4, INDJ5
	INTEGER			   INDJ6, INDJ7, INDJ8, RSCNT, CSCNT
	INTEGER			   F,WSIZE, ICOL, IROW, NCHASE,ORIGJ
	INTEGER			   WSIZE2, DEBUG1, LR,LC
	LOGICAL			   SENTV, INCINDJ1, INCINDJ2, INCINDJ5
	EXTERNAL		   MPI_WTIME
	DOUBLE PRECISION   MPI_WTIME
*     ..
*     .. External Subroutines ..
*     ..
      EXTERNAL           BLACS_BARRIER, PDLACP3

*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT, MOD
*     ..
*     .. External Functions ..
*     ..
      INTEGER            INDXG2L, INDXG2P
      EXTERNAL           INDXG2L, INDXG2P

	IF (IFIRST.GE.ILAST) RETURN
	IF (J.GE.ILAST) RETURN

	ICTXT = DESCA(CTXT_)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )  
	NSBULGE = 3
      LT = 1
      RT = NBULGE*3*(NSBULGE + 2 + DIST+NB) + NBULGE*2
	DEBUG1 = 0
      IF (IFIRST+1.EQ.ILAST) THEN
         RETURN
      END IF
      IF (J+1.EQ.ILAST) THEN
         RETURN
      END IF
		

		ORIGJ = J
	    T1 = MPI_WTIME()




	    IAROW = INDXG2P(J, NB, 0, 0, NPROW)
	    IACOL = INDXG2P(J, NB, 0, 0, NPCOL)
		WSIZE = ((NSBULGE+1 + DIST) + 2) * IBULGE + 2
		WSIZE = MIN(WSIZE, ILAST-J+1)
		WSIZE = MIN(WSIZE, NB)
		WSIZE2 = MIN(NB, ILAST-J+1)
		WSIZE2 = MIN(WSIZE2, WSIZE+MOVEDIST)
	    IACOL2 = INDXG2P(J+WSIZE2-1, NB, 0, 0, NPCOL)
		IAROW2 = INDXG2P(J+WSIZE2-1, NB, 0, 0, NPROW)



		CALL UPDABQZ(A, B, Q, Z,
     $			LDA, LDB, LDQ, LDZ,
     $			TMPW, TMPW2,
     $			LDTMP,
     $			U, V,
     $			LDUV, LDUV,
     $			DESCA, DESCB, DESCQ, DESCZ,
     $			ORIGJ, ILASTM,
     $			ILO, NB, N, WSIZE2,
     $			ILQ, ILZ,WSIZE2)	 	

		T2 = MPI_WTIME()
		TT4 = TT4 + (T2 - T1)
		T1 = MPI_WTIME()

		JJ = J + (NSBULGE+1+DIST+MOVEDIST) 

	    IF (LDIST.GT.0.AND.JJ.LE.(ILAST)) RETURN
		IF (JJ.LE.ILAST) RETURN
			
		T1 = MPI_WTIME()
		JJ = ILAST - 1
		
		IAROW = INDXG2P(JJ, NB, 0, 0, NPROW)
		IACOL = INDXG2P(JJ, NB, 0, 0, NPCOL)
	    IACOL2 = INDXG2P(JJ+1, NB, 0, 0, NPCOL)
		IAROW2 = INDXG2P(JJ+1, NB, 0, 0, NPROW)

		CALL PDLACP3(3, JJ-1, A, DESCA, TMPW, LDTMP, IAROW, IACOL, 0)
		CALL PDLACP3(2, JJ, B, DESCB, TMPW2, LDTMP, IAROW, IACOL, 0)
		IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
			TEMP = TMPW( 2, 1 )
			CALL DLARTG( TEMP, TMPW( 3, 1 ), C, S, TMPW( 2, 1 ) )
			TMPW( 3, 1 ) = ZERO
	        CS(1) = C
		    CS(2) = S	
*			Update rows of A
			CALL MYROT( 'R', TMPW, LDTMP, 
     $				2, 2, 1,1,
     $				2, 3, CS(1))

*			Update rows of B
 			CALL MYROT( 'R', TMPW2, LDTMP, 
     $				1, 1, 1,1,
     $				1, 2, CS(1))

							
			TEMP = TMPW2( 2, 2 )
			CALL DLARTG( TEMP, TMPW2( 2, 1 ), C, S, TMPW2( 2, 2 ))
			TMPW2( 2, 1 ) = ZERO
			CS(3) = C
			CS(4) = S

*			Update cols of A
			CALL MYROT( 'C', TMPW, LDTMP, 
     $				2, 2, 1,1,
     $				2, 3, CS(3))

*			Update col of B
			CALL MYROT( 'C', TMPW2, LDTMP, 
     $				1, 1, 1,1,
     $				1, 1, CS(3))


		END IF
		
		CALL PDLACP3(3, JJ-1, A, DESCA, TMPW, LDTMP, IAROW, IACOL, -1)
		CALL PDLACP3(2, JJ, B, DESCB, TMPW2, LDTMP, IAROW, IACOL, -1)
		

*          ..
*          .. Bcast Rotations from the left.and right
*          ..		   
		IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN

			Call BSend(4, 1, ROWWISE, CS, 4, ICTXT)
		
			Call BSend(4, 1, COLUMNWISE, CS, 4, ICTXT)

	    ELSE IF (MYROW.EQ.IAROW) THEN

			Call BRecv(4, 1, ROWWISE, CS, 4, 
     $			IAROW, IACOL, ICTXT)

		ELSE IF (MYCOL.EQ.IACOL) THEN

			Call BRecv(4, 1, COLUMNWISE, CS, 4, 
     $			IAROW, IACOL, ICTXT)

		END IF
		IF (IAROW.NE.IAROW2) THEN
		   DO 3500 IC = 1, NPCOL 
			 IF ((IC-1).EQ.IACOL) GOTO 3500
		     IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.(IC-1))) THEN
					CALL Send(4, 1, DOWN, CS, 4,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
			 ELSE IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.(IC-1))) THEN
					CALL Recv(4, 1, UP, CS, 4,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
			 END IF
 3500		   CONTINUE
	    END IF

		IF (IACOL.NE.IACOL2) THEN
		   DO 3600 IR = 1, NPROW 
			 IF ((IR-1).EQ.IAROW) GOTO 3600
		     IF ((MYCOL.EQ.IACOL).AND.(MYROW.EQ.(IR-1))) THEN
					CALL Send(4, 1, RIGHT, CS, 4,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
			 ELSE IF ((MYCOL.EQ.IACOL2).AND.(MYROW.EQ.(IR-1))) THEN
					CALL Recv(4, 1, LEFT, CS, 4,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
			 END IF
 3600		   CONTINUE
	    END IF 
		
		CALL PDMYROT( 'R', A, DESCA, LDA, NB, JJ, JJ, 
     $		1, 1, JJ+2, ILASTM,
     $        CS, TMPW,LDTMP)  


		CALL PDMYROT( 'R', B, DESCB, LDB, NB, JJ, JJ, 
     $		1, 1, JJ+2, ILASTM,
     $		CS, TMPW,LDTMP)

		IF( ILQ ) THEN
			CALL PDMYROT( 'X', Q, DESCQ, LDQ, NB, JJ, JJ, 
     $			1, 1, 1, N,
     $			CS, TMPW,LDTMP)
		END IF


	    CALL PDMYROT( 'C', A, DESCA, LDA, NB, JJ, JJ, 
     $		1, 1, IFRSTM, ILAST-2,
     $		CS(3), TMPW,LDTMP)  

 
		CALL PDMYROT( 'C', B, DESCB, LDB, NB, JJ, JJ, 
     $		1, 1, IFRSTM, ILAST-2,
     $		CS(3), TMPW,LDTMP) 

		IF( ILZ ) THEN
			CALL PDMYROT( 'C', Z, DESCZ, LDZ, NB, JJ, JJ, 
     $		1, 1, 1, N,
     $		CS(3), TMPW,LDTMP)
		END IF  

		
 
		T2 = MPI_WTIME()
	    TT5 = TT5 + (T2 - T1)

		

 3000 CONTINUE



      RETURN
*
*     End of CHASESHIFTS
*
 800  FORMAT(40f7.2)
 850  FORMAT(40F7.2)
      END



      SUBROUTINE CHASES_ROUNDOF(A, B, Q, Z,
     $     DESCA, DESCB, DESCQ, DESCZ,
     $     LDA, LDB, LDQ, LDZ,
     $     TMPW, TMPW2,
     $	 LDTMP,
     $     IBULGE, NBULGE, 
     $     IFIRST, ILAST, ILASTM, IFRSTM,
     $	 SAFMIN, ILO, NB, N,      
     $     ILQ, ILZ,
     $     DIST, LDIST, TALPHAR, TALPHAI, TBETA,
     $     WORK, JITER,ATOL,U, V, LDUV,
     $	 TT1, TT2, TT3, TT4, TT5,MOVEDIST,J)
      IMPLICIT NONE
*
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDQ, LDZ
      INTEGER            IBULGE, NBULGE, DIST
      INTEGER            IFIRST, ILAST, ILASTM, IFRSTM
      INTEGER            ILO, NB, N, JITER, LDIST,MOVEDIST
      LOGICAL            ILQ, ILZ
	INTEGER            UCNT, DCNT,LDTMP,LDUV
	DOUBLE PRECISION   TT1, TT2, TT3, TT4, TT5
	DOUBLE PRECISION   SAFMIN,ATOL
*     ..
*     .. Array Arguments ..
*     ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
      DOUBLE PRECISION   Q( LDQ, * ), Z( LDZ, * )
      INTEGER            DESCA(*), DESCB(*)
      INTEGER			   DESCQ(*), DESCZ(*)
      DOUBLE PRECISION   WORK( * )
      DOUBLE PRECISION   TMPW(LDTMP, *), TMPW2(LDTMP, *)
	DOUBLE PRECISION   U(LDUV,*), V(LDUV, *)
	DOUBLE PRECISION   TALPHAI(*), TALPHAR(*), TBETA(*)
*
*     Purpose
*     =======
*
*     CHASESHIFTS chases shifts numstep steps downwards the diagonal.
*     
*
*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )

      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
      INTEGER            ROWWISE, COLUMNWISE, ALL
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2, ALL=3)


*     ..
*     .. Local Arrays ..
*     ..      
      DOUBLE PRECISION   CS((NBULGE+1)*4)
	DOUBLE PRECISION   VV(4)
	
	LOGICAL            NDONE(NBULGE)
	LOGICAL            NDONE2(NBULGE)


*     ..
*     .. Local Scalars ..
*     ..
      INTEGER            J, JJ, IACOL, IAROW, IAROW2, IACOL2
	INTEGER			   OLDIAROW, OLDIACOL, IC, IR
      INTEGER            IB, RT, LT, RET, NNB, JI, SCNT
      INTEGER            ICTXT, MYROW, MYCOL, NPROW, NPCOL 
      DOUBLE PRECISION   D1, D2, D3, D4, D5, C, S, TEMP, TAU
	DOUBLE PRECISION   TAUU, TAUV, V1(3), V2(3)
	DOUBLE PRECISION   T1, T2

	INTEGER            I, INDX, NCOLS, JC, JR, DISPL,NSBULGE
	INTEGER			   LROW1, LCOL1, LROW2, LCOL2,ICHASE
	INTEGER            INDJ1, INDJ2, INDJ3, INDJ4, INDJ5
	INTEGER			   INDJ6, INDJ7, INDJ8, RSCNT, CSCNT
	INTEGER			   F,WSIZE, ICOL, IROW, NCHASE,ORIGJ
	INTEGER			   WSIZE2, DEBUG1, LR,LC
	LOGICAL			   SENTV, INCINDJ1, INCINDJ2, INCINDJ5
	EXTERNAL		   MPI_WTIME
	DOUBLE PRECISION   MPI_WTIME
*     ..
*     .. External Subroutines ..
*     ..
      EXTERNAL           BLACS_BARRIER, PDLACP3

*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT, MOD
*     ..
*     .. External Functions ..
*     ..
      INTEGER            INDXG2L, INDXG2P
      EXTERNAL           INDXG2L, INDXG2P

	IF (IFIRST.GE.ILAST) RETURN
	IF (J.GE.ILAST) RETURN

	ICTXT = DESCA(CTXT_)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )  
	NSBULGE = 3
      LT = 1
      RT = NBULGE*3*(NSBULGE + 2 + DIST+NB) + NBULGE*2
	DEBUG1 = 0
      IF (IFIRST+1.EQ.ILAST) THEN
         RETURN
      END IF
      IF (J+1.EQ.ILAST) THEN
         RETURN
      END IF		



		T1 = MPI_WTIME()
		JJ = ILAST - 1
		
		IAROW = INDXG2P(JJ, NB, 0, 0, NPROW)
		IACOL = INDXG2P(JJ, NB, 0, 0, NPCOL)
	    IACOL2 = INDXG2P(JJ+1, NB, 0, 0, NPCOL)
		IAROW2 = INDXG2P(JJ+1, NB, 0, 0, NPROW)

		CALL PDLACP3(3, JJ-1, A, DESCA, TMPW, LDTMP, IAROW, IACOL, 0)
		CALL PDLACP3(2, JJ, B, DESCB, TMPW2, LDTMP, IAROW, IACOL, 0)
		IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
			TEMP = TMPW( 2, 1 )
			CALL DLARTG( TEMP, TMPW( 3, 1 ), C, S, TMPW( 2, 1 ) )
			TMPW( 3, 1 ) = ZERO
	        CS(1) = C
		    CS(2) = S	
*			Update rows of A
			CALL MYROT( 'R', TMPW, LDTMP, 
     $				2, 2, 1,1,
     $				2, 3, CS(1))

*			Update rows of B
 			CALL MYROT( 'R', TMPW2, LDTMP, 
     $				1, 1, 1,1,
     $				1, 2, CS(1))

							
			TEMP = TMPW2( 2, 2 )
			CALL DLARTG( TEMP, TMPW2( 2, 1 ), C, S, TMPW2( 2, 2 ))
			TMPW2( 2, 1 ) = ZERO
			CS(3) = C
			CS(4) = S

*			Update cols of A
			CALL MYROT( 'C', TMPW, LDTMP, 
     $				2, 2, 1,1,
     $				2, 3, CS(3))

*			Update col of B
			CALL MYROT( 'C', TMPW2, LDTMP, 
     $				1, 1, 1,1,
     $				1, 1, CS(3))


		END IF
		
		CALL PDLACP3(3, JJ-1, A, DESCA, TMPW, LDTMP, IAROW, IACOL, -1)
		CALL PDLACP3(2, JJ, B, DESCB, TMPW2, LDTMP, IAROW, IACOL, -1)
		

*          ..
*          .. Bcast Rotations from the left.and right
*          ..		   
		IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN

			Call BSend(4, 1, ROWWISE, CS, 4, ICTXT)
		
			Call BSend(4, 1, COLUMNWISE, CS, 4, ICTXT)

	    ELSE IF (MYROW.EQ.IAROW) THEN

			Call BRecv(4, 1, ROWWISE, CS, 4, 
     $			IAROW, IACOL, ICTXT)

		ELSE IF (MYCOL.EQ.IACOL) THEN

			Call BRecv(4, 1, COLUMNWISE, CS, 4, 
     $			IAROW, IACOL, ICTXT)

		END IF
		IF (IAROW.NE.IAROW2) THEN
		   DO 3500 IC = 1, NPCOL 
			 IF ((IC-1).EQ.IACOL) GOTO 3500
		     IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.(IC-1))) THEN
					CALL Send(4, 1, DOWN, CS, 4,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
			 ELSE IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.(IC-1))) THEN
					CALL Recv(4, 1, UP, CS, 4,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
			 END IF
 3500		   CONTINUE
	    END IF

		IF (IACOL.NE.IACOL2) THEN
		   DO 3600 IR = 1, NPROW 
			 IF ((IR-1).EQ.IAROW) GOTO 3600
		     IF ((MYCOL.EQ.IACOL).AND.(MYROW.EQ.(IR-1))) THEN
					CALL Send(4, 1, RIGHT, CS, 4,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
			 ELSE IF ((MYCOL.EQ.IACOL2).AND.(MYROW.EQ.(IR-1))) THEN
					CALL Recv(4, 1, LEFT, CS, 4,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
			 END IF
 3600		   CONTINUE
	    END IF 
		
		CALL PDMYROT( 'R', A, DESCA, LDA, NB, JJ, JJ, 
     $		1, 1, JJ+2, ILASTM,
     $        CS, TMPW,LDTMP)  


		CALL PDMYROT( 'R', B, DESCB, LDB, NB, JJ, JJ, 
     $		1, 1, JJ+2, ILASTM,
     $		CS, TMPW,LDTMP)

		IF( ILQ ) THEN
			CALL PDMYROT( 'X', Q, DESCQ, LDQ, NB, JJ, JJ, 
     $			1, 1, 1, N,
     $			CS, TMPW,LDTMP)
		END IF


	    CALL PDMYROT( 'C', A, DESCA, LDA, NB, JJ, JJ, 
     $		1, 1, IFRSTM, ILAST-2,
     $		CS(3), TMPW,LDTMP)  

 
		CALL PDMYROT( 'C', B, DESCB, LDB, NB, JJ, JJ, 
     $		1, 1, IFRSTM, ILAST-2,
     $		CS(3), TMPW,LDTMP) 

		IF( ILZ ) THEN
			CALL PDMYROT( 'C', Z, DESCZ, LDZ, NB, JJ, JJ, 
     $		1, 1, 1, N,
     $		CS(3), TMPW,LDTMP)
		END IF  

		
 
		T2 = MPI_WTIME()
	    TT5 = TT5 + (T2 - T1)

		

 3000 CONTINUE


      RETURN
*
*     End of CHASESHIFTS
*
 800  FORMAT(40f7.2)
 850  FORMAT(40F7.2)
      END
      SUBROUTINE CREATEINITBUMP(TMPA, TMPB,
     $	 ESHIFT,
     $     NB, 
     $	 TALPHAR, TALPHAI, TBETA, V, LDTMP)

      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
*     ..
      INTEGER            NB
	INTEGER            ESHIFT
*     ..
*     .. Array Arguments ..
*     ..
      DOUBLE PRECISION   TALPHAI( * ), TALPHAR( * ), TBETA( * )
	INTEGER			   LDTMP
      DOUBLE PRECISION   TMPA(LDTMP,LDTMP), TMPB(LDTMP,LDTMP)
	DOUBLE PRECISION   V( * )
	


*


*     Purpose
*     =======

*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )

      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
*     ..
      INTEGER            IACOL, IAROW, INFO
      INTEGER            ICTXT, MYROW, MYCOL, NPROW, NPCOL
	DOUBLE PRECISION   TAU


*     ..
*     .. Local Arrays
*     ..


*     ..
*     .. External Subroutines ..
*     ..
      EXTERNAL           DLARFG
*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC          ABS

	   CALL QZFCOL(3, 2, TMPA, TALPHAR(ESHIFT-1), TALPHAI(ESHIFT-1),
     $		TBETA(ESHIFT-1), TMPA, LDTMP, TMPB, LDTMP, V, INFO)
               
         CALL DLARFG( 3, V( 1 ), V( 2 ), 1, TAU )
         V( 1 ) = ONE
         V( 4 ) = TAU


	RETURN
	END
      SUBROUTINE FINDSTART(A,
     $     DESCA, LDA,
     $     IFIRST, ILAST,
     $     ILO, IHI, GOTOL,
     $     ATOL)
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
*     ..
      INTEGER          IFIRST, ILAST
      INTEGER          ILO,GOTOL, LDA,IHI
      DOUBLE PRECISION ATOL, BTOL

*     ..
*     .. Array Arguments ..
*     ..
      DOUBLE PRECISION A(LDA,*)
      INTEGER          DESCA(*)

*
*     Purpose
*     =======
*
*     FINDSTART finds a suitable place to start the iteration on. 
*     It also chases zero diag elements in B downward the diagonal.
*     
*

*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D0+0 )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
*     ..
      INTEGER            J,NH
	DOUBLE PRECISION   ALPHA,ALPHA1,ALPHA2,SMLNUM
	DOUBLE PRECISION   ULP, TST1,UNFL,OVFL


*     ..
*     .. External Subroutines ..
*     ..
      EXTERNAL           PDELGET, PDELSET, DLABAD

*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC          ABS

*     ..
*     .. External Functions ..
*     ..
	DOUBLE PRECISION	DLAMCH
	EXTERNAL			DLAMCH	


	NH = IHI-ILO+1
      ULP = DLAMCH( 'Precision' ) 
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL DLABAD( UNFL, OVFL )
      ULP = DLAMCH( 'Precision' ) 
      SMLNUM = UNFL*( NH / ULP )

	
      DO 10 J = ILAST - 1, ILO, -1
*
*        Test 1: for A(j,j-1)=0 or j=ILO
*                 
         IF( J.EQ.ILO ) THEN
		  IFIRST = J
			GO TO 110
         ELSE
		  CALL PDELGET('A', ' ', ALPHA, A, J,J-1,DESCA)
		 
		 
		  IF (ALPHA.EQ.ZERO) THEN		

			 
			 ALPHA = ZERO
* 			 CALL PDELSET(A, J,J-1, DESCA, ALPHA)
			
			 IFIRST = J
			 GO TO 110			 

		  END IF
		  CALL PDELGET('A', ' ', ALPHA1, A, J,J,DESCA)
		  CALL PDELGET('A', ' ', ALPHA2, A, J-1,J-1,DESCA)
            
		  TST1 = ABS( ALPHA1 ) + ABS( ALPHA2 )
*		  IF (TST1.EQ.ZERO) WRITE(*,*)'Zero'
	      IF ( ABS( ALPHA).LE.MAX( ULP*TST1, SMLNUM ) ) THEN
			 ALPHA = ZERO
 			 CALL PDELSET(A, J,J-1, DESCA, ALPHA)
			
			 IFIRST = J
			 GO TO 110			 


		  END IF
         END IF

*     
*        Neither test passed -- try next J
*
 10   CONTINUE
      
      GOTOL=120
      GOTO 120
      
 70   CONTINUE
      GOTOL=70
      GOTO 120
 80   CONTINUE
      GOTOL=80
      GOTO 120
      
 110  CONTINUE
      GOTOL=110
      GOTO 120
      
 120  CONTINUE
      RETURN
*     ..
*     .. END OF FINDSTART
*     ..
      END
*                     R C
      SUBROUTINE Send(M,N,DEST,A,LDA,MYROW,MYCOL,NPCOL,NPROW,CTX)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER N, M, DEST, LDA, MYROW, MYCOL, NPCOL, NPROW, CTX

*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*)

*     .. Parameters ..
      INTEGER DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )

*     .. Local scalars ..
      INTEGER DEST_ROW, DEST_COL

*     ..
*     .. Start execution
*     ..
      IF (DEST.EQ.UP) THEN
         DEST_ROW = MYROW - 1
         IF (DEST_ROW.LT.0) THEN
            DEST_ROW = NPROW - 1
         ENDIF
         DEST_COL = MYCOL
      ELSE IF (DEST.EQ.DOWN) THEN
         DEST_ROW = MYROW + 1
         IF (DEST_ROW.EQ.NPROW) THEN
            DEST_ROW = 0
         ENDIF
         DEST_COL = MYCOL    
      ELSE IF (DEST.EQ.RIGHT) THEN
         DEST_COL = MYCOL + 1
         IF (DEST_COL.EQ.NPCOL) THEN
            DEST_COL = 0
         ENDIF
         DEST_ROW = MYROW
      ELSE
         DEST_COL = MYCOL - 1
         IF (DEST_COL.LT.0) THEN
            DEST_COL = NPCOL - 1
         ENDIF
         DEST_ROW = MYROW     
      ENDIF
      CALL DGESD2D(CTX,M,N,A,LDA,DEST_ROW,DEST_COL)

      RETURN
      END

*                     R C
      SUBROUTINE Sendtmp(M,N,DEST,A,LDA,MYROW,MYCOL,NPCOL,NPROW,CTX)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER N, M, DEST, LDA, MYROW, MYCOL, NPCOL, NPROW, CTX

*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*)

*     .. Parameters ..
      INTEGER DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )

*     .. Local scalars ..
      INTEGER DEST_ROW, DEST_COL

*     ..
*     .. Start execution
*     ..
      IF (DEST.EQ.UP) THEN
         DEST_ROW = MYROW - 1
         IF (DEST_ROW.LT.0) THEN
            DEST_ROW = NPROW - 1
         ENDIF
         DEST_COL = MYCOL
      ELSE IF (DEST.EQ.DOWN) THEN
         DEST_ROW = MYROW + 1
         IF (DEST_ROW.EQ.NPROW) THEN
            DEST_ROW = 0
         ENDIF
         DEST_COL = MYCOL    
      ELSE IF (DEST.EQ.RIGHT) THEN
         DEST_COL = MYCOL + 1
         IF (DEST_COL.EQ.NPCOL) THEN
            DEST_COL = 0
         ENDIF
         DEST_ROW = MYROW
      ELSE
         DEST_COL = MYCOL - 1
         IF (DEST_COL.LT.0) THEN
            DEST_COL = NPCOL - 1
         ENDIF
         DEST_ROW = MYROW     
      ENDIF
      CALL DGESD2D(CTX,M,N,A,1,DEST_ROW,DEST_COL)

      RETURN
      END
*                     R C
      SUBROUTINE Recv(M,N,DEST,A,LDA,MYROW,MYCOL,NPCOL,NPROW,CTX)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER N, M, DEST, LDA, MYROW, MYCOL, NPCOL, NPROW, CTX

*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*)

*     .. Parameters ..
      INTEGER DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )

*     .. Local scalars ..
      INTEGER SRC_ROW, SRC_COL


*     ..
*     .. Start execution
*     ..
      IF (DEST.EQ.UP) THEN
         SRC_ROW = MYROW - 1
         IF (SRC_ROW.LT.0) THEN
            SRC_ROW = NPROW - 1
         ENDIF
         SRC_COL = MYCOL
      ELSE IF (DEST.EQ.DOWN) THEN
         SRC_ROW = MYROW + 1
         IF (SRC_ROW.EQ.NPROW) THEN
            SRC_ROW = 0
         ENDIF
         SRC_COL = MYCOL
      ELSE IF (DEST.EQ.RIGHT) THEN
         SRC_COL = MYCOL + 1
         IF (SRC_COL.EQ.NPCOL) THEN
            SRC_COL = 0
         ENDIF
         SRC_ROW = MYROW
      ELSE
         SRC_COL = MYCOL - 1
         IF (SRC_COL.LT.0) THEN
            SRC_COL = NPCOL - 1
         ENDIF
         SRC_ROW = MYROW
      ENDIF

      CALL DGERV2D(CTX,M,N,A,LDA,SRC_ROW,SRC_COL)
      RETURN
      END

*                     R C
      SUBROUTINE Recvtmp(M,N,DEST,A,LDA,MYROW,MYCOL,NPCOL,NPROW,CTX)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER N, M, DEST, LDA, MYROW, MYCOL, NPCOL, NPROW, CTX

*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*)

*     .. Parameters ..
      INTEGER DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )

*     .. Local scalars ..
      INTEGER SRC_ROW, SRC_COL


*     ..
*     .. Start execution
*     ..
      IF (DEST.EQ.UP) THEN
         SRC_ROW = MYROW - 1
         IF (SRC_ROW.LT.0) THEN
            SRC_ROW = NPROW - 1
         ENDIF
         SRC_COL = MYCOL
      ELSE IF (DEST.EQ.DOWN) THEN
         SRC_ROW = MYROW + 1
         IF (SRC_ROW.EQ.NPROW) THEN
            SRC_ROW = 0
         ENDIF
         SRC_COL = MYCOL
      ELSE IF (DEST.EQ.RIGHT) THEN
         SRC_COL = MYCOL + 1
         IF (SRC_COL.EQ.NPCOL) THEN
            SRC_COL = 0
         ENDIF
         SRC_ROW = MYROW
      ELSE
         SRC_COL = MYCOL - 1
         IF (SRC_COL.LT.0) THEN
            SRC_COL = NPCOL - 1
         ENDIF
         SRC_ROW = MYROW
      ENDIF

      CALL DGERV2D(CTX,M,N,A,1,SRC_ROW,SRC_COL)
      RETURN
      END

*                      R C
      SUBROUTINE BRecv(M,N,DIRECTION,A,LDA,SRC_ROW,SRC_COL,CTX)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER N, M, LDA, DIRECTION, SRC_ROW, SRC_COL, CTX

*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*)

*     .. Parameters ..
      INTEGER ROWWISE, COLUMNWISE, ALL
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2, ALL=3)


*     ..
*     .. Start execution
*     ..
      IF (DIRECTION.EQ.ROWWISE) THEN
         CALL DGEBR2D(CTX,'Row','1',M,N,A,LDA,SRC_ROW,SRC_COL)
      ELSE IF (DIRECTION.EQ.COLUMNWISE) THEN
         CALL DGEBR2D(CTX,'Column','1',M,N,A,LDA,SRC_ROW,SRC_COL)
      ELSE
         CALL DGEBR2D(CTX,'All','1',M,N,A,LDA,SRC_ROW,SRC_COL)
      ENDIF

      RETURN
      END

*                      R C
      SUBROUTINE BSend(M,N,DIRECTION,A,LDA,CTX)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER N, M, LDA, DIRECTION, CTX

*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*)

*     .. Parameters ..
      INTEGER ROWWISE, COLUMNWISE, ALL
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2, ALL=3)


*     ..
*     .. Start execution
*     ..
      IF (DIRECTION.EQ.ROWWISE) THEN
         CALL DGEBS2D(CTX,'Row','1',M,N,A,LDA)
      ELSE IF (DIRECTION.EQ.COLUMNWISE) THEN
         CALL DGEBS2D(CTX,'Column','1',M,N,A,LDA)
      ELSE
         CALL DGEBS2D(CTX,'All','1',M,N,A,LDA)
      ENDIF

      RETURN
      END
      SUBROUTINE INSHIFTS(A, B, Q, Z,
     $     ALPHAI, ALPHAR, BETA,
     $     LDA, LDB, LDQ, LDZ,
     $     TMPW, TMPW2,
     $	 LDTMP,
     $     DESCA, DESCB, DESCQ, DESCZ,
     $     IFRSTM, IFIRST, ILAST, ILASTM,
     $     ILO, NB, N, MAXSHIFTS, SDIST, LDIST, 
     $     NBULGES, IITER, ESHIFT, MAXIT,
     $	 ILQ, ILZ, ILSCHR, SAFMIN,
     $	 TALPHAR, TALPHAI, TBETA, WORK, JITER,
     $	 ATOL,PU,PV, LDUV)



      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
*     ..
      INTEGER            LDA, LDB, LDQ, LDZ
      INTEGER            IFRSTM, IFIRST, ILAST, ILASTM
      INTEGER            ILO, NB, N, MAXSHIFTS, SDIST, LDIST
      INTEGER            NBULGES, LDTMP
      INTEGER            IITER, MAXIT, JITER, RET
	INTEGER			   LDUV
      DOUBLE PRECISION   ESHIFT
	DOUBLE PRECISION   SAFMIN,ATOL
      LOGICAL            ILQ, ILZ, ILSCHR
*     ..
*     .. Array Arguments ..
*     ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
      DOUBLE PRECISION   Q( LDQ, * ), Z( LDZ, * )
      DOUBLE PRECISION   ALPHAI( * ), ALPHAR( * )
      DOUBLE PRECISION   WORK( * )
      DOUBLE PRECISION   TMPW(LDTMP, *), TMPW2(LDTMP,*)
	DOUBLE PRECISION   PU(LDUV,*), PV(LDUV, *)
      DOUBLE PRECISION   TALPHAI(*), TALPHAR(*), TBETA(*)
      DOUBLE PRECISION   BETA( * )
      INTEGER            DESCA(*), DESCB(*)
      INTEGER            DESCQ(*), DESCZ(*)


*


*     Purpose
*     =======
*
*     INSHIFT introduces a maximum of MAXSHIFTS in A and B.
*     If more than one shift is introduced the shifts which already
*	have been introduced are chased SDIST steps down the diagonal 
*	to make room for the new shift
*
*	We use either Double Implicit Shifts or Single Shifts dependning
*	on the size of our workarea (defined by IFIRST..ILAST). 
*	If IFIRST+1 = ILAST then we use single shifts, otherwise we use
*	double implicit shifts.
*
*	Notice that if we use single shift we do not introduce more than 
*	one shift.	
*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION   HALF, ZERO, ONE, SAFETY
      PARAMETER          ( HALF = 0.5D+0, ZERO = 0.0D+0, ONE = 1.0D+0,
     $                   SAFETY = 1.0D+2 )
      INTEGER            DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
      INTEGER            ROWWISE, COLUMNWISE, ALL
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2, ALL=3)
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
*     ..
      INTEGER            IACOL, IAROW
      INTEGER            IACOL1, IAROW1
      INTEGER            IACOL2, IAROW2, IC, IR
      INTEGER            ISHIFT, J, JJ, LR, LC
      INTEGER            SCNT, DCNT, RCNT,NN
      INTEGER            NUMSHIFTS,I
      INTEGER            IBULGE, LS, LDV
	INTEGER			   ICTXT, MYROW, MYCOL, NPROW, NPCOL

	DOUBLE PRECISION   D1, D2, D3, D4, D5
      DOUBLE PRECISION   TAU
      
	INTEGER			   LDT 
	PARAMETER		   (LDT = 4)
	DOUBLE PRECISION   TMPA(LDT,LDT), TMPB(LDT,LDT)
	
      INTEGER            INFO, WSIZE,NSBULGE

*     ..
*     .. Local Arrays
*     ..

      DOUBLE PRECISION   V(4)


*     ..
*     .. External Subroutines ..
*     ..
      EXTERNAL			PDLACP3,INVHSE,KDHGEQZ,UPDABQZ
	EXTERNAL			CHASESHIFTS,CREATEINITBUMP
	EXTERNAL			SEND,RECV, BSEND, BRECV,PDLAREF
*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC			MIN, MOD
*     ..
*     .. External Functions ..
*     ..
      INTEGER				INDXG2L, INDXG2P
      EXTERNAL			INDXG2L, INDXG2P


	ICTXT = DESCA(CTXT_)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )  

	LDV = 4
      NBULGES = 0
	NSBULGE = 3

 1	CONTINUE
      J = IFIRST

      IAROW = INDXG2P(J, NB, 0, 0, NPROW)
      IACOL = INDXG2P(J, NB, 0, 0, NPCOL)


*     .. Determine where to calc shifts from
*     .. Use a maximum of NB/2 bulges, that is NB shifts
      NUMSHIFTS = MIN(MAXSHIFTS, ILAST-J+1)

		
*	NUMSHIFTS = MIN(NUMSHIFTS, NB-MOD(J,NB)+1)
*	IF (NUMSHIFTS.LE.1) THEN
*		NUMSHIFTS = MIN(MAXSHIFTS, ILAST-J+1)
*	END IF
*     .. Total number of bulges
      NBULGES = NUMSHIFTS/2
*     .. Number of bulges going so far
      IBULGE = 0 

      IAROW = INDXG2P(J, NB, 0, 0, NPROW)
      IACOL = INDXG2P(J, NB, 0, 0, NPCOL) 

	
      CALL DLASET( 'All', LDUV, LDUV, ZERO, ONE, PU,
     $                      LDUV )
	CALL DLASET( 'All', LDUV, LDUV, ZERO, ONE, PV,
     $                      LDUV )


	WSIZE = ((NSBULGE+1 + SDIST) + 2) * NBULGES + 2

	WSIZE = MIN(WSIZE, ILAST-J+1)
	WSIZE = MIN(WSIZE, NB)
	
      IAROW2 = INDXG2P(J+WSIZE-1, NB, 0, 0, NPROW)
      IACOL2 = INDXG2P(J+WSIZE-1, NB, 0, 0, NPCOL) 

	IF (WSIZE.GT.MOD(J+WSIZE-1,NB)) THEN

		CALL PDLACP3(WSIZE, J, A, DESCA, 
     $		TMPW, LDTMP, IAROW, IACOL, 0)
		CALL PDLACP3(WSIZE, J, B, DESCB, 
     $		TMPW2, LDTMP, IAROW, IACOL, 0)           
	ELSE IF (MYROW.EQ.IAROW.AND.MYCOL.EQ.IACOL) THEN

		 LR = INDXG2L(J,NB, 0,0,NPROW)
		 LC = INDXG2L(J,NB, 0,0,NPCOL)
		 
		 CALL DLACPY('All',WSIZE,WSIZE,A(LR,LC),
     $		LDA,TMPW,LDTMP)
		 CALL DLACPY('All',WSIZE,WSIZE,B(LR,LC),
     $		LDB,TMPW2,LDTMP)
		
	END IF

	   NN = MIN(16, NB)
	   IF (((ILAST-IFIRST).LE.(NN-1))) THEN
	

	         NN = ILAST - IFIRST + 1

c			 CALL PDLACP3(NN,J,A,DESCA,TMPW,LDTMP,-1,-1,0)
c		     CALL PDLACP3(NN,J,B,DESCB,TMPW2,LDTMP,-1,-1,0)

			 
	         CALL KDHGEQZ('S','I','I',NN,1,NN, 
     $	        TMPW,LDTMP,TMPW2,LDTMP, 
     $		    ALPHAR(J),ALPHAI(J),
     $			BETA(J),
     $	        PU,LDUV,PV,LDUV, 
     $			WORK,4*NN,INFO,NB)


*		     .. Restore updated parts of A and B

			IF (WSIZE.GT.MOD(J+WSIZE-1,NB)) THEN

				CALL PDLACP3(WSIZE, J, A, DESCA, 
     $				TMPW, LDTMP, IAROW, IACOL, -1)
				CALL PDLACP3(WSIZE, J, B, DESCB, 
     $				TMPW2, LDTMP, IAROW, IACOL, -1)
           
			ELSE IF (MYROW.EQ.IAROW.AND.MYCOL.EQ.IACOL) THEN

				LR = INDXG2L(J,NB, 0,0,NPROW)
				LC = INDXG2L(J,NB, 0,0,NPCOL)
		 
				CALL DLACPY('All',WSIZE,WSIZE,TMPW,LDTMP,
     $				A(LR,LC),LDA)
				CALL DLACPY('All',WSIZE,WSIZE,TMPW2,LDTMP,
     $				B(LR,LC),LDB)
		
			END IF


			 CALL UPDABQZ(A, B, Q, Z,
     $			LDA, LDB, LDQ, LDZ,
     $			TMPW, TMPW2,
     $			LDTMP,
     $			PU, PV,
     $			LDUV, LDUV,
     $			DESCA, DESCB, DESCQ, DESCZ,
     $			IFIRST, ILASTM,
     $			ILO, NB, N, NN,
     $			ILQ, ILZ,NN)	 				



			 ILAST = ILAST - NN
	         GOTO 1000
	   END IF


*     .. Gather elements to introduce shifts
	
      DO 10 ISHIFT=NUMSHIFTS, 1, -2	
		IF (IBULGE.GE.NBULGES) GOTO 10
		  
* 	   JJ = NUMSHIFTS - ISHIFT + 2
	   JJ = ISHIFT+1
	   IF (MYROW.EQ.IAROW.AND.MYCOL.EQ.IACOL) THEN
*
*	.. Introduce a bump
*	   
		CALL CREATEINITBUMP(TMPW,TMPW2, 
     $			JJ,
     $			LDTMP,
     $			TALPHAR, TALPHAI, TBETA, V,LDTMP)

	    TAU = V(4) 


*	    Update A, B and PU
		CALL DLARFX( 'Left', 3, WSIZE, V, TAU, 
     $             TMPW, LDTMP, WORK )		
		CALL DLARFX( 'Left', 3, WSIZE, V, TAU, 
     $             TMPW2, LDTMP, WORK )
          CALL DLARFX( 'Right', LDUV, 3, V, TAU,
     $        PU, LDUV, WORK )


*         Zero j-th column of B
*	   
          CALL INVHSE(3, TMPW2, TMPW2, LDTMP, TAU,V,INFO)	
		V(4) = TAU

		
*	    Update A, B and PV
          CALL DLARFX( 'Right', 3, 3, V, TAU,
     $		TMPW2, LDTMP, WORK )
		TMPW2(2,1) = ZERO
		TMPW2(3,1) = ZERO

          CALL DLARFX( 'Right', 4, 3, V, TAU,
     $         TMPW, LDTMP, WORK)
          CALL DLARFX( 'Right', LDUV, 3, V, TAU,
     $         PV, LDUV, WORK )


	   END IF
*        .. Number of bulges going so far
         IBULGE = IBULGE + 1
	   IF (MYROW.EQ.IAROW.AND.MYCOL.EQ.IACOL) THEN
	   IF (IBULGE.NE.NBULGES) THEN
*			 .. Hunt bulges one step forward to prepare for next
	        CALL CHASESINTRO(A, B, Q, Z,
     $	        DESCA, DESCB, DESCQ, DESCZ,
     $		    LDA, LDB, LDQ, LDZ,
     $			TMPW, TMPW2, LDTMP,
     $	        IBULGE, NBULGES,
     $			IFIRST, ILAST, ILASTM, IFRSTM,
     $			SAFMIN, ILO, NB, N,
     $			ILQ, ILZ, SDIST, LDIST,
     $			TALPHAR, TALPHAI, TBETA, 
     $			WORK, JITER,ATOL,PU,PV, LDUV,
     $			ZERO,ZERO,ZERO,ZERO,ZERO)			
			
	   END IF

	   END IF


 10   CONTINUE

C	Restore A and B
	IF (WSIZE.GT.MOD(J+WSIZE-1,NB)) THEN

		CALL PDLACP3(WSIZE, J, A, DESCA, 
     $		TMPW, LDTMP, IAROW, IACOL, -1)
		CALL PDLACP3(WSIZE, J, B, DESCB, 
     $		TMPW2, LDTMP, IAROW, IACOL, -1)
           
	ELSE IF (MYROW.EQ.IAROW.AND.MYCOL.EQ.IACOL) THEN

		 LR = INDXG2L(J,NB, 0,0,NPROW)
		 LC = INDXG2L(J,NB, 0,0,NPCOL)
		 
		 CALL DLACPY('All',WSIZE,WSIZE,TMPW,LDTMP,
     $		A(LR,LC),LDA)
		 CALL DLACPY('All',WSIZE,WSIZE,TMPW2,LDTMP,
     $		B(LR,LC),LDB)
		
	END IF


C		Distribute U (=Q) and V(=Z)
C	    U both directions(needed for Q accumulations), V only columnwize


		IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
			
			IF (ILQ.AND.NPROW.GT.1) THEN
     				Call BSend(WSIZE, WSIZE, COLUMNWISE, PU,
     $				LDUV, ICTXT)
			END IF
			IF (NPCOL.GT.1) THEN
				Call BSend(WSIZE, WSIZE, ROWWISE, PU,
     $				LDUV, ICTXT)
			END IF

			IF (NPROW.GT.1) THEN
     				Call BSend(WSIZE,WSIZE, COLUMNWISE, PV,
     $				LDUV, ICTXT)
			END IF

		ELSE IF (MYROW.EQ.IAROW.AND.NPCOL.GT.1) THEN
			Call BRecv(WSIZE,WSIZE, ROWWISE, PU,
     $				LDUV, IAROW,IACOL, ICTXT)
			
		ELSE IF (MYCOL.EQ.IACOL.AND.NPROW.GT.1) THEN
			IF (ILQ) THEN
				Call BRecv(WSIZE,WSIZE, COLUMNWISE, PU,
     $				LDUV, IAROW,IACOL, ICTXT)
			END IF
			Call BRecv(WSIZE,WSIZE, COLUMNWISE, PV,
     $				LDUV, IAROW,IACOL, ICTXT)
		END IF

		IF (IAROW.NE.IAROW2) THEN
		   DO 3501 IC = 1, NPCOL 
			 IF ((IC-1).EQ.IACOL) GOTO 3501
		     IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.(IC-1))) THEN
				CALL Send(WSIZE,WSIZE, DOWN, PV, LDUV,
     $				MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				IF (ILQ) THEN
					CALL Send(WSIZE,WSIZE, DOWN, PU, LDUV,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				END IF
			 ELSE IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.(IC-1))) THEN
				CALL Recv(WSIZE,WSIZE, UP, PV, LDUV,
     $				MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				IF (ILQ) THEN
					CALL Recv(WSIZE,WSIZE, UP, PU, LDUV,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				END IF
			 END IF
 3501		   CONTINUE
	    END IF

		IF (IACOL.NE.IACOL2) THEN
		   DO 3601 IR = 1, NPROW 
			 IF ((IR-1).EQ.IAROW) GOTO 3601
		     IF ((MYCOL.EQ.IACOL).AND.(MYROW.EQ.(IR-1))) THEN
				CALL Send(WSIZE,WSIZE, RIGHT, PU, LDUV,
     $				MYROW, MYCOL, NPCOL, NPROW, ICTXT)
			 ELSE IF ((MYCOL.EQ.IACOL2).AND.(MYROW.EQ.(IR-1))) THEN
				CALL Recv(WSIZE,WSIZE, LEFT, PU, LDUV,
     $				MYROW, MYCOL, NPCOL, NPROW, ICTXT)

			 END IF
 3601		   CONTINUE
	    END IF 

*		  .. End of distribute U and V

	CALL UPDABQZ(A, B, Q, Z,
     $			LDA, LDB, LDQ, LDZ,
     $			TMPW, TMPW2,
     $			LDTMP,
     $			PU, PV,
     $			LDUV, LDUV,
     $			DESCA, DESCB, DESCQ, DESCZ,
     $			J, ILASTM,
     $			ILO, NB, N, WSIZE,
     $			ILQ, ILZ,WSIZE)
 1000 CONTINUE

      RETURN
*     End of INSHIFTS

	END 


      SUBROUTINE INVHSE( N, PARA, T, LDT, TAU, V, INFO )
      IMPLICIT NONE
C
C     Computes the first column of inv(T).
C     To be done: version based on RQ-factorization.
C     It is assumed that N <= 13. (possibly better to do TMP -> DWORK)
C
C     .. Parameters ..
      INTEGER           BULGMX
      PARAMETER         ( BULGMX = 12 )
      INTEGER           LDTMP
      PARAMETER         ( LDTMP = BULGMX+1 )
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDT, N
      DOUBLE PRECISION  TAU
C     .. Array Arguments ..
      DOUBLE PRECISION  PARA(*), T(LDT,*), V(*)
C     .. Local Scalars ..
      INTEGER           I
C     .. Local Arrays ..
      DOUBLE PRECISION  TMP(LDTMP,LDTMP)
      INTEGER           IPIV(LDTMP)
C     .. External Subroutines ..
      EXTERNAL          DGESV, DLACPY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, INT, MAX, SQRT
C
C     .. Executable Statements ..
C
      INFO = 0
      IF ( N.EQ.0 )
     $   RETURN
C
C     Gauss with partial pivoting is ok, but better use RQ.
C
      CALL DLACPY( 'All', N, N, T, LDT, TMP, LDTMP )
      V(1) = ONE
      DO 10 I = 2, N
         V(I) = ZERO
   10 CONTINUE
      CALL DGESV( N, 1, TMP, LDTMP, IPIV, V, LDTMP, INFO )
      IF ( INFO.NE.0 ) THEN
         PRINT*, 'INVHSE FAILED!'
         STOP
      END IF
      CALL DLARFG( N, V(1), V(2), 1, TAU )
      V(1) = ONE

C
C     Last line of INVHSE
C
      END
      SUBROUTINE KDHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, 
     $                   A, LDA, B, LDB,
     $                   ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, WORK,
     $                   LWORK, INFO, NB)
      implicit none
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ, JOB
      INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, LWORK, N
      INTEGER            NB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
     $                   B( LDB, * ), BETA( * ), Q( LDQ, * ), WORK( * ),
     $                   Z( LDZ, * )
      DOUBLE PRECISION   D1, D2, D3, D4, D5
*     ..
*
*  Purpose
*  =======
*
*  KDHGEQZ implements a single-/double-shift version of the QZ method for
*  finding the generalized eigenvalues
*
*  w(j)=(ALPHAR(j) + i*ALPHAI(j))/BETAR(j)   of the equation
*
*       det( A - w(i) B ) = 0
*
*  In addition, the pair A,B may be reduced to generalized Schur form:
*  B is upper triangular, and A is block upper triangular, where the
*  diagonal blocks are either 1-by-1 or 2-by-2, the 2-by-2 blocks having
*  complex generalized eigenvalues (see the description of the argument
*  JOB.)
*
*  If JOB='S', then the pair (A,B) is simultaneously reduced to Schur
*  form by applying one orthogonal tranformation (usually called Q) on
*  the left and another (usually called Z) on the right.  The 2-by-2
*  upper-triangular diagonal blocks of B corresponding to 2-by-2 blocks
*  of A will be reduced to positive diagonal matrices.  (I.e.,
*  if A(j+1,j) is non-zero, then B(j+1,j)=B(j,j+1)=0 and B(j,j) and
*  B(j+1,j+1) will be positive.)
*
*  If JOB='E', then at each iteration, the same transformations
*  are computed, but they are only applied to those parts of A and B
*  which are needed to compute ALPHAR, ALPHAI, and BETAR.
*
*  If JOB='S' and COMPQ and COMPZ are 'V' or 'I', then the orthogonal
*  transformations used to reduce (A,B) are accumulated into the arrays
*  Q and Z s.t.:
*
*       Q(in) A(in) Z(in)* = Q(out) A(out) Z(out)*
*       Q(in) B(in) Z(in)* = Q(out) B(out) Z(out)*
*
*  Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix
*       Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),
*       pp. 241--256.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          = 'E': compute only ALPHAR, ALPHAI, and BETA.  A and B will
*                 not necessarily be put into generalized Schur form.
*          = 'S': put A and B into generalized Schur form, as well
*                 as computing ALPHAR, ALPHAI, and BETA.
*
*  COMPQ   (input) CHARACTER*1
*          = 'N': do not modify Q.
*          = 'V': multiply the array Q on the right by the transpose of
*                 the orthogonal tranformation that is applied to the
*                 left side of A and B to reduce them to Schur form.
*          = 'I': like COMPQ='V', except that Q will be initialized to
*                 the identity first.
*
*  COMPZ   (input) CHARACTER*1
*          = 'N': do not modify Z.
*          = 'V': multiply the array Z on the right by the orthogonal
*                 tranformation that is applied to the right side of
*                 A and B to reduce them to Schur form.
*          = 'I': like COMPZ='V', except that Z will be initialized to
*                 the identity first.
*
*  N       (input) INTEGER
*          The order of the matrices A, B, Q, and Z.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          It is assumed that A is already upper triangular in rows and
*          columns 1:ILO-1 and IHI+1:N.
*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
*          On entry, the N-by-N upper Hessenberg matrix A.  Elements
*          below the subdiagonal must be zero.
*          If JOB='S', then on exit A and B will have been
*             simultaneously reduced to generalized Schur form.
*          If JOB='E', then on exit A will have been destroyed.
*             The diagonal blocks will be correct, but the off-diagonal
*             portion will be meaningless.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max( 1, N ).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
*          On entry, the N-by-N upper triangular matrix B.  Elements
*          below the diagonal must be zero.  2-by-2 blocks in B
*          corresponding to 2-by-2 blocks in A will be reduced to
*          positive diagonal form.  (I.e., if A(j+1,j) is non-zero,
*          then B(j+1,j)=B(j,j+1)=0 and B(j,j) and B(j+1,j+1) will be
*          positive.)
*          If JOB='S', then on exit A and B will have been
*             simultaneously reduced to Schur form.
*          If JOB='E', then on exit B will have been destroyed.
*             Elements corresponding to diagonal blocks of A will be
*             correct, but the off-diagonal portion will be meaningless.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max( 1, N ).
*
*  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
*          ALPHAR(1:N) will be set to real parts of the diagonal
*          elements of A that would result from reducing A and B to
*          Schur form and then further reducing them both to triangular
*          form using unitary transformations s.t. the diagonal of B
*          was non-negative real.  Thus, if A(j,j) is in a 1-by-1 block
*          (i.e., A(j+1,j)=A(j,j+1)=0), then ALPHAR(j)=A(j,j).
*          Note that the (real or complex) values
*          (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the
*          generalized eigenvalues of the matrix pencil A - wB.
*
*  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
*          ALPHAI(1:N) will be set to imaginary parts of the diagonal
*          elements of A that would result from reducing A and B to
*          Schur form and then further reducing them both to triangular
*          form using unitary transformations s.t. the diagonal of B
*          was non-negative real.  Thus, if A(j,j) is in a 1-by-1 block
*          (i.e., A(j+1,j)=A(j,j+1)=0), then ALPHAR(j)=0.
*          Note that the (real or complex) values
*          (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the
*          generalized eigenvalues of the matrix pencil A - wB.
*
*  BETA    (output) DOUBLE PRECISION array, dimension (N)
*          BETA(1:N) will be set to the (real) diagonal elements of B
*          that would result from reducing A and B to Schur form and
*          then further reducing them both to triangular form using
*          unitary transformations s.t. the diagonal of B was
*          non-negative real.  Thus, if A(j,j) is in a 1-by-1 block
*          (i.e., A(j+1,j)=A(j,j+1)=0), then BETA(j)=B(j,j).
*          Note that the (real or complex) values
*          (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the
*          generalized eigenvalues of the matrix pencil A - wB.
*          (Note that BETA(1:N) will always be non-negative, and no
*          BETAI is necessary.)
*
*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ, N)
*          If COMPQ='N', then Q will not be referenced.
*          If COMPQ='V' or 'I', then the transpose of the orthogonal
*             transformations which are applied to A and B on the left
*             will be applied to the array Q on the right.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  LDQ >= 1.
*          If COMPQ='V' or 'I', then LDQ >= N.
*
*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)
*          If COMPZ='N', then Z will not be referenced.
*          If COMPZ='V' or 'I', then the orthogonal transformations
*             which are applied to A and B on the right will be applied
*             to the array Z on the right.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1.
*          If COMPZ='V' or 'I', then LDZ >= N.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= 4*N.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          = 1,...,N: the QZ iteration did not converge.  (A,B) is not
*                     in Schur form, but ALPHAR(i), ALPHAI(i), and
*                     BETA(i), i=INFO+1,...,N should be correct.
*          = N+1,...,2*N: the shift calculation failed.  (A,B) is not
*                     in Schur form, but ALPHAR(i), ALPHAI(i), and
*                     BETA(i), i=INFO-N+1,...,N should be correct.
*          > 2*N:     various "impossible" errors.
*
*  Further Details
*  ===============
*
*  Iteration counters:
*
*  JITER  -- counts iterations.
*  IITER  -- counts iterations run since ILAST was last
*            changed.  This is therefore reset only when a 1-by-1 or
*            2-by-2 block deflates off the bottom.
*
*  =====================================================================
*
*     .. Parameters ..
*    $                     SAFETY = 1.0E+0 )
      DOUBLE PRECISION   HALF, ZERO, ONE, SAFETY
      DOUBLE PRECISION   TIME1, TM2, TIME2
      PARAMETER          ( HALF = 0.5D+0, ZERO = 0.0D+0, ONE = 1.0D+0,
     $                   SAFETY = 1.0D+2 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILAZR2, ILAZRO, ILPIVT, ILQ, ILSCHR, ILZ
      LOGICAL            NDONE
      INTEGER            ICOMPQ, ICOMPZ, IFIRST, IFRSTM, IITER, ILAST,
     $                   ILASTM, IN, ISCHUR, ISTART, J, JC, JCH, JITER,
     $                   JR, MAXIT, DT, NNB, DT2, RT
      DOUBLE PRECISION   A11, A12, A1I, A1R, A21, A22, A2I, A2R, AD11,
     $                   AD11L, AD12, AD12L, AD21, AD21L, AD22, AD22L,
     $                   AD32L, AN, ANORM, ASCALE, ATOL, B11, B1A, B1I,
     $                   B1R, B22, B2A, B2I, B2R, BN, BNORM, BSCALE,
     $                   BTOL, C, C11I, C11R, C12, C21, C22I, C22R, CL,
     $                   CQ, CR, CZ, ESHIFT, S, S1, S1INV, S2, SAFMAX,
     $                   SAFMIN, SCALE, SL, SQI, SQR, SR, SZI, SZR, T,
     $                   TAU, TEMP, TEMP2, TEMPI, TEMPR, U1, U12, U12L,
     $                   U2, ULP, VS, W11, W12, W21, W22, WABS, WI, WR,
     $                   WR2
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   V( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANHS, DLAPY2, DLAPY3, DSECND
      EXTERNAL           LSAME, DLAMCH, DLANHS, DLAPY2, DLAPY3, DSECND
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAG2, DLARFG, DLARTG, DLASET, DLASV2, DROT,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
      INTEGER DQZ, NCOLS, JJ, NREF, NBL
*     ..
*     .. Executable Statements ..
*      NB = ILAENV( 1, '   ', ' ', N, N, -1, -1 )
*      write(*,*)'NB =',NB
*
*     Decode JOB, COMPQ, COMPZ
*
      RT = 3*N+1
      DT = 0
      DT2 = 0
      TIME1 = 0.d0
      TIME2 = 0.d0
      IF( LSAME( JOB, 'E' ) ) THEN
         ILSCHR = .FALSE.
         ISCHUR = 1
      ELSE IF( LSAME( JOB, 'S' ) ) THEN
         ILSCHR = .TRUE.
         ISCHUR = 2
      ELSE
         ISCHUR = 0
      END IF
*
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'V' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      ELSE
         ICOMPQ = 0
      END IF
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      ELSE
         ICOMPZ = 0
      END IF
*
*     Check Argument Values
*
      INFO = 0
      IF( ISCHUR.EQ.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPQ.EQ.0 ) THEN
         INFO = -2
      ELSE IF( ICOMPZ.EQ.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -5
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.N ) THEN
         INFO = -8
      ELSE IF( LDB.LT.N ) THEN
         INFO = -10
      ELSE IF( LDQ.LT.1 .OR. ( ILQ .AND. LDQ.LT.N ) ) THEN
         INFO = -15
      ELSE IF( LDZ.LT.1 .OR. ( ILZ .AND. LDZ.LT.N ) ) THEN
         INFO = -17
      ELSE IF( LWORK.LT.MAX( 1, 4*N ) ) THEN
         INFO = -19
      ELSE IF (NB .LT. 5) THEN 
         INFO = -20
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DHGEQZ', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*



      IF( N.LE.0 ) THEN
         WORK( 1 ) = DBLE( 1 )
         RETURN
      END IF
*
*     Initialize Q and Z
*
      IF( ICOMPQ.EQ.3 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )
      IF( ICOMPZ.EQ.3 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
*
*     Machine Constants
*
      IN = IHI + 1 - ILO
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      ULP = DLAMCH( 'E' )*DLAMCH( 'B' )
      ANORM = DLANHS( 'F', IN, A( ILO, ILO ), LDA, WORK )
      BNORM = DLANHS( 'F', IN, B( ILO, ILO ), LDB, WORK )
      ATOL = MAX( SAFMIN, ULP*ANORM )
      BTOL = MAX( SAFMIN, ULP*BNORM )
      ASCALE = ONE / MAX( SAFMIN, ANORM )
      BSCALE = ONE / MAX( SAFMIN, BNORM )
*
*     Set Eigenvalues IHI+1:N
*
      DO 30 J = IHI + 1, N
         IF( B( J, J ).LT.ZERO ) THEN
            IF( ILSCHR ) THEN
               DO 10 JR = 1, J
                  A( JR, J ) = -A( JR, J )
                  B( JR, J ) = -B( JR, J )
   10          CONTINUE
            ELSE
               A( J, J ) = -A( J, J )
               B( J, J ) = -B( J, J )
            END IF
            IF( ILZ ) THEN
               DO 20 JR = 1, N
                  Z( JR, J ) = -Z( JR, J )
   20          CONTINUE
            END IF
         END IF
         ALPHAR( J ) = A( J, J )
         ALPHAI( J ) = ZERO
         BETA( J ) = B( J, J )
   30 CONTINUE
*
*     If IHI < ILO, skip QZ steps
*
      IF( IHI.LT.ILO )
     $   GO TO 380
*
*     MAIN QZ ITERATION LOOP
*
*     Initialize dynamic indices
*
*     Eigenvalues ILAST+1:N have been found.
*        Column operations modify rows IFRSTM:whatever.
*        Row operations modify columns whatever:ILASTM.
*
*     If only eigenvalues are being computed, then
*        IFRSTM is the row of the last splitting row above row ILAST;
*        this is always at least ILO.
*     IITER counts iterations since the last eigenvalue was found,
*        to tell when to use an extraordinary shift.
*     MAXIT is the maximum number of QZ sweeps allowed.
*
      ILAST = IHI
      IF( ILSCHR ) THEN
         IFRSTM = 1
         ILASTM = N
      ELSE
         IFRSTM = ILO
         ILASTM = IHI
      END IF
      IITER = 0
      ESHIFT = ZERO
      MAXIT = 30*( IHI-ILO+1 )
*
      DO 360 JITER = 1, MAXIT
*
*        Split the matrix if possible.
*
*        Two tests:
*           1: A(j,j-1)=0  or  j=ILO
*           2: B(j,j)=0
*
	   
         IF( ILAST.EQ.ILO ) THEN
*
*           Special case: j=ILAST
*
            GO TO 80
         ELSE
            IF( ABS( A( ILAST, ILAST-1 ) ).LE.ATOL ) THEN
               A( ILAST, ILAST-1 ) = ZERO
               GO TO 80
            END IF
         END IF
*
         IF( ABS( B( ILAST, ILAST ) ).LE.BTOL ) THEN
            B( ILAST, ILAST ) = ZERO
            GO TO 70
         END IF
*
*        General case: j<ILAST
*
         DO 60 J = ILAST - 1, ILO, -1
*
*           Test 1: for A(j,j-1)=0 or j=ILO
*
            IF( J.EQ.ILO ) THEN
               ILAZRO = .TRUE.
            ELSE
               IF( ABS( A( J, J-1 ) ).LE.ATOL ) THEN
                  A( J, J-1 ) = ZERO
                  ILAZRO = .TRUE.
               ELSE
                  ILAZRO = .FALSE.
               END IF
            END IF
*
*           Test 2: for B(j,j)=0
*
            IF( ABS( B( J, J ) ).LT.BTOL ) THEN
               B( J, J ) = ZERO
*
*              Test 1a: Check for 2 consecutive small subdiagonals in A
*
               ILAZR2 = .FALSE.
               IF( .NOT.ILAZRO ) THEN
                  TEMP = ABS( A( J, J-1 ) )
                  TEMP2 = ABS( A( J, J ) )
                  TEMPR = MAX( TEMP, TEMP2 )
                  IF( TEMPR.LT.ONE .AND. TEMPR.NE.ZERO ) THEN
                     TEMP = TEMP / TEMPR
                     TEMP2 = TEMP2 / TEMPR
                  END IF
                  IF( TEMP*( ASCALE*ABS( A( J+1, J ) ) ).LE.TEMP2*
     $                ( ASCALE*ATOL ) )ILAZR2 = .TRUE.
               END IF
*
*              If both tests pass (1 & 2), i.e., the leading diagonal
*              element of B in the block is zero, split a 1x1 block off
*              at the top. (I.e., at the J-th row/column) The leading
*              diagonal element of the remainder can also be zero, so
*              this may have to be done repeatedly.
*
               IF( ILAZRO .OR. ILAZR2 ) THEN
                  DO 40 JCH = J, ILAST - 1
                     TEMP = A( JCH, JCH )
                     CALL DLARTG( TEMP, A( JCH+1, JCH ), C, S,
     $                            A( JCH, JCH ) )
                     A( JCH+1, JCH ) = ZERO
                     CALL DROT( ILASTM-JCH, A( JCH, JCH+1 ), LDA,
     $                          A( JCH+1, JCH+1 ), LDA, C, S )
                     CALL DROT( ILASTM-JCH, B( JCH, JCH+1 ), LDB,
     $                          B( JCH+1, JCH+1 ), LDB, C, S )
                     IF( ILQ )
     $                  CALL DROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1,
     $                             C, S )
                     IF( ILAZR2 )
     $                  A( JCH, JCH-1 ) = A( JCH, JCH-1 )*C
                     ILAZR2 = .FALSE.
                     IF( ABS( B( JCH+1, JCH+1 ) ).GE.BTOL ) THEN
                        IF( JCH+1.GE.ILAST ) THEN
                           GO TO 80
                        ELSE
                           IFIRST = JCH + 1
                           GO TO 110
                        END IF
                     END IF
                     B( JCH+1, JCH+1 ) = ZERO
   40             CONTINUE
                  GO TO 70
               ELSE
*
*                 Only test 2 passed -- chase the zero to B(ILAST,ILAST)
*                 Then process as in the case B(ILAST,ILAST)=0
*
                  DO 50 JCH = J, ILAST - 1
                     TEMP = B( JCH, JCH+1 )
                     CALL DLARTG( TEMP, B( JCH+1, JCH+1 ), C, S,
     $                            B( JCH, JCH+1 ) )
                     B( JCH+1, JCH+1 ) = ZERO
                     IF( JCH.LT.ILASTM-1 )
     $                  CALL DROT( ILASTM-JCH-1, B( JCH, JCH+2 ), LDB,
     $                             B( JCH+1, JCH+2 ), LDB, C, S )
                     CALL DROT( ILASTM-JCH+2, A( JCH, JCH-1 ), LDA,
     $                          A( JCH+1, JCH-1 ), LDA, C, S )
                     IF( ILQ )
     $                  CALL DROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1,
     $                             C, S )
                     TEMP = A( JCH+1, JCH )
                     CALL DLARTG( TEMP, A( JCH+1, JCH-1 ), C, S,
     $                            A( JCH+1, JCH ) )
                     A( JCH+1, JCH-1 ) = ZERO
                     CALL DROT( JCH+1-IFRSTM, A( IFRSTM, JCH ), 1,
     $                          A( IFRSTM, JCH-1 ), 1, C, S )
                     CALL DROT( JCH-IFRSTM, B( IFRSTM, JCH ), 1,
     $                          B( IFRSTM, JCH-1 ), 1, C, S )
                     IF( ILZ )
     $                  CALL DROT( N, Z( 1, JCH ), 1, Z( 1, JCH-1 ), 1,
     $                             C, S )
   50             CONTINUE
                  GO TO 70
               END IF
            ELSE IF( ILAZRO ) THEN
*
*              Only test 1 passed -- work on J:ILAST
*
               IFIRST = J
               GO TO 110
            END IF
*
*           Neither test passed -- try next J
*
   60    CONTINUE
*
*        (Drop-through is "impossible")
*
         INFO = N + 1
         GO TO 420
*
*        B(ILAST,ILAST)=0 -- clear A(ILAST,ILAST-1) to split off a
*        1x1 block.
*
   70    CONTINUE
         TEMP = A( ILAST, ILAST )
         CALL DLARTG( TEMP, A( ILAST, ILAST-1 ), C, S,
     $                A( ILAST, ILAST ) )
         A( ILAST, ILAST-1 ) = ZERO
         CALL DROT( ILAST-IFRSTM, A( IFRSTM, ILAST ), 1,
     $              A( IFRSTM, ILAST-1 ), 1, C, S )
         CALL DROT( ILAST-IFRSTM, B( IFRSTM, ILAST ), 1,
     $              B( IFRSTM, ILAST-1 ), 1, C, S )
         IF( ILZ )
     $      CALL DROT( N, Z( 1, ILAST ), 1, Z( 1, ILAST-1 ), 1, C, S )
*
*        A(ILAST,ILAST-1)=0 -- Standardize B, set ALPHAR, ALPHAI,
*                              and BETA
*
   80    CONTINUE
         IF( B( ILAST, ILAST ).LT.ZERO ) THEN
            IF( ILSCHR ) THEN
               DO 90 J = IFRSTM, ILAST
                  A( J, ILAST ) = -A( J, ILAST )
                  B( J, ILAST ) = -B( J, ILAST )
   90          CONTINUE
            ELSE
               A( ILAST, ILAST ) = -A( ILAST, ILAST )
               B( ILAST, ILAST ) = -B( ILAST, ILAST )
            END IF
            IF( ILZ ) THEN
               DO 100 J = 1, N
                  Z( J, ILAST ) = -Z( J, ILAST )
  100          CONTINUE
            END IF
         END IF
         ALPHAR( ILAST ) = A( ILAST, ILAST )
         ALPHAI( ILAST ) = ZERO
         BETA( ILAST ) = B( ILAST, ILAST )
*
*        Go to next block -- exit if finished.
*
         ILAST = ILAST - 1
         IF( ILAST.LT.ILO )
     $      GO TO 380
*
*        Reset counters
*
         IITER = 0
         ESHIFT = ZERO
         IF( .NOT.ILSCHR ) THEN
            ILASTM = ILAST
            IF( IFRSTM.GT.ILAST )
     $         IFRSTM = ILO
         END IF
         GO TO 350
*
*        QZ step
*
*        This iteration only involves rows/columns IFIRST:ILAST. We
*        assume IFIRST < ILAST, and that the diagonal of B is non-zero.
*
  110    CONTINUE
         IITER = IITER + 1
         DT = DT + 1
         IF( .NOT.ILSCHR ) THEN
            IFRSTM = IFIRST
         END IF
*
*        Compute single shifts.
*
*        At this point, IFIRST < ILAST, and the diagonal elements of
*        B(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
*        magnitude)
*
         IF( ( IITER / 10 )*10.EQ.IITER ) THEN
*
*           Exceptional shift.  Chosen for no particularly good reason.
*           (Single shift only.)
*
            IF( ( DBLE( MAXIT )*SAFMIN )*ABS( A( ILAST-1, ILAST ) ).LT.
     $          ABS( B( ILAST-1, ILAST-1 ) ) ) THEN
               ESHIFT = ESHIFT + A( ILAST-1, ILAST ) /
     $                  B( ILAST-1, ILAST-1 )
            ELSE
               ESHIFT = ESHIFT + ONE / ( SAFMIN*DBLE( MAXIT ) )
            END IF
            S1 = ONE
            WR = ESHIFT
*
         ELSE
*
*           Shifts based on the generalized eigenvalues of the
*           bottom-right 2x2 block of A and B. The first eigenvalue
*           returned by DLAG2 is the Wilkinson shift (AEP p.512),
*
            CALL DLAG2( A( ILAST-1, ILAST-1 ), LDA,
     $                  B( ILAST-1, ILAST-1 ), LDB, SAFMIN*SAFETY, S1,
     $                  S2, WR, WR2, WI )
*
            TEMP = MAX( S1, SAFMIN*MAX( ONE, ABS( WR ), ABS( WI ) ) )
            IF( WI.NE.ZERO )
     $         GO TO 200
         END IF
*
*        Fiddle with shift to avoid overflow
*
         TEMP = MIN( ASCALE, ONE )*( HALF*SAFMAX )
         IF( S1.GT.TEMP ) THEN
            SCALE = TEMP / S1
         ELSE
            SCALE = ONE
         END IF
*
         TEMP = MIN( BSCALE, ONE )*( HALF*SAFMAX )
         IF( ABS( WR ).GT.TEMP )
     $      SCALE = MIN( SCALE, TEMP / ABS( WR ) )
         S1 = SCALE*S1
         WR = SCALE*WR
*
*        Now check for two consecutive small subdiagonals.
*
         DO 120 J = ILAST - 1, IFIRST + 1, -1
            ISTART = J
            TEMP = ABS( S1*A( J, J-1 ) )
            TEMP2 = ABS( S1*A( J, J )-WR*B( J, J ) )
            TEMPR = MAX( TEMP, TEMP2 )
            IF( TEMPR.LT.ONE .AND. TEMPR.NE.ZERO ) THEN
               TEMP = TEMP / TEMPR
               TEMP2 = TEMP2 / TEMPR
            END IF
            IF( ABS( ( ASCALE*A( J+1, J ) )*TEMP ).LE.( ASCALE*ATOL )*
     $          TEMP2 )GO TO 130
  120    CONTINUE
*
         ISTART = IFIRST
  130    CONTINUE
*
*        Do an implicit single-shift QZ sweep.
*
*        Initial Q
*
         TEMP = S1*A( ISTART, ISTART ) - WR*B( ISTART, ISTART )
         TEMP2 = S1*A( ISTART+1, ISTART )
         CALL DLARTG( TEMP, TEMP2, C, S, TEMPR )
*
*        Sweep
*
         TM2 = DSECND()
         J = ISTART
*
         DO 140 JC = J, ILASTM
            TEMP = C*A( J, JC ) + S*A( J+1, JC )
            A( J+1, JC ) = -S*A( J, JC ) + C*A( J+1, JC )
            A( J, JC ) = TEMP
            TEMP2 = C*B( J, JC ) + S*B( J+1, JC )
            B( J+1, JC ) = -S*B( J, JC ) + C*B( J+1, JC )
            B( J, JC ) = TEMP2
  140    CONTINUE
         IF( ILQ ) THEN
             DO 150 JR = 1, N
                TEMP = C*Q( JR, J ) + S*Q( JR, J+1 )
                Q( JR, J+1 ) = -S*Q( JR, J ) + C*Q( JR, J+1 )
                Q( JR, J ) = TEMP
  150        CONTINUE
         ENDIF
*
         TEMP = B( J+1, J+1 )
         CALL DLARTG( TEMP, B( J+1, J ), C, S, B( J+1, J+1 ) )
         B( J+1, J ) = ZERO
*
         DO 160 JR = IFRSTM, MIN( J+2, ILAST )
            TEMP = C*A( JR, J+1 ) + S*A( JR, J )
            A( JR, J ) = -S*A( JR, J+1 ) + C*A( JR, J )
            A( JR, J+1 ) = TEMP
  160    CONTINUE
         DO 170 JR = IFRSTM, J
            TEMP = C*B( JR, J+1 ) + S*B( JR, J )
            B( JR, J ) = -S*B( JR, J+1 ) + C*B( JR, J )
            B( JR, J+1 ) = TEMP
  170    CONTINUE
         IF( ILZ ) THEN
            DO 180 JR = 1, N
               TEMP = C*Z( JR, J+1 ) + S*Z( JR, J )
               Z( JR, J ) = -S*Z( JR, J+1 ) + C*Z( JR, J )
               Z( JR, J+1 ) = TEMP
  180       CONTINUE
         END IF
*
*        New Blocked SIngle Code Starts Here
*
         JJ = IFIRST
         NDONE = .TRUE.
         NREF = 1
  190    CONTINUE
            NNB = MIN(NB, ILAST-JJ+1)
            NCOLS = 2
            NBL = NNB - 3
*
*           Move bulge and accumulate transformations Q and Z
*
*           Last block?
            IF ((JJ + NNB -1) .EQ. ILAST) THEN
                NCOLS = 1
                NDONE = .FALSE.
                NBL = ILASTM - (JJ + 2)
            ENDIF
*
*           Transformations from the left
*
            IF (JJ-IFIRST .GT. 0) THEN
                CALL MYROT('Rows', A(1,JJ+3),
     $              LDA, IFIRST+1, 1, 1, JJ-IFIRST,
     $              1, NBL, WORK)
                CALL MYROT('Rows', B(1,JJ+3),
     $              LDB, IFIRST+1, 1, 1, JJ-IFIRST,
     $              1, NBL, WORK)
            ENDIF
            CALL GIV(NNB, A(JJ, JJ), LDA, B(JJ, JJ), LDB,
     $               WORK(NREF), WORK(RT), NCOLS)
*
*           Transformations from the left (last iteration)
*
            IF (.NOT. NDONE) THEN
                CALL MYROT('Rows', A(JJ,JJ+NNB), LDA,
     $               2, 1, 1, NNB-(NCOLS+1), 1, ILASTM-JJ-NNB+1,
     $               WORK(NREF) )
                CALL MYROT('Rows', B(JJ,JJ+NNB), LDB,
     $               2, 2, 1, NNB-(NCOLS+1), 1, ILASTM-JJ-NNB+1,
     $               WORK(NREF))
            ENDIF
 
            IF (ILQ) THEN
               CALL MYROT('X', Q(1,JJ), LDQ, 
     $                 1, 2, 1, NNB-(NCOLS+1),
     $                 1, N, WORK(NREF) )
            ENDIF
*
*           Transformations from the right
*
 
            IF (JJ-IFRSTM .GT. 0) THEN
               CALL MYROT('Cols', A(IFRSTM,JJ), LDA,
     $             1, 2, 1, NNB-(NCOLS+1), 1, JJ-IFRSTM,
     $             WORK(RT))
               CALL MYROT('Cols', B(IFRSTM,JJ), LDB,
     $             1, 2, 1, NNB-(NCOLS+1), 1, JJ-IFRSTM,
     $             WORK(RT))
            ENDIF
 
            IF (ILZ) THEN
               CALL MYROT('Cols', Z(1,JJ), LDZ, 1, 2, 1, 
     $                NNB-(NCOLS+1), 1, N, WORK(RT))
            ENDIF
            NREF = NREF + (NNB-(NCOLS+1)) * 2
            JJ = JJ + NB - 3
         IF (NDONE)  GOTO 190
*
*        New Blocked Code Ends Here
*
         TIME2 = TIME2 +  DSECND() - TM2
*
         GO TO 350
*
*        Use Francis double-shift
*
*        Note: the Francis double-shift should work with real shifts,
*              but only if the block is at least 3x3.
*              This code may break if this point is reached with
*              a 2x2 block with real eigenvalues.
*
  200    CONTINUE
         IF( IFIRST+1.EQ.ILAST ) THEN
*
*           Special case -- 2x2 block with complex eigenvectors
*
*           Step 1: Standardize, that is, rotate so that
*
*                       ( B11  0  )
*                   B = (         )  with B11 non-negative.
*                       (  0  B22 )
*
            CALL DLASV2( B( ILAST-1, ILAST-1 ), B( ILAST-1, ILAST ),
     $                   B( ILAST, ILAST ), B22, B11, SR, CR, SL, CL )
*
            IF( B11.LT.ZERO ) THEN
               CR = -CR
               SR = -SR
               B11 = -B11
               B22 = -B22
            END IF
*
            CALL DROT( ILASTM+1-IFIRST, A( ILAST-1, ILAST-1 ), LDA,
     $                 A( ILAST, ILAST-1 ), LDA, CL, SL )
            CALL DROT( ILAST+1-IFRSTM, A( IFRSTM, ILAST-1 ), 1,
     $                 A( IFRSTM, ILAST ), 1, CR, SR )
*
            IF( ILAST.LT.ILASTM )
     $         CALL DROT( ILASTM-ILAST, B( ILAST-1, ILAST+1 ), LDB,
     $                    B( ILAST, ILAST+1 ), LDA, CL, SL )
            IF( IFRSTM.LT.ILAST-1 )
     $         CALL DROT( IFIRST-IFRSTM, B( IFRSTM, ILAST-1 ), 1,
     $                    B( IFRSTM, ILAST ), 1, CR, SR )
*
            IF( ILQ )
     $         CALL DROT( N, Q( 1, ILAST-1 ), 1, Q( 1, ILAST ), 1, CL,
     $                    SL )
            IF( ILZ )
     $         CALL DROT( N, Z( 1, ILAST-1 ), 1, Z( 1, ILAST ), 1, CR,
     $                    SR )
*
            B( ILAST-1, ILAST-1 ) = B11
            B( ILAST-1, ILAST ) = ZERO
            B( ILAST, ILAST-1 ) = ZERO
            B( ILAST, ILAST ) = B22
*
*           If B22 is negative, negate column ILAST
*
            IF( B22.LT.ZERO ) THEN
               DO 210 J = IFRSTM, ILAST
                  A( J, ILAST ) = -A( J, ILAST )
                  B( J, ILAST ) = -B( J, ILAST )
  210          CONTINUE
*
               IF( ILZ ) THEN
                  DO 220 J = 1, N
                     Z( J, ILAST ) = -Z( J, ILAST )
  220             CONTINUE
               END IF
            END IF
*
*           Step 2: Compute ALPHAR, ALPHAI, and BETA (see refs.)
*
*           Recompute shift
*
            CALL DLAG2( A( ILAST-1, ILAST-1 ), LDA,
     $                  B( ILAST-1, ILAST-1 ), LDB, SAFMIN*SAFETY, S1,
     $                  TEMP, WR, TEMP2, WI )
*
*           If standardization has perturbed the shift onto real line,
*           do another (real single-shift) QR step.
*
            IF( WI.EQ.ZERO )
     $         GO TO 350
            S1INV = ONE / S1
*
*           Do EISPACK (QZVAL) computation of alpha and beta
*
            A11 = A( ILAST-1, ILAST-1 )
            A21 = A( ILAST, ILAST-1 )
            A12 = A( ILAST-1, ILAST )
            A22 = A( ILAST, ILAST )
*
*           Compute complex Givens rotation on right
*           (Assume some element of C = (sA - wB) > unfl )
*                            __
*           (sA - wB) ( CZ   -SZ )
*                     ( SZ    CZ )
*
            C11R = S1*A11 - WR*B11
            C11I = -WI*B11
            C12 = S1*A12
            C21 = S1*A21
            C22R = S1*A22 - WR*B22
            C22I = -WI*B22
*
            IF( ABS( C11R )+ABS( C11I )+ABS( C12 ).GT.ABS( C21 )+
     $          ABS( C22R )+ABS( C22I ) ) THEN
               T = DLAPY3( C12, C11R, C11I )
               CZ = C12 / T
               SZR = -C11R / T
               SZI = -C11I / T
            ELSE
               CZ = DLAPY2( C22R, C22I )
               IF( CZ.LE.SAFMIN ) THEN
                  CZ = ZERO
                  SZR = ONE
                  SZI = ZERO
               ELSE
                  TEMPR = C22R / CZ
                  TEMPI = C22I / CZ
                  T = DLAPY2( CZ, C21 )
                  CZ = CZ / T
                  SZR = -C21*TEMPR / T
                  SZI = C21*TEMPI / T
               END IF
            END IF
*
*           Compute Givens rotation on left
*
*           (  CQ   SQ )
*           (  __      )  A or B
*           ( -SQ   CQ )
*
            AN = ABS( A11 ) + ABS( A12 ) + ABS( A21 ) + ABS( A22 )
            BN = ABS( B11 ) + ABS( B22 )
            WABS = ABS( WR ) + ABS( WI )
            IF( S1*AN.GT.WABS*BN ) THEN
               CQ = CZ*B11
               SQR = SZR*B22
               SQI = -SZI*B22
            ELSE
               A1R = CZ*A11 + SZR*A12
               A1I = SZI*A12
               A2R = CZ*A21 + SZR*A22
               A2I = SZI*A22
               CQ = DLAPY2( A1R, A1I )
               IF( CQ.LE.SAFMIN ) THEN
                  CQ = ZERO
                  SQR = ONE
                  SQI = ZERO
               ELSE
                  TEMPR = A1R / CQ
                  TEMPI = A1I / CQ
                  SQR = TEMPR*A2R + TEMPI*A2I
                  SQI = TEMPI*A2R - TEMPR*A2I
               END IF
            END IF
            T = DLAPY3( CQ, SQR, SQI )
            CQ = CQ / T
            SQR = SQR / T
            SQI = SQI / T
*
*           Compute diagonal elements of QBZ
*
            TEMPR = SQR*SZR - SQI*SZI
            TEMPI = SQR*SZI + SQI*SZR
            B1R = CQ*CZ*B11 + TEMPR*B22
            B1I = TEMPI*B22
            B1A = DLAPY2( B1R, B1I )
            B2R = CQ*CZ*B22 + TEMPR*B11
            B2I = -TEMPI*B11
            B2A = DLAPY2( B2R, B2I )
*
*           Normalize so beta > 0, and Im( alpha1 ) > 0
*
            BETA( ILAST-1 ) = B1A
            BETA( ILAST ) = B2A
            ALPHAR( ILAST-1 ) = ( WR*B1A )*S1INV
            ALPHAI( ILAST-1 ) = ( WI*B1A )*S1INV
            ALPHAR( ILAST ) = ( WR*B2A )*S1INV
            ALPHAI( ILAST ) = -( WI*B2A )*S1INV
*
*           Step 3: Go to next block -- exit if finished.
*
            ILAST = IFIRST - 1
            IF( ILAST.LT.ILO )
     $         GO TO 380
*
*           Reset counters
*
            IITER = 0
            ESHIFT = ZERO
            IF( .NOT.ILSCHR ) THEN
               ILASTM = ILAST
               IF( IFRSTM.GT.ILAST )
     $            IFRSTM = ILO
            END IF
            GO TO 350
         ELSE
            TM2 = DSECND()
            DT2 = DT2 + 1
*
*           Usual case: 3x3 or larger block, using Francis implicit
*                       double-shift
*
*                                    2
*           Eigenvalue equation is  w  - c w + d = 0,
*
*                                         -1 2        -1
*           so compute 1st column of  (A B  )  - c A B   + d
*           using the formula in QZIT (from EISPACK)
*
*           We assume that the block is at least 3x3
*
            AD11 = ( ASCALE*A( ILAST-1, ILAST-1 ) ) /
     $             ( BSCALE*B( ILAST-1, ILAST-1 ) )
            AD21 = ( ASCALE*A( ILAST, ILAST-1 ) ) /
     $             ( BSCALE*B( ILAST-1, ILAST-1 ) )
            AD12 = ( ASCALE*A( ILAST-1, ILAST ) ) /
     $             ( BSCALE*B( ILAST, ILAST ) )
            AD22 = ( ASCALE*A( ILAST, ILAST ) ) /
     $             ( BSCALE*B( ILAST, ILAST ) )
            U12 = B( ILAST-1, ILAST ) / B( ILAST, ILAST )
            AD11L = ( ASCALE*A( IFIRST, IFIRST ) ) /
     $              ( BSCALE*B( IFIRST, IFIRST ) )
            AD21L = ( ASCALE*A( IFIRST+1, IFIRST ) ) /
     $              ( BSCALE*B( IFIRST, IFIRST ) )
            AD12L = ( ASCALE*A( IFIRST, IFIRST+1 ) ) /
     $              ( BSCALE*B( IFIRST+1, IFIRST+1 ) )
            AD22L = ( ASCALE*A( IFIRST+1, IFIRST+1 ) ) /
     $              ( BSCALE*B( IFIRST+1, IFIRST+1 ) )
            AD32L = ( ASCALE*A( IFIRST+2, IFIRST+1 ) ) /
     $              ( BSCALE*B( IFIRST+1, IFIRST+1 ) )
            U12L = B( IFIRST, IFIRST+1 ) / B( IFIRST+1, IFIRST+1 )
*
            V( 1 ) = ( AD11-AD11L )*( AD22-AD11L ) - AD12*AD21 +
     $               AD21*U12*AD11L + ( AD12L-AD11L*U12L )*AD21L
            V( 2 ) = ( ( AD22L-AD11L )-AD21L*U12L-( AD11-AD11L )-
     $               ( AD22-AD11L )+AD21*U12 )*AD21L
            V( 3 ) = AD32L*AD21L
*
            ISTART = IFIRST
*
            CALL DLARFG( 3, V( 1 ), V( 2 ), 1, TAU )
            V( 1 ) = ONE
*
*           Sweep XXX
*
            J = ISTART

*
*              All but last elements: use 3x3 Householder transforms.
*
*              Zero (j-1)st column of A
*
               DO 230 JC = J, ILASTM
                  TEMP = TAU*( A( J, JC )+V( 2 )*A( J+1, JC )+V( 3 )*
     $                   A( J+2, JC ) )
                  A( J, JC ) = A( J, JC ) - TEMP
                  A( J+1, JC ) = A( J+1, JC ) - TEMP*V( 2 )
                  A( J+2, JC ) = A( J+2, JC ) - TEMP*V( 3 )
                  TEMP2 = TAU*( B( J, JC )+V( 2 )*B( J+1, JC )+V( 3 )*
     $                    B( J+2, JC ) )
                  B( J, JC ) = B( J, JC ) - TEMP2
                  B( J+1, JC ) = B( J+1, JC ) - TEMP2*V( 2 )
                  B( J+2, JC ) = B( J+2, JC ) - TEMP2*V( 3 )
  230          CONTINUE
               IF( ILQ ) THEN
                  DO 240 JR = 1, N
                     TEMP = TAU*( Q( JR, J )+V( 2 )*Q( JR, J+1 )+V( 3 )*
     $                      Q( JR, J+2 ) )
                     Q( JR, J ) = Q( JR, J ) - TEMP
                     Q( JR, J+1 ) = Q( JR, J+1 ) - TEMP*V( 2 )
                     Q( JR, J+2 ) = Q( JR, J+2 ) - TEMP*V( 3 )
  240             CONTINUE
               END IF
*
*              Zero j-th column of B (see DLAGBC for details)
*
*              Swap rows to pivot
*
               ILPIVT = .FALSE.
               TEMP = MAX( ABS( B( J+1, J+1 ) ), ABS( B( J+1, J+2 ) ) )
               TEMP2 = MAX( ABS( B( J+2, J+1 ) ), ABS( B( J+2, J+2 ) ) )
               IF( MAX( TEMP, TEMP2 ).LT.SAFMIN ) THEN
                  SCALE = ZERO
                  U1 = ONE
                  U2 = ZERO
                  GO TO 250
               ELSE IF( TEMP.GE.TEMP2 ) THEN
                  W11 = B( J+1, J+1 )
                  W21 = B( J+2, J+1 )
                  W12 = B( J+1, J+2 )
                  W22 = B( J+2, J+2 )
                  U1 = B( J+1, J )
                  U2 = B( J+2, J )
               ELSE
                  W21 = B( J+1, J+1 )
                  W11 = B( J+2, J+1 )
                  W22 = B( J+1, J+2 )
                  W12 = B( J+2, J+2 )
                  U2 = B( J+1, J )
                  U1 = B( J+2, J )
               END IF
*
*              Swap columns if nec.
*
               IF( ABS( W12 ).GT.ABS( W11 ) ) THEN
                  ILPIVT = .TRUE.
                  TEMP = W12
                  TEMP2 = W22
                  W12 = W11
                  W22 = W21
                  W11 = TEMP
                  W21 = TEMP2
               END IF
*
*              LU-factor
*
               TEMP = W21 / W11
               U2 = U2 - TEMP*U1
               W22 = W22 - TEMP*W12
               W21 = ZERO
*
*              Compute SCALE
*
               SCALE = ONE
               IF( ABS( W22 ).LT.SAFMIN ) THEN
                  SCALE = ZERO
                  U2 = ONE
                  U1 = -W12 / W11
                  GO TO 250
               END IF
               IF( ABS( W22 ).LT.ABS( U2 ) )
     $            SCALE = ABS( W22 / U2 )
               IF( ABS( W11 ).LT.ABS( U1 ) )
     $            SCALE = MIN( SCALE, ABS( W11 / U1 ) )
*
*              Solve
*
               U2 = ( SCALE*U2 ) / W22
               U1 = ( SCALE*U1-W12*U2 ) / W11
*
  250          CONTINUE
               IF( ILPIVT ) THEN
                  TEMP = U2
                  U2 = U1
                  U1 = TEMP
               END IF
*
*              Compute Householder Vector
*
               T = SQRT( SCALE**2+U1**2+U2**2 )
               TAU = ONE + SCALE / T
               VS = -ONE / ( SCALE+T )
               V( 1 ) = ONE
               V( 2 ) = VS*U1
               V( 3 ) = VS*U2
*
*              Apply transformations from the right.
*
               DO 260 JR = IFRSTM, MIN( J+3, ILAST )
                  TEMP = TAU*( A( JR, J )+V( 2 )*A( JR, J+1 )+V( 3 )*
     $                   A( JR, J+2 ) )
                  A( JR, J ) = A( JR, J ) - TEMP
                  A( JR, J+1 ) = A( JR, J+1 ) - TEMP*V( 2 )
                  A( JR, J+2 ) = A( JR, J+2 ) - TEMP*V( 3 )
  260          CONTINUE
               DO 270 JR = IFRSTM, J + 2
                  TEMP = TAU*( B( JR, J )+V( 2 )*B( JR, J+1 )+V( 3 )*
     $                   B( JR, J+2 ) )
                  B( JR, J ) = B( JR, J ) - TEMP
                  B( JR, J+1 ) = B( JR, J+1 ) - TEMP*V( 2 )
                  B( JR, J+2 ) = B( JR, J+2 ) - TEMP*V( 3 )
  270          CONTINUE
               IF( ILZ ) THEN
                  DO 280 JR = 1, N
                     TEMP = TAU*( Z( JR, J )+V( 2 )*Z( JR, J+1 )+V( 3 )*
     $                      Z( JR, J+2 ) )
                     Z( JR, J ) = Z( JR, J ) - TEMP
                     Z( JR, J+1 ) = Z( JR, J+1 ) - TEMP*V( 2 )
                     Z( JR, J+2 ) = Z( JR, J+2 ) - TEMP*V( 3 )
  280             CONTINUE
               END IF
               B( J+1, J ) = ZERO
               B( J+2, J ) = ZERO
*
*              New Blocked Code Starts Here
*
               JJ = IFIRST
               NDONE = .TRUE.
               NREF = 1
1000           CONTINUE
                  NNB = MIN(NB, ILAST-JJ+1)
                  NCOLS = 3
                  NBL = NNB - 4
*
*                 Move bulge and accumulate transformations Q and Z
*
*
*                 Last block? 
                  IF ((JJ + NNB -1) .EQ. ILAST) THEN
                        NCOLS = 2
                        NDONE = .FALSE.
                        NBL = ILASTM - (JJ + 3) 
                  ENDIF                             
*
*                 Transformations from the left
*
                  IF (JJ-IFIRST .GT. 0) THEN
                     CALL MYDLAREF('Rows', A(1,JJ+4), 
     $                    LDA, .FALSE., Q,
     $                    LDQ, .TRUE., IFIRST+1, 1, 1, JJ-IFIRST, 
     $                    1, NBL, 1, 1, WORK, 
     $                    D1, D2, D3, D4, D5)
                     CALL MYDLAREF('Rows', B(1,JJ+4), 
     $                    LDB, .FALSE., Q,
     $                    LDQ, .TRUE., IFIRST+1, 1, 1, JJ-IFIRST, 
     $                    1, NBL, 1, 1, WORK, 
     $                    D1, D2, D3, D4, D5)
                  ENDIF
                  CALL QZ2(NNB, A(JJ, JJ), LDA, B(JJ, JJ), LDB, 
     $                 WORK(NREF), WORK(RT), NCOLS)
                  IF (.NOT. NDONE) THEN
                     CALL MYDLAREF('Rows', A(JJ,JJ+NNB), LDA, .FALSE.,Q,
     $                    LDQ, .TRUE., 2, 1, 1, NNB-(NCOLS+1),
     $                    1, ILASTM-JJ-NNB+1,
     $                    1, 1, WORK(NREF), D1, D2, D3, D4, D5)
                     CALL MYDLAREF('Rows', B(JJ,JJ+NNB), LDB, .FALSE.,Q,
     $                    LDQ, .TRUE., 2, 2, 1, NNB-(NCOLS+1),
     $                    1, ILASTM-JJ-NNB+1,
     $                    1, 1, WORK(NREF), D1, D2, D3, D4, D5)
                  ENDIF

                  IF (ILQ) THEN
                     CALL MYDLAREF('Cols', Q(1,JJ), LDQ, .FALSE., Q,
     $                    LDQ, .TRUE., 1, 2, 1, NNB-(NCOLS+1), 
     $                    1, N, 1, 1, WORK(NREF), 
     $                    D1, D2, D3, D4, D5)
                  ENDIF
*
*                 Transformations from the right
*

                  IF (JJ-IFRSTM .GT. 0) THEN
                     CALL MYDLAREF('Cols', A(IFRSTM,JJ), LDA, .FALSE.,Q,
     $                    LDQ, .TRUE., 1, 2, 1, NNB-(NCOLS+1), 1, 
     $                    JJ-IFRSTM,
     $                    1, 1, WORK(RT), D1, D2, D3, D4, D5)
                     CALL MYDLAREF('Cols', B(IFRSTM,JJ), LDB, .FALSE.,Q,
     $                    LDQ, .TRUE., 1, 2, 1, NNB-(NCOLS+1), 1, JJ-
     $                    IFRSTM,
     $                    1, 1, WORK(RT), D1, D2, D3, D4, D5)
                  ENDIF

                  IF (ILZ) THEN
                     CALL MYDLAREF('Cols', Z(1,JJ), LDZ, .FALSE., Z,
     $                    LDQ, .TRUE., 1, 2, 1, NNB-(NCOLS+1), 
     $                    1, N,
     $                    1, 1, WORK(RT), D1, D2, D3, D4, D5)
                  ENDIF
                  NREF = NREF + (NNB-(NCOLS+1)) * 3
                  JJ = JJ + NB - 4
               IF (NDONE)  GOTO 1000
*
*           New Blocked Code Ends Here
*
*
*           Last elements: Use Givens rotations
*
*           Rotations from the left
*
            J = ILAST - 1
            TEMP = A( J, J-1 )
            CALL DLARTG( TEMP, A( J+1, J-1 ), C, S, A( J, J-1 ) )
            A( J+1, J-1 ) = ZERO
*
            DO 300 JC = J, ILASTM
               TEMP = C*A( J, JC ) + S*A( J+1, JC )
               A( J+1, JC ) = -S*A( J, JC ) + C*A( J+1, JC )
               A( J, JC ) = TEMP
               TEMP2 = C*B( J, JC ) + S*B( J+1, JC )
               B( J+1, JC ) = -S*B( J, JC ) + C*B( J+1, JC )
               B( J, JC ) = TEMP2
  300       CONTINUE
            IF( ILQ ) THEN
               DO 310 JR = 1, N
                  TEMP = C*Q( JR, J ) + S*Q( JR, J+1 )
                  Q( JR, J+1 ) = -S*Q( JR, J ) + C*Q( JR, J+1 )
                  Q( JR, J ) = TEMP
  310          CONTINUE
            END IF
*
*           Rotations from the right.
*
            TEMP = B( J+1, J+1 )
            CALL DLARTG( TEMP, B( J+1, J ), C, S, B( J+1, J+1 ) )
            B( J+1, J ) = ZERO
*
            DO 320 JR = IFRSTM, ILAST
               TEMP = C*A( JR, J+1 ) + S*A( JR, J )
               A( JR, J ) = -S*A( JR, J+1 ) + C*A( JR, J )
               A( JR, J+1 ) = TEMP
  320       CONTINUE
            DO 330 JR = IFRSTM, ILAST - 1
               TEMP = C*B( JR, J+1 ) + S*B( JR, J )
               B( JR, J ) = -S*B( JR, J+1 ) + C*B( JR, J )
               B( JR, J+1 ) = TEMP
  330       CONTINUE
            IF( ILZ ) THEN
               DO 340 JR = 1, N
                  TEMP = C*Z( JR, J+1 ) + S*Z( JR, J )
                  Z( JR, J ) = -S*Z( JR, J+1 ) + C*Z( JR, J )
                  Z( JR, J+1 ) = TEMP
  340          CONTINUE
            END IF
            TIME1 = TIME1 + DSECND() - TM2 
*
*           End of Double-Shift code
*
         END IF
*
         GO TO 350
*
*        End of iteration loop
*
  350    CONTINUE
  360 CONTINUE
*
*     Drop-through = non-convergence
*
  370 CONTINUE
      INFO = ILAST
      GO TO 420
*
*     Successful completion of all QZ steps
*
  380 CONTINUE
*
*     Set Eigenvalues 1:ILO-1
*
      DO 410 J = 1, ILO - 1
         IF( B( J, J ).LT.ZERO ) THEN
            IF( ILSCHR ) THEN
               DO 390 JR = 1, J
                  A( JR, J ) = -A( JR, J )
                  B( JR, J ) = -B( JR, J )
  390          CONTINUE
            ELSE
               A( J, J ) = -A( J, J )
               B( J, J ) = -B( J, J )
            END IF
            IF( ILZ ) THEN
               DO 400 JR = 1, N
                  Z( JR, J ) = -Z( JR, J )
  400          CONTINUE
            END IF
         END IF
         ALPHAR( J ) = A( J, J )
         ALPHAI( J ) = ZERO
         BETA( J ) = B( J, J )
  410 CONTINUE
*
*     Normal Termination
*
*	WRITE(*,*)'Number of QZ iterations:',JITER
      INFO = 0
	
*
*     Exit (other than argument error) -- return optimal workspace size
*
  420 CONTINUE
      WORK( 1 ) = DBLE( N )
      RETURN
*
*     End of DHGEQZ
*
      END

      SUBROUTINE MYDLAREF( TYPE, A, LDA, WANTZ, Z, LDZ, BLOCK, IROW1,
     $                    ICOL1, ISTART, ISTOP, ITMP1, ITMP2, LILOZ,
     $                    LIHIZ, VECS, V2, V3, T1, T2, T3 )
*
*  -- ScaLAPACK routine (version 1.4 ALPHA) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 17, 1996
*
*     .. Scalar Arguments ..
      LOGICAL            BLOCK, WANTZ
      CHARACTER          TYPE
      INTEGER            ICOL1, IROW1, ISTART, ISTOP, ITMP1, ITMP2, LDA,
     $                   LDZ, LIHIZ, LILOZ
      DOUBLE PRECISION   T1, T2, T3, V2, V3
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), VECS( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  MYDLAREF applies one or several Householder reflectors of size 3
*     to one or two matrices (if column is specified) on either their
*     rows or columns.
*
*  Arguments
*  =========
*
*  TYPE    (global input) CHARACTER*1
*          If 'R': Apply reflectors to the rows of the matrix
*              (apply from left)
*          Otherwise: Apply reflectors to the columns of the matrix
*          Unchanged on exit.
*
*  A       (global input/output) DOUBLE PRECISION array, (LDA,*)
*          On entry, the matrix to receive the reflections.
*          The updated matrix on exit.
*
*  LDA     (local input) INTEGER
*          On entry, the leading dimension of A.  Unchanged on exit.
*
*  WANTZ   (global input) LOGICAL
*          If .TRUE., then apply any column reflections to Z as well.
*          If .FALSE., then do no additional work on Z.
*
*  Z       (global input/output) DOUBLE PRECISION array, (LDZ,*)
*          On entry, the second matrix to receive column reflections.
*          This is changed only if WANTZ is set.
*
*  LDZ     (local input) INTEGER
*          On entry, the leading dimension of Z.  Unchanged on exit.
*
*  BLOCK   (global input) LOGICAL
*          If .TRUE., then apply several reflectors at once and read
*             their data from the VECS array.
*          If .FALSE., apply the single reflector given by V2, V3,
*             T1, T2, and T3.
*
*  IROW1   (local input/output) INTEGER
*          On entry, the local row element of A.
*          Undefined on output.
*
*
*  ICOL1   (local input/output) INTEGER
*          On entry, the local column element of A.
*          Undefined on output.
*
*  ISTART  (global input) INTEGER
*          Specifies the "number" of the first reflector.  This is
*              used as an index into VECS if BLOCK is set.
*              ISTART is ignored if BLOCK is .FALSE..
*
*  ISTOP   (global input) INTEGER
*          Specifies the "number" of the last reflector.  This is
*              used as an index into VECS if BLOCK is set.
*              ISTOP is ignored if BLOCK is .FALSE..
*
*  ITMP1   (local input) INTEGER
*          Starting range into A.  For rows, this is the local
*              first column.  For columns, this is the local first row.
*
*  ITMP2   (local input) INTEGER
*          Ending range into A.  For rows, this is the local last
*              column.  For columns, this is the local last row.
*
*  LILOZ
*  LIHIZ   (local input) INTEGER
*          These serve the same purpose as ITMP1,ITMP2 but for Z
*              when WANTZ is set.
*
*  VECS    (global input) DOUBLE PRECISION array of size 3*N (matrix
*                                                             size)
*          This holds the size 3 reflectors one after another and this
*              is only accessed when BLOCK is .TRUE.
*
*  V2
*  V3
*  T1
*  T2
*  T3      (global input/output) DOUBLE PRECISION
*          This holds information on a single size 3 Householder
*              reflector and is read when BLOCK is .FALSE., and
*              overwritten when BLOCK is .TRUE.
*
*  Implemented by:  G. Henry, November 17, 1996
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            J, K, IR1, IC1
      DOUBLE PRECISION   H11, H22, SUM, T12, T13, T22, T23, T32, T33,
     $                   V22, V23, V32, V33

*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
      IF (( ISTART .GT. ISTOP ) .OR. (ITMP1 .GT. ITMP2 )) THEN
         RETURN
      END IF

      IR1 = IROW1 
      IC1 = ICOL1
      IF(TYPE.EQ.'R' ) THEN
         IF( BLOCK ) THEN
               DO 30 J = ITMP1, ITMP2
            IR1 = IROW1
            DO 40 K = ISTART, ISTOP
               V2 = VECS( ( K-1 )*3+1 )
               V3 = VECS( ( K-1 )*3+2 )
               T1 = VECS( ( K-1 )*3+3 )
               T2 = T1*V2
               T3 = T1*V3
                  SUM = A( IR1, J ) + V2*A( IR1+1, J ) +
     $                  V3*A( IR1+2, J )
                  A( IR1, J ) = A( IR1, J ) - SUM*T1
                  A( IR1+1, J ) = A( IR1+1, J ) - SUM*T2
                  A( IR1+2, J ) = A( IR1+2, J ) - SUM*T3
               IR1 = IR1 + 1
   40       CONTINUE
   30          CONTINUE
         ELSE
            DO 50 J = ITMP1, ITMP2
               SUM = A( IR1, J ) + V2*A( IR1+1, J ) +
     $               V3*A( IROW1+2, J )
               A( IR1, J ) = A( IROW1, J ) - SUM*T1
               A( IR1+1, J ) = A( IR1+1, J ) - SUM*T2
               A( IR1+2, J ) = A( IR1+2, J ) - SUM*T3
   50       CONTINUE
         END IF
      ELSE
*
*        Do column transforms
*
         IF( BLOCK ) THEN
            DO 110 K = ISTART, ISTOP
               V2 = VECS( ( K-1 )*3+1 )
               V3 = VECS( ( K-1 )*3+2 )
               T1 = VECS( ( K-1 )*3+3 )
               T2 = T1*V2
               T3 = T1*V3
               DO 90 J = ITMP1, ITMP2
                  SUM = A( J, IC1 ) + V2*A( J, IC1+1 ) +
     $                  V3*A( J, IC1+2 )
                  A( J, IC1 ) = A( J, IC1 ) - SUM*T1
                  A( J, IC1+1 ) = A( J, IC1+1 ) - SUM*T2
                  A( J, IC1+2 ) = A( J, IC1+2 ) - SUM*T3
   90          CONTINUE
               IF( WANTZ ) THEN
                  DO 100 J = LILOZ, LIHIZ
                     SUM = Z( J, IC1 ) + V2*Z( J, IC1+1 ) +
     $                     V3*Z( J, IC1+2 )
                     Z( J, IC1 ) = Z( J, IC1 ) - SUM*T1
                     Z( J, IC1+1 ) = Z( J, IC1+1 ) - SUM*T2
                     Z( J, IC1+2 ) = Z( J, IC1+2 ) - SUM*T3
  100             CONTINUE
               END IF
               IC1 = IC1 + 1
  110       CONTINUE
         ELSE
            DO 120 J = ITMP1, ITMP2
               SUM = A( J, IC1 ) + V2*A( J, IC1+1 ) +
     $               V3*A( J, IC1+2 )
               A( J, IC1 ) = A( J, IC1 ) - SUM*T1
               A( J, IC1+1 ) = A( J, IC1+1 ) - SUM*T2
               A( J, IC1+2 ) = A( J, IC1+2 ) - SUM*T3
  120       CONTINUE
         END IF
      END IF
      RETURN
*
*     End of MYDLAREF
*
      END
*

      SUBROUTINE MYROT( TYPE, A, LDA, IROW1, ICOL1, ISTART,
     $                   ISTOP, ITMP1, ITMP2, CS)
      implicit none
*     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            ICOL1, IROW1, ISTART, ISTOP, ITMP1, ITMP2, LDA
      DOUBLE PRECISION   IRW
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), CS( * )
*     ..
*
*  Purpose
*  =======
*
*  MYROT applies one or several Givens rotations to a matrix
*     on either their rows or columns.
*
*  Arguments
*  =========
*
*  TYPE    (global input) CHARACTER*1
*          If 'R': Apply reflectors to the rows of the matrix
*              (apply from left)
*          Otherwise: Apply reflectors to the columns of the matrix
*          Unchanged on exit.
*
*  A       (global input/output) DOUBLE PRECISION array, (LDA,*)
*          On entry, the matrix to receive the reflections.
*          The updated matrix on exit.
*
*  LDA     (local input) INTEGER
*          On entry, the leading dimension of A.  Unchanged on exit.
*
*  IROW1   (local input/output) INTEGER
*          On entry, the local row element of A.
*          Undefined on output.
*
*
*  ICOL1   (local input/output) INTEGER
*          On entry, the local column element of A.
*          Undefined on output.
*
*  ISTART  (global input) INTEGER
*          Specifies the "number" of the first cspair.  This is
*              used as an index into VECS.
*
*  ISTOP   (global input) INTEGER
*          Specifies the "number" of the last reflector.
*
*  ITMP1   (local input) INTEGER
*          Starting range into A.  For rows, this is the local
*              first column.  For columns, this is the local first row.
*
*  ITMP2   (local input) INTEGER
*          Ending range into A.  For rows, this is the local last
*              column.  For columns, this is the local last row.
*
*  CS    (global input) DOUBLE PRECISION array of size 2*N (matrix
*                                                             size)
*          This holds the C and S one after another and this
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            J, K, JC, JR, II
      DOUBLE PRECISION   TEMP, TEMP2, C, S
      DOUBLE PRECISION C1, C2, C3, C4, S1, S2, S3, S4
      DOUBLE PRECISION C5, C6, C7, C8, S5, S6, S7, S8
*     ..

*     ..

*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     
      IF (( ISTART .GT. ISTOP ) .OR. (ITMP1 .GT. ITMP2 )) THEN
         RETURN
      END IF
      IF( TYPE.EQ. 'R' )  THEN
            DO 10 JC = ITMP1, ITMP2
         J = IROW1
         DO 20 K = ISTART, ISTOP
            C = CS( ( K-1 )*2+1 )
            S = CS( ( K-1 )*2+2 )
               TEMP = C*A( J, JC ) + S*A( J+1, JC )
               A( J+1, JC ) = -S*A( J, JC ) + C*A( J+1, JC )
               A( J, JC ) = TEMP
            J = J + 1
20       CONTINUE
10          CONTINUE
      ELSEIF( TYPE.EQ.'C' )  THEN
         J = ICOL1
         DO 40 K = ISTART, ISTOP
            C = CS( ( K-1 )*2+1 )
            S = CS( ( K-1 )*2+2 )
            DO 30 JR = ITMP1, ITMP2
               TEMP = C*A( JR, J+1 ) + S*A( JR, J )
               A( JR, J ) = -S*A( JR, J+1 ) + C*A( JR, J )
               A( JR, J+1 ) = TEMP
30          CONTINUE
            J = J + 1
40       CONTINUE
      ELSE
        J = ICOL1
        DO 80 K = ISTART, ISTOP        
            C = CS( ( K-1 )*2+1 )
            S = CS( ( K-1 )*2+2 )
            DO 90 JR = ITMP1, ITMP2
               TEMP = C*A( JR, J ) + S*A( JR, J+1 )
               A( JR, J+1 ) = -S*A( JR, J ) + C*A( JR, J+1 )
               A( JR, J ) = TEMP
90          CONTINUE
            J = J + 1
80      CONTINUE
      END IF
      RETURN
*
*     End of MYROT
*
      END
      SUBROUTINE QZ2( N, A, LDA, B, LDB, VL, VR, NCOL)
	  IMPLICIT NONE
*
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, N, NCOL
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), VL( * ), VR( * )
*     ..
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILPIVT
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   V( 3 )
*     ..
*     .. External Functions ..
 
      DOUBLE PRECISION   DLAMCH, DLANHS, DLAPY2, DLAPY3
	  DOUBLE PRECISION   TIME1, TIME2, TIME3, TM0, TM1, TM2, TM3, TMDUM
      EXTERNAL           DLAMCH, DLANHS, DLAPY2, DLAPY3
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAG2, DLARFG, DLARTG, DLASET, DLASV2, DROT,
     $                   XERBLA
	  INTEGER            J, JC, JR, JMP
      DOUBLE PRECISION   SCALE, T, VS, U2, 
     $                   SAFMIN, SAFMAX,
     $                   TAU, TEMP, TEMP2, U1, W11, W12, W21, W22

*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Machine Constants
*
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
	  JMP = 0
      DO 290 J = 2, N - NCOL 
*
*        Use 3x3 Householder transforms.
*
*        Zero (j-1)st column of A
*
             V( 1 ) = A( J, J-1 )
             V( 2 ) = A( J+1, J-1 )
             V( 3 ) = A( J+2, J-1 )
*
             CALL DLARFG( 3, A( J, J-1 ), V( 2 ), 1, TAU )
             V( 1 ) = ONE
             A( J+1, J-1 ) = ZERO
             A( J+2, J-1 ) = ZERO
*
*            Save Householder transforms.
*
			 VL( JMP + 1) = V( 2 )
			 VL( JMP + 2) = V( 3 )
			 VL( JMP + 3) = TAU 
*
             DO 230 JC = J, N
                TEMP = TAU*( A( J, JC )+V( 2 )*A( J+1, JC )+V( 3 )*
     $                   A( J+2, JC ) )
                A( J, JC ) = A( J, JC ) - TEMP
                A( J+1, JC ) = A( J+1, JC ) - TEMP*V( 2 )
                A( J+2, JC ) = A( J+2, JC ) - TEMP*V( 3 )
                TEMP2 = TAU*( B( J, JC )+V( 2 )*B( J+1, JC )+V( 3 )*
     $                  B( J+2, JC ) )
                B( J, JC ) = B( J, JC ) - TEMP2
                B( J+1, JC ) = B( J+1, JC ) - TEMP2*V( 2 )
                B( J+2, JC ) = B( J+2, JC ) - TEMP2*V( 3 )
  230        CONTINUE
*
*            Zero j-th column of B (see DLAGBC for details)
*
*            Swap rows to pivot
*
             ILPIVT = .FALSE.
             TEMP = MAX( ABS( B( J+1, J+1 ) ), ABS( B( J+1, J+2 ) ) )
             TEMP2 = MAX( ABS( B( J+2, J+1 ) ), ABS( B( J+2, J+2 ) ) )
             IF( MAX( TEMP, TEMP2 ).LT.SAFMIN ) THEN
                  SCALE = ZERO
                  U1 = ONE
                  U2 = ZERO
                  GO TO 250
             ELSE IF( TEMP.GE.TEMP2 ) THEN
                  W11 = B( J+1, J+1 )
                  W21 = B( J+2, J+1 )
                  W12 = B( J+1, J+2 )
                  W22 = B( J+2, J+2 )
                  U1 = B( J+1, J )
                  U2 = B( J+2, J )
             ELSE
                  W21 = B( J+1, J+1 )
                  W11 = B( J+2, J+1 )
                  W22 = B( J+1, J+2 )
                  W12 = B( J+2, J+2 )
                  U2 = B( J+1, J )
                  U1 = B( J+2, J )
             END IF
*
*            Swap columns if nec.
*
             IF( ABS( W12 ).GT.ABS( W11 ) ) THEN
                  ILPIVT = .TRUE.
                  TEMP = W12
                  TEMP2 = W22
                  W12 = W11
                  W22 = W21
                  W11 = TEMP
                  W21 = TEMP2
             END IF
*
*            LU-factor
*
             TEMP = W21 / W11
             U2 = U2 - TEMP*U1
             W22 = W22 - TEMP*W12
             W21 = ZERO
*
*            Compute SCALE
*
             SCALE = ONE
             IF( ABS( W22 ).LT.SAFMIN ) THEN
                SCALE = ZERO
                U2 = ONE
                U1 = -W12 / W11
                GO TO 250
             END IF
             IF( ABS( W22 ).LT.ABS( U2 ) )
     $            SCALE = ABS( W22 / U2 )
             IF( ABS( W11 ).LT.ABS( U1 ) )
     $            SCALE = MIN( SCALE, ABS( W11 / U1 ) )
*
*            Solve
*
              U2 = ( SCALE*U2 ) / W22
              U1 = ( SCALE*U1-W12*U2 ) / W11
*
  250         CONTINUE
              IF( ILPIVT ) THEN
                  TEMP = U2
                  U2 = U1
                  U1 = TEMP
              END IF
*
*             Compute Householder Vector
*
              T = SQRT( SCALE**2+U1**2+U2**2 )
              TAU = ONE + SCALE / T
              VS = -ONE / ( SCALE+T )
              V( 1 ) = ONE
              V( 2 ) = VS*U1
              V( 3 ) = VS*U2
*
*             Save Householder transforms.
*
			  VR( JMP + 1) = V( 2 )
			  VR( JMP + 2) = V( 3 )
			  VR( JMP + 3) = TAU 
*
*             Apply transformations from the right.
*
              DO 260 JR = 1, MIN( J+3, N )
                  TEMP = TAU*( A( JR, J )+V( 2 )*A( JR, J+1 )+V( 3 )*
     $                   A( JR, J+2 ) )
                  A( JR, J ) = A( JR, J ) - TEMP
                  A( JR, J+1 ) = A( JR, J+1 ) - TEMP*V( 2 )
                  A( JR, J+2 ) = A( JR, J+2 ) - TEMP*V( 3 )
  260         CONTINUE
              DO 270 JR = 1, J + 2
                  TEMP = TAU*( B( JR, J )+V( 2 )*B( JR, J+1 )+V( 3 )*
     $                   B( JR, J+2 ) )
                 B( JR, J ) = B( JR, J ) - TEMP
                 B( JR, J+1 ) = B( JR, J+1 ) - TEMP*V( 2 )
                 B( JR, J+2 ) = B( JR, J+2 ) - TEMP*V( 3 )
  270         CONTINUE
               B( J+1, J ) = ZERO
               B( J+2, J ) = ZERO
			   JMP = JMP + 3
  290 CONTINUE
      RETURN
*
*     End of DHGEQZ
*
      END
      SUBROUTINE GIV(N, A, LDA, B, LDB, LCS, RCS, NCOLS)
      IMPLICIT NONE
      INTEGER N, LDA, LDB, NCOLS
      DOUBLE PRECISION A(LDA, *), B(LDB, *), LCS( * ), RCS( * )
*
*     Moves a bulge using Givens Rotations
*
      INTEGER J, JC, JR, I 
      DOUBLE PRECISION   TEMP, TEMP2, ZERO
      PARAMETER          ( ZERO = 0.0D+0 )

      I = 1
      DO 190 J = 2, N - NCOLS
         TEMP = A( J, J-1 )
         CALL DLARTG( TEMP, A( J+1, J-1 ), LCS(I), LCS(I+1), A(J,J-1))
         A( J+1, J-1 ) = ZERO
*
         DO 140 JC = J, N
            TEMP = LCS(I)*A( J, JC ) + LCS(I+1)*A( J+1, JC )
            A( J+1, JC ) = -LCS(I+1)*A( J, JC ) + LCS(I)*A( J+1, JC )
            A( J, JC ) = TEMP
            TEMP2 = LCS(I)*B( J, JC ) + LCS(I+1)*B( J+1, JC )
            B( J+1, JC ) = -LCS(I+1)*B( J, JC ) + LCS(I)*B( J+1, JC )
            B( J, JC ) = TEMP2
  140    CONTINUE
*
         TEMP = B( J+1, J+1 )
         CALL DLARTG( TEMP, B( J+1,J), RCS(I), RCS(I+1), B( J+1, J+1))
         B( J+1, J ) = ZERO
*
         DO 160 JR = 1, MIN( J+2, N )
            TEMP = RCS(I)*A( JR, J+1 ) + RCS(I+1)*A( JR, J )
            A( JR, J ) = -RCS(I+1)*A( JR, J+1 ) + RCS(I)*A( JR, J )
            A( JR, J+1 ) = TEMP
  160    CONTINUE
         DO 170 JR = 1, J
            TEMP = RCS(I)*B( JR, J+1 ) + RCS(I+1)*B( JR, J )
            B( JR, J ) = -RCS(I+1)*B( JR, J+1 ) + RCS(I)*B( JR, J )
            B( JR, J+1 ) = TEMP
  170    CONTINUE
         I = I + 2
  190 CONTINUE
      RETURN
      END


      SUBROUTINE PDGATTER( M, N, A, IA, JA, DESCA, IRPRNT, ICPRNT,
     $                     AA, LD, WORK)
*
*
*     A Dackland Routine: (Modified PDLAPRNT)
*
*     November 20 1996
*
*     .. Scalar Arguments ..
      INTEGER            IA, ICPRNT, IRPRNT, JA, M, N, LD
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
	  DOUBLE PRECISION   A( * ), WORK( * ), AA(LD, *)
*     ..
*
*
*  Purpose
*  =======
*
*  PDGATTER collects the parts of a distributed matrix sub( A )
*  denoting A(IA:IA+M-1,JA:JA+N-1). The local pieces are sent and
*  collected by the process of coordinates (IRPRNT, ICPRNT).
*
*  Notes
*  =====
*
*  A description vector is associated with each 2D block-cyclicly dis-
*  tributed matrix.  This vector stores the information required to
*  establish the mapping between a matrix entry and its corresponding
*  process and memory location.
*
*  In the following comments, the character _ should be read as
*  "of the distributed matrix".  Let A be a generic term for any 2D
*  block cyclicly distributed matrix.  Its description vector is DESCA:
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DT_A   (global) DESCA( DT_ )   The descriptor type.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the distributed
*                                 matrix A.
*  N_A    (global) DESCA( N_ )    The number of columns in the distri-
*                                 buted matrix A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of A.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of A.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the matrix A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of A is distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array storing the local blocks of the
*                                 distributed matrix A.
*                                 LLD_A >= MAX(1,LOCp(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCp( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCq( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCp() and LOCq() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCp( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCq( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*
*  Arguments
*  =========
*
*  M       (global input) INTEGER
*          The number of rows to be operated on i.e the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on i.e the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input) DOUBLE PRECISION pointer into the local memory to a
*          local array of dimension (LLD_A, LOCq(JA+N-1) ) containing
*          the local pieces of the distributed matrix sub( A ).
*
*  IA      (global input) INTEGER
*          A's global row index, which points to the beginning of the
*          submatrix which is to be operated on.
*
*  JA      (global input) INTEGER
*          A's global column index, which points to the beginning of
*          the submatrix which is to be operated on.
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  IRPRNT  (global input) INTEGER
*          The row index of the gattering process.
*
*  ICPRNT  (global input) INTEGER
*          The column index of the gattering process.
*
*  AA      (local Output) The gattered Matrix 
*
*  LD      (local input) The Leading Dimension of AA
*
*  WORK    (local workspace) DOUBLE PRECISION
*          Working array of minimum size equal to MB_A.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            H, I, IACOL, IAROW, IB, ICTXT, ICURCOL,
     $                   ICURROW, II, IIA, IN, J, JB, JJ, JJA, JN, K,
     $                   LDA, MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_BARRIER, BLACS_GRIDINFO, INFOG2L,
     $                   DGERV2D, DGESD2D
*     ..
*     .. External Functions ..
      INTEGER            ICEIL
      EXTERNAL           ICEIL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $              IIA, JJA, IAROW, IACOL )
      ICURROW = IAROW
      ICURCOL = IACOL
      II = IIA
      JJ = JJA
      LDA = DESCA( LLD_ )
*
*     Handle the first block of column separately
*
      JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
      JB = JN-JA+1
      DO 60 H = 0, JB-1
         IN = MIN( ICEIL( IA, DESCA( MB_ ) ) * DESCA( MB_ ), IA+M-1 )
         IB = IN-IA+1
         IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
            IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
               DO 10 K = 0, IB-1
                  AA(IA+K, JA+H) = A( II+K+(JJ+H-1)*LDA )
   10          CONTINUE
            END IF
         ELSE
            IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
               CALL DGESD2D( ICTXT, IB, 1, A( II+(JJ+H-1)*LDA ), LDA,
     $                       IRPRNT, ICPRNT )
            ELSE IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
               CALL DGERV2D( ICTXT, IB, 1, WORK, DESCA( MB_ ),
     $                       ICURROW, ICURCOL )
               DO 20 K = 1, IB
                  AA( IA+K-1, JA+H )=  WORK( K )
   20          CONTINUE
            END IF
         END IF
         IF( MYROW.EQ.ICURROW )
     $      II = II + IB
         ICURROW = MOD( ICURROW+1, NPROW )
         CALL BLACS_BARRIER( ICTXT, 'All' )
*
*        Loop over remaining block of rows
*
         DO 50 I = IN+1, IA+M-1, DESCA( MB_ )
            IB = MIN( DESCA( MB_ ), IA+M-I )
            IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
               IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                  DO 30 K = 0, IB-1
                     AA( I+K, JA+H )= A( II+K+(JJ+H-1)*LDA )
   30             CONTINUE
               END IF
            ELSE
               IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
                  CALL DGESD2D( ICTXT, IB, 1, A( II+(JJ+H-1)*LDA ),
     $                          LDA, IRPRNT, ICPRNT )
               ELSE IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                  CALL DGERV2D( ICTXT, IB, 1, WORK, DESCA( MB_ ),
     $                          ICURROW, ICURCOL )
                  DO 40 K = 1, IB
                      AA( I+K-1, JA+H) = WORK( K )
   40             CONTINUE
               END IF
            END IF
            IF( MYROW.EQ.ICURROW )
     $         II = II + IB
            ICURROW = MOD( ICURROW+1, NPROW )
            CALL BLACS_BARRIER( ICTXT, 'All' )
   50    CONTINUE
*
        II = IIA
        ICURROW = IAROW
   60 CONTINUE
*
      IF( MYCOL.EQ.ICURCOL )
     $   JJ = JJ + JB
      ICURCOL = MOD( ICURCOL+1, NPCOL )
      CALL BLACS_BARRIER( ICTXT, 'All' )
*
*     Loop over remaining column blocks
*
      DO 130 J = JN+1, JA+N-1, DESCA( NB_ )
         JB = MIN(  DESCA( NB_ ), JA+N-J )
         DO 120 H = 0, JB-1
            IN = MIN( ICEIL( IA, DESCA( MB_ ) ) * DESCA( MB_ ), IA+M-1 )
            IB = IN-IA+1
            IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
               IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                  DO 70 K = 0, IB-1
                     AA( IA+K, J+H ) =  A( II+K+(JJ+H-1)*LDA )
   70             CONTINUE
               END IF
            ELSE
               IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
                  CALL DGESD2D( ICTXT, IB, 1, A( II+(JJ+H-1)*LDA ),
     $                          LDA, IRPRNT, ICPRNT )
               ELSE IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                  CALL DGERV2D( ICTXT, IB, 1, WORK, DESCA( MB_ ),
     $                          ICURROW, ICURCOL )
                  DO 80 K = 1, IB
                     AA( IA+K-1, J+H ) =  WORK( K )
   80             CONTINUE
               END IF
            END IF
            IF( MYROW.EQ.ICURROW )
     $         II = II + IB
            ICURROW = MOD( ICURROW+1, NPROW )
            CALL BLACS_BARRIER( ICTXT, 'All' )
*
*           Loop over remaining block of rows
*
            DO 110 I = IN+1, IA+M-1, DESCA( MB_ )
               IB = MIN( DESCA( MB_ ), IA+M-I )
               IF( ICURROW.EQ.IRPRNT .AND. ICURCOL.EQ.ICPRNT ) THEN
                  IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                     DO 90 K = 0, IB-1
                        AA( I+K, J+H) = A( II+K+(JJ+H-1)*LDA )
   90                CONTINUE
                  END IF
               ELSE
                  IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
                     CALL DGESD2D( ICTXT, IB, 1, A( II+(JJ+H-1)*LDA ),
     $                             LDA, IRPRNT, ICPRNT )
                   ELSE IF( MYROW.EQ.IRPRNT .AND. MYCOL.EQ.ICPRNT ) THEN
                     CALL DGERV2D( ICTXT, IB, 1, WORK, DESCA( MB_ ),
     $                             ICURROW, ICURCOL )
                     DO 100 K = 1, IB
                        AA( I+K-1, J+H) = WORK( K )
  100                CONTINUE
                  END IF
               END IF
               IF( MYROW.EQ.ICURROW )
     $            II = II + IB
               ICURROW = MOD( ICURROW+1, NPROW )
               CALL BLACS_BARRIER( ICTXT, 'All' )
  110       CONTINUE
*
            II = IIA
            ICURROW = IAROW
  120    CONTINUE
*
         IF( MYCOL.EQ.ICURCOL )
     $      JJ = JJ + JB
         ICURCOL = MOD( ICURCOL+1, NPCOL )
         CALL BLACS_BARRIER( ICTXT, 'All' )
*
  130 CONTINUE
*
 9999 FORMAT(A,'(',I6,',',I6,')=',D30.18)
*
      RETURN
*
*     End of PDGATTER
*
      END
      SUBROUTINE PDGGHRD(COMPQ, COMPZ, N, ILO, IHI, A, DESCA, 
     $     B, DESCB, Q, DESCQ, Z, DESCZ, WORK, LWORK, INFO)      
      IMPLICIT NONE
*     2001-10-14 Bjorn Adlerborn, interface for step 1 and 2 of the reduction
*     ..
*     .. Scalar Arguments ..
*     ..
      INTEGER            N, INFO, ILO, IHI, LWORK
      CHARACTER          COMPQ, COMPZ
*     ..
*     .. Array Arguments ..
*     ..
      DOUBLE PRECISION   A(*), B(*), WORK(*)
      DOUBLE PRECISION   Q(*), Z(*)
      INTEGER            DESCA(*), DESCB(*)
      INTEGER            DESCQ(*), DESCZ(*)

*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0)
      

      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      
*     ..
*     .. Local Scalars ..
*     ..
      INTEGER             IPCSROW, IPCSCOL, IPTMP, NOBLK
      INTEGER             IPTAU, IPT, IPW
      INTEGER             MYROW, MYCOL, NPROW, NPCOL, ICTXT
      INTEGER             LDA, LDB, LDQ, LDZ
      LOGICAL             ILZ, ILQ
      INTEGER             ICOMPQ, ICOMPZ, MAXNOJ
	INTEGER				NPROCS, IAM, NB
*     ..
*     .. External Subroutines ..
*     ..
      EXTERNAL            PXERBLA, PDLASET, BLACS_GRIDINFO
      EXTERNAL            MPI_WTIME
      DOUBLE PRECISION   MPI_WTIME,T1, T2
*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC           MAX       
*     ..
*     .. External Subroutines
*     ..
      EXTERNAL           LSAME
      LOGICAL            LSAME

*     ..
*     .. Executable Statements ..
*     ..
	NB = DESCA(NB_)
      ICTXT = DESCA(CTXT_)
	CALL BLACS_PINFO( IAM, NPROCS )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      NOBLK = MAX(NPROW, 2)
	IF (LWORK.EQ.-1) THEN
		WORK(1) = 2*N + MAX( 6*N, N*(NB+1) )
		RETURN
	END IF

*     ..
*     .. Decode COMPQ
*     ..
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'V' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      ELSE
         ICOMPQ = 0
      END IF

*     ..
*     .. Decode COMPZ
*     ..
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      ELSE
         ICOMPZ = 0
      END IF
*     ..
*     .. Test the input parameters
*     ..
      INFO = 0
      IF( ICOMPQ.LE.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPZ.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -4
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -5
	ELSE IF (LWORK .LT. 2*N + MAX( 6*N, N*(NB+1) )) THEN
		INFO = - 6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDGGHRD', -INFO )
         RETURN
      END IF
*
*     Initialize Q and Z if desired.
*
      IF( ICOMPQ.EQ.3 )
     $   CALL PDLASET( 'Full', N, N, ZERO, ONE, Q, 1, 1, DESCQ )
      IF( ICOMPZ.EQ.3 )
     $   CALL PDLASET( 'Full', N, N, ZERO, ONE, Z, 1, 1, DESCZ )
*
*     Quick return if possible
*
      IF( N.LE.1 ) 
     $     RETURN

      MAXNOJ = 2
      IPT = 1
      IPTAU = IPT + NB*NB
      IPW = IPTAU + N
      
      IPCSROW = 1
      IPCSCOL = IPCSROW + MAXNOJ*2*N
      IPTMP = IPCSCOL + MAXNOJ*2*N
*     ..
*     .. Execute the steps
*     ..
      T1 = MPI_WTIME()
      CALL STEP1(ILQ, ILZ, NOBLK, N, ILO, IHI, INFO, A, DESCA, 
     $     B, DESCB, Q, DESCQ, Z, DESCZ, WORK(IPTAU), WORK(IPW), 
     $     NB, WORK(IPT), LWORK)

      T2 = MPI_WTIME()

      IF (INFO.NE.0)  THEN
C          WRITE(*,*)'STEP1 returned=',INFO
          RETURN
      END IF 
      IF (MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
C          WRITE(*,*)'HRT Time=',T2-T1

      END IF
      T1 = MPI_WTIME()
      CALL STEP2(ILQ, ILZ, N, ILO, IHI, A, DESCA, B, DESCB, NB,
     $     INFO, WORK(IPCSROW), WORK(IPCSCOL), WORK(IPTMP), 
     $     MAXNOJ, Q, Z, DESCQ, DESCZ, LWORK)
      T2 = MPI_WTIME()
      IF (INFO.NE.0) THEN
C         WRITE(*,*)'STEP2 returned=',INFO
      END IF
      IF (MYROW.EQ.0.AND.MYCOL.EQ.0) THEN
C          WRITE(*,*)'HR Time=',T2-T1
      END IF

 999  RETURN
      END



      SUBROUTINE PDHGEQZ(JOB, COMPQ, COMPZ, N, ILO, IHI, 
     $                   A, DESCA, B, DESCB,
     $                   ALPHAR, ALPHAI, BETA, Q, DESCQ, Z, DESCZ,
     $				   ILOQ, IHIQ, ILOZ, IHIZ,	
     $                   MAXBULGES, WORK,
     $                   LWORK, INFO)
      implicit none
*     ..
*     .. Scalar Arguments ..
*     ..
      CHARACTER          COMPQ, COMPZ, JOB
      INTEGER            IHI, ILO, INFO, LWORK, N
      INTEGER            MAXBULGES
	INTEGER			   ILOQ, IHIQ, ILOZ, IHIZ
*     ..
*     .. Array Arguments ..
*     ..
      DOUBLE PRECISION   A(*), ALPHAI(*), ALPHAR(*),
     $                   B(*), BETA(*), Q(*), WORK(*),
     $                   Z(*)
      INTEGER            DESCA(*), DESCB(*), DESCQ(*), DESCZ(*)
*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION   HALF, ZERO, ONE, TWO, THREE, SAFETY
      PARAMETER          ( HALF = 0.5D+0, ZERO = 0.0D+0, ONE = 1.0D+0,
     $                   TWO = 2.0D+0, THREE = 3.0D+0, SAFETY = 1.0D+2 )

      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            ROWWISE, COLUMNWISE, ALL
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2, ALL=3)
      INTEGER	           DOWN, UP, RIGHT, LEFT
      PARAMETER		   ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
*     ..
*     .. Local Scalars ..
*     ..
      DOUBLE PRECISION   TIME1, TIME2
      LOGICAL            ILQ, ILSCHR, ILZ
      INTEGER            ICOMPQ, ICOMPZ, IFIRST, IFRSTM, IITER, ILAST,
     $                   ILASTM, IN, ISCHUR, J, JITER,
     $                   JR, MAXIT, DT, DT2, 
     $                   LDA, LDB, LDZ, LDQ, ICTXT,
     $                   NPROW, NPCOL, MYROW, MYCOL, IACOL, IAROW,
     $                   JJ, JJR, LR, LC,IAROW2, IACOL2,
     $                   I, MAXSHIFTS, NPROCS,IR,IC,
     $                   LDIST, SDIST, NBULGES,DEBUG0,
     $                   GOTOL, FJ, LCOL1, LCOL2

      INTEGER          NSBULG, NSMIN, NSMAX, AWSIZE
      PARAMETER        ( NSBULG = 2, NSMIN = 20, NSMAX = 20,
     $                            AWSIZE = 200 )

      LOGICAL AGGRRQ
      PARAMETER         (AGGRRQ = .TRUE.)
      DOUBLE PRECISION PARA(9)



      DOUBLE PRECISION EXTFAC, EXTADD, NOSWP
      PARAMETER        ( EXTFAC = 1.5D0, EXTADD = -4.0D0,
     $                           NOSWP = 1.5D-1 )

      INTEGER            DEBUG1, DEBUG2, IAM, DEBUG3, DEBUG4, I1, I2,II
	DOUBLE PRECISION   UMAXBULGES 
      DOUBLE PRECISION   ANORM, ASCALE, ATOL,
     $                   BNORM, BSCALE,
     $                   BTOL, C, 
     $                   ESHIFT, S, SAFMAX,
     $                   SAFMIN, 
     $                   TEMP, D1, D2,D3, D4,
     $                   ULP, TST1, SMLNUM
	DOUBLE PRECISION   SCAL
	INTEGER			   PHESS, PTRIU, IPV, IPU, P, AWS, IERR, KBOT
	INTEGER			   PDW

	DOUBLE PRECISION   T1, T2, T3, T4, T5, T6, T7, T8, TMPT1, TMPT2
	DOUBLE PRECISION   CHT1, CHT2, CHT3, CHT4, CHT5





	INTEGER				WORKMIN, LDTMP
*     ..
*     .. Local Arrays ..
*     ..
      DOUBLE PRECISION    CS(2)
	

*     ..
*     .. External Functions ..
*     ..
      LOGICAL             LSAME
      EXTERNAL            LSAME
      DOUBLE PRECISION    PDLAMCH, PDLANHS
      EXTERNAL            PDLAMCH, PDLANHS
      EXTERNAL            INDXG2L, INDXG2P
      INTEGER             INDXG2L, INDXG2P
	INTEGER				NB, NSBULGE
	DOUBLE PRECISION	MPI_WTIME
	EXTERNAL			MPI_WTIME
	LOGICAL				REPEATEARLY	
	LOGICAL				RUNTTQZ
	INTEGER				OLDILAST,OLDAWS
	INTEGER				IPALPHAI, IPALPHAR, IPBETA
	INTEGER				IPTMP, IPTMP2, IPWORK
	

*     ..
*     .. External Subroutines ..
*     ..
      EXTERNAL           DLARTG, PDLASET
      EXTERNAL           PXERBLA, PDLACP3
      EXTERNAL           KDHGEQZ
	EXTERNAL			ICEIL
	INTEGER				ICEIL
*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT, MOD
*     ..
*     .. Executable Statements ..
*     ..

*     .. Might put this as parameters
	RUNTTQZ = .FALSE.
	NB = DESCA(NB_)
	NSBULGE = 3
	INFO = 0
	DEBUG0 = 0
	DEBUG1 = 0
      DEBUG2 = 0
      DEBUG3 = 0
	DEBUG4 = 0


      SDIST = NB-3
      LDIST = NB-3
      SDIST = 0
      LDIST = 0

	T1 = 0
	T2 = 0
	T3 = 0
	T4 = 0
	T5 = 0
	CHT1 = 0
	CHT2 = 0
	CHT3 = 0
	CHT4 = 0
	CHT5 = 0
	MAXSHIFTS = MAXBULGES*2

*	Adjust maximum number of shifts,
*	No more than NB can be calculated	
	IF (MAXSHIFTS.GT.NB) THEN
		MAXSHIFTS = NB
	END IF
	MAXBULGES = MAXSHIFTS / 2

	LDTMP = NB+ MAXBULGES*(SDIST+3) 

	WORKMIN = (LDTMP)*(LDTMP) * 6
	WORKMIN = WORKMIN + LDTMP *3
	WORKMIN = WORKMIN + NB*N + 16
		

	
	IF (LWORK.EQ.-1) THEN
		WORK(1) = WORKMIN	
		RETURN
	END IF


	
      FJ = 0



	


	
 


*     ..
*     .. Decode JOB
*     ..
      DT = 0
      DT2 = 0
      TIME1 = 0.d0
      TIME2 = 0.d0
      IF( LSAME( JOB, 'E' ) ) THEN
         ILSCHR = .FALSE.
         ISCHUR = 1
      ELSE IF( LSAME( JOB, 'S' ) ) THEN
         ILSCHR = .TRUE.
         ISCHUR = 2
      ELSE
         ISCHUR = 0
      END IF

*     ..
*     .. Decode COMPQ
*     ..
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'V' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      ELSE
         ICOMPQ = 0
      END IF

*     ..
*     .. Decode COMPZ
*     ..
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      ELSE
         ICOMPZ = 0
      END IF

      LDA = DESCA(LLD_)
      LDB = DESCB(LLD_)
      IF (ILZ) THEN
         LDZ = DESCZ(LLD_)
      ELSE
         LDZ = 1
      END IF
      IF (ILQ) THEN
         LDQ = DESCQ(LLD_)
      ELSE
         LDQ = 1
      END IF

*
*     Check Argument Values
*
      INFO = 0
      IF( ISCHUR.EQ.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPQ.EQ.0 ) THEN
         INFO = -2
      ELSE IF( ICOMPZ.EQ.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -5
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -6
      ELSE IF( LWORK.LT.MAX( 1, WORKMIN )) THEN
         INFO = -7
      ELSE IF( NB.LT.5 ) THEN
         INFO = -8
	ELSE IF (SDIST.GT.NB-4) THEN
	   INFO = -9	
	ELSE IF (NB.GT.N) THEN
	   INFO = -10	
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA(ICTXT, 'PDHGEQZ', -INFO )
         RETURN
      END IF
*     ..
*     .. Quick return if possible
*     ..
      IF( N.LE.0 ) THEN
         WORK( 1 ) = DBLE( 1 )
         RETURN
      END IF
      ICTXT = DESCA(CTXT_)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )  
	IAM = MYROW + (NPCOL-1)*MYCOL
	
*	.. Adjust number of bulges
*	.. Make sure they fit within a NB sized window
	UMAXBULGES = (DBLE(NB - 2) / DBLE(NSBULGE+3+SDIST))
	IF (MAXBULGES.GT.UMAXBULGES) THEN
		MAXBULGES = UMAXBULGES
		MAXSHIFTS = MAXBULGES * 2
	END IF
	

	IF (DEBUG0.EQ.1.AND.IAM.EQ.0) THEN
		WRITE(*,*)'Number of bulges=',MAXBULGES
	END IF

	IPALPHAI = 1
	IPALPHAR = IPALPHAI + LDTMP 
	IPBETA = IPALPHAR + LDTMP 
	IPTMP = IPBETA + LDTMP
	
	IPTMP2 = IPTMP + (LDTMP)*(LDTMP)

	IPU = IPTMP2 + (LDTMP)*(LDTMP)
	IPV = IPU + (LDTMP)*(LDTMP)

	IPWORK = IPV + (LDTMP)*(LDTMP)
*     ..
*     ..
*     ..

*     ..
*     .. Initialize Q and Z if desired
*     ..
      IF( ICOMPQ.EQ.3 ) THEN
	   CALL PDLASET( 'Full', N, N, ZERO, ONE, Q, 1, 1, DESCQ )
	END IF
      IF( ICOMPZ.EQ.3 ) THEN
	   CALL PDLASET( 'Full', N, N, ZERO, ONE, Z, 1, 1, DESCZ )
	END IF

*     ..
*     .. Machine Constants, calculated by all processors
*     ..
      IN = IHI + 1 - ILO
      SAFMIN = PDLAMCH(ICTXT, 'S')
      SAFMAX = ONE / SAFMIN
*	CALL PDLABAD(ICTXT, SAFMIN, SAFMAX)
      ULP = PDLAMCH(ICTXT, 'P')
      ANORM = PDLANHS( 'F', IN, A, 1, 1, DESCA, WORK )
      BNORM = PDLANHS( 'F', IN, B, 1, 1, DESCB, WORK )
      ATOL = MAX( SAFMIN, ULP*ANORM )
      BTOL = MAX( SAFMIN, ULP*BNORM )
      ASCALE = ONE / MAX( SAFMIN, ANORM )
      BSCALE = ONE / MAX( SAFMIN, BNORM )
	SMLNUM = SAFMIN*( N/ULP )

*     ..
*     .. Set Eigenvalues IHI+1:N
*     ..
      DO 30 J = IHI + 1, N
         IACOL = INDXG2P(J, NB, 0, 0, NPCOL)
         JJ = INDXG2L(J, NB, 0, 0, NPCOL)
         JJR = INDXG2L(J, NB, 0, 0, NPROW)
         IF (MYCOL.EQ.IACOL) THEN
            IF ( B( JJR + (JJ-1)*LDB ).LT.ZERO ) THEN
               IF ( ILSCHR ) THEN
                  DO 10 JR = 1, J
                     IAROW = INDXG2P(JR, NB, 0, 0, NPROW)
                     JJR = INDXG2L(JR, NB, 0, 0, NPROW)
                     IF (MYROW.EQ.IAROW)THEN
                        A( JJR + (JJ-1)*LDA ) = -A( JJR + (JJ-1)*LDA )
                        B( JJR + (JJ-1)*LDB ) = -B( JJR + (JJ-1)*LDB )
                     END IF
 10               CONTINUE
               ELSE
                  IAROW = INDXG2P(J, NB, 0, 0, NPROW)
                  IF (MYROW.EQ.IAROW) THEN
                     A( JJR + (JJ-1)*LDA ) = -A( JJR + (JJ-1)*LDA )
                     B( JJR + (JJ-1)*LDB ) = -B( JJR + (JJ-1)*LDB )
                  ENDIF
               END IF
               IF( ILZ ) THEN
                  DO 20 JR = 1, N
                     IAROW = INDXG2P(JR, NB, 0, 0, NPROW)
                     JJR = INDXG2L(JR, NB, 0, 0, NPROW)
                     IF (MYROW.EQ.IAROW) THEN
                        Z( JJR + (JJ-1)*LDZ ) = -Z( JJR + (JJ-1)*LDZ )
                     ENDIF
 20               CONTINUE
               END IF
            END IF
            IAROW = INDXG2P(J, NB, 0, 0, NPROW)
            IF (MYROW.EQ.IAROW) THEN
               ALPHAR( J ) = A( JJR + (JJ-1)*LDA )
               ALPHAI( J ) = ZERO
               BETA( J ) = B( JJR + (JJ-1)*LDB )
            END IF
         END IF
 30   CONTINUE
*
*     If IHI < ILO, skip QZ steps
*
      IF( IHI.LT.ILO )
     $   GO TO 380
*
*     MAIN QZ ITERATION LOOP
*
*     Initialize dynamic indices
*
*     Eigenvalues ILAST+1:N have been found.
*        Column operations modify rows IFRSTM:whatever.
*        Row operations modify columns whatever:ILASTM.
*
*     If only eigenvalues are being computed, then
*        IFRSTM is the row of the last splitting row above row ILAST;
*        this is always at least ILO.
*     IITER counts iterations since the last eigenvalue was found,
*        to tell when to use an extraordinary shift.
*     MAXIT is the maximum number of QZ sweeps allowed.
*
	PARA(9) = 1
	PARA(8) = 1
      ILAST = IHI
      IF( ILSCHR ) THEN
         IFRSTM = 1
         ILASTM = N
      ELSE
         IFRSTM = ILO
         ILASTM = IHI
      END IF
      IITER = 1
      ESHIFT = ZERO
      
      MAXIT = 30*( IHI-ILO+1 )
	IFIRST = 1

*	MAXIT =	1
      IF ((DEBUG1.EQ.1)) THEN
         WRITE(*,*)'PQZ step0, MAXIT=', MAXIT, 
     $        ' IFRSTM=', IFRSTM, ' ILASTM=', ILASTM,
     $        ' ILO=', ILO, ' IHI=', IHI,
     $		' MYROW=', MYROW, 'MYCOL=',MYCOL,  
     $		' ATOL=', ATOL, 'BTOL=', BTOL
      ENDIF   
      DO 360 JITER = 1, MAXIT
3        CONTINUE   
	   
	   TMPT1 = MPI_WTIME()	
	   AWS = MOD(ILAST,NB)
	   AWS = NB

*	   IF(AWS.LE.50) THEN
*		   IF (NPROW.EQ.1.AND.NPCOL.EQ.1) THEN
*				AWS = MIN(LDTMP,200)
*		   ELSE
*				AWS = NB
*		   END IF
*	   END IF
	   AWS = MIN(AWS, ILAST-IFIRST+1)
	   
	   MAXSHIFTS = MIN(AWS,MAXBULGES * 2)		   
	   IF (AWS.LT.2) GOTO 4


	   P = ILAST-AWS+1
	   
	   IAROW = INDXG2P(P, NB, 0, 0, NPROW)
         IACOL = INDXG2P(P, NB, 0, 0, NPCOL)

	   IAROW2 = INDXG2P(ILAST, NB, 0, 0, NPROW)
         IACOL2 = INDXG2P(ILAST, NB, 0, 0, NPCOL)

         CALL DLASET( 'All', LDTMP,LDTMP, 
     $		ZERO, ONE, WORK(IPU), LDTMP)
         CALL DLASET( 'All', LDTMP, LDTMP, 
     $		ZERO, ONE, WORK(IPV), LDTMP)

	   IF ( P.GT.1 ) THEN
		  CALL PDELGET('A', ' ', SCAL, A, P,P-1,DESCA)
	   ELSE  IF ( P.EQ.1 ) THEN
            SCAL = ZERO
	   END IF
	   	
	   IF (AWS.GT.MOD(ILAST,NB)) THEN
*			CALL PDCOPY(AWS,P,A,DESCA,1,WORK(IPTMP),LDTMP, 
*     $			MYROW,MYCOL, NPROW, NPCOL, ICTXT, IAROW,IACOL)
*			CALL PDCOPY(AWS,P,B,DESCB,1,WORK(IPTMP2),LDTMP, 
*     $			MYROW,MYCOL, NPROW, NPCOL, ICTXT, IAROW,IACOL)


	   CALL PDLACP3(AWS, P, A, DESCA, WORK(IPTMP), 
     $			LDTMP,IAROW,IACOL,0)
	   CALL PDLACP3(AWS, P, B, DESCB, WORK(IPTMP2), 
     $			LDTMP,IAROW,IACOL,0)
     	   END IF
	   IF (MYROW.EQ.IAROW.AND.MYCOL.EQ.IACOL) THEN 
*		 Copy A and B to tmparrays
		 LR = INDXG2L(P,NB, 0,0,NPROW)
		 LC = INDXG2L(P,NB, 0,0,NPCOL)
		 IF (AWS.LE.MOD(ILAST,NB)) THEN
		 CALL DLACPY('All',AWS,AWS,A(LR+(LC-1)*LDA),
     $		LDA,WORK(IPTMP),LDTMP)
		 CALL DLACPY('All',AWS,AWS,B(LR+(LC-1)*LDB),
     $		LDB,WORK(IPTMP2),LDTMP)
		 END IF

	         CALL QZEARLY( AWS, SCAL, PARA, WORK(IPTMP),LDTMP,
     $				 WORK(IPTMP2),LDTMP,
     $                 KBOT, ALPHAR(P), ALPHAI(P), BETA(P), 
     $                 WORK(IPU), LDTMP, 
     $				 WORK(IPV), LDTMP, 
     $				 WORK(IPWORK),
     $                 LWORK-IPWORK+1, IERR )

		  IF (AWS.LE.MOD(ILAST,NB).AND.KBOT.LE.AWS) THEN
	        
*			.. Restore A and B
			CALL DLACPY('All',AWS,AWS,WORK(IPTMP),LDTMP,
     $			A(LR+(LC-1)*LDA),LDA)
			CALL DLACPY('All',AWS,AWS,WORK(IPTMP2),LDTMP,
     $			B(LR+(LC-1)*LDB),LDB)
		  END IF

	      WORK(IPWORK) = DBLE(KBOT)
		  
		  Call DCOPY(AWS,ALPHAR(P),1,WORK(IPWORK+LDTMP),1)
		  Call DCOPY(AWS,ALPHAI(P),1,WORK(IPWORK+LDTMP*2),1)
		  Call DCOPY(AWS,BETA(P),1,WORK(IPWORK+LDTMP*3),1)
		 
	      Call BSend(AWS, 4, ALL, WORK(IPWORK), LDTMP, ICTXT)

	   ELSE
		  Call BRecv(AWS, 4, ALL, WORK(IPWORK), LDTMP, 
     $			IAROW, IACOL, ICTXT)
		  KBOT = INT( WORK(IPWORK) )
		  
		  Call DCOPY(AWS,WORK(IPWORK+LDTMP),1, ALPHAR(P),1)
		  Call DCOPY(AWS,WORK(IPWORK+LDTMP*2),1,ALPHAI(P),1)
		  Call DCOPY(AWS,WORK(IPWORK+LDTMP*3),1,BETA(P),1)

	   END IF

	   IF (KBOT.LT.AWS.AND.AWS.GT.MOD(ILAST,NB)) THEN
*			CALL PDCOPY(AWS,P,A,DESCA,0,WORK(IPTMP),LDTMP, 
*     $			MYROW,MYCOL, NPROW, NPCOL, ICTXT, IAROW,IACOL)
*			CALL PDCOPY(AWS,P,B,DESCB,0,WORK(IPTMP2),LDTMP, 
*     $			MYROW,MYCOL, NPROW, NPCOL, ICTXT, IAROW,IACOL)

		   CALL PDLACP3(AWS, P, A, DESCA, WORK(IPTMP), 
     $			LDTMP,IAROW,IACOL,-1)
		   CALL PDLACP3(AWS, P, B, DESCB, WORK(IPTMP2), 
     $			LDTMP,IAROW,IACOL,-1)
     	   END IF

	



	   IF (DEBUG1.EQ.1.AND.MYROW.EQ.IAROW.AND.MYCOL.EQ.IACOL) THEN	
			PRINT *,'KTTQZ:', AWS-KBOT,' deflations found',JITER
	   END IF	
	   	   
	   IF (IERR.NE.0.AND.IERR.NE.18.AND.
     $			MYROW.EQ.IAROW.AND.MYCOL.EQ.IACOL) THEN
			WRITE(*,*)'QZEARLY FAILED:',IERR
	   ELSEIF(IERR.EQ.18.AND.MYROW.EQ.IAROW.AND.MYCOL.EQ.IACOL) THEN
			INFO = -22
	   END IF


	   IF (KBOT.LT.AWS) THEN
					  	  


		
C		Distribute U (=Q) and V(=Z)
C	    U both directions(needed for Q accumulations), V only columnwize


		IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
			
			IF (ILQ.AND.NPROW.GT.1) THEN
     				Call BSend(AWS, AWS, COLUMNWISE, WORK(IPU),
     $				LDTMP, ICTXT)
			END IF
			IF (NPCOL.GT.1) THEN
				Call BSend(AWS, AWS, ROWWISE, WORK(IPU),
     $				LDTMP, ICTXT)
			END IF

			IF (NPROW.GT.1) THEN
     				Call BSend(AWS,AWS, COLUMNWISE, WORK(IPV),
     $				LDTMP, ICTXT)
			END IF

		ELSE IF (MYROW.EQ.IAROW.AND.NPCOL.GT.1) THEN
			Call BRecv(AWS,AWS, ROWWISE, WORK(IPU),
     $				LDTMP, IAROW,IACOL, ICTXT)
			
		ELSE IF (MYCOL.EQ.IACOL.AND.NPROW.GT.1) THEN
			IF (ILQ) THEN
				Call BRecv(AWS,AWS, COLUMNWISE, WORK(IPU),
     $				LDTMP, IAROW,IACOL, ICTXT)
			END IF
			Call BRecv(AWS,AWS, COLUMNWISE, WORK(IPV),
     $				LDTMP, IAROW,IACOL, ICTXT)
		END IF

		IF (IAROW.NE.IAROW2) THEN
		   DO 3501 IC = 1, NPCOL 
			 IF ((IC-1).EQ.IACOL) GOTO 3501
		     IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.(IC-1))) THEN
				CALL Send(AWS,AWS, DOWN, WORK(IPV), LDTMP,
     $				MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				IF (ILQ) THEN
					CALL Send(AWS,AWS, DOWN, WORK(IPU), LDTMP,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				END IF
			 ELSE IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.(IC-1))) THEN
				CALL Recv(AWS,AWS, UP, WORK(IPV), LDTMP,
     $				MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				IF (ILQ) THEN
					CALL Recv(AWS,AWS, UP, WORK(IPU), LDTMP,
     $					MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				END IF
			 END IF
 3501		   CONTINUE
	    END IF

		IF (IACOL.NE.IACOL2) THEN
		   DO 3601 IR = 1, NPROW 
			 IF ((IR-1).EQ.IAROW) GOTO 3601
		     IF ((MYCOL.EQ.IACOL).AND.(MYROW.EQ.(IR-1))) THEN
				CALL Send(AWS,AWS, RIGHT, WORK(IPU), LDTMP,
     $				MYROW, MYCOL, NPCOL, NPROW, ICTXT)
			 ELSE IF ((MYCOL.EQ.IACOL2).AND.(MYROW.EQ.(IR-1))) THEN
				CALL Recv(AWS,AWS, LEFT, WORK(IPU), LDTMP,
     $				MYROW, MYCOL, NPCOL, NPROW, ICTXT)

			 END IF
 3601		   CONTINUE
	    END IF 

*		  .. End of distribute U and V

		  IF (P.GT.1) THEN
			CALL PDELSET(A, P,P-1, DESCA, SCAL * WORK(IPU))
		  END IF

*		  ..  Apply WORK(PU)=Q and WORK(PV)=Z to A, B, Q, Z with dgemm
		  CALL UPDABQZ(A, B, Q, Z,
     $			LDA, LDB, LDQ, LDZ,
     $			WORK(IPTMP),  WORK(IPTMP2),
     $			LDTMP,	
     $			WORK(IPU), WORK(IPV),
     $			LDTMP,LDTMP,
     $			DESCA, DESCB, DESCQ,DESCZ,
     $			P, ILASTM,
     $			ILO, NB, N, AWS,
     $			ILQ, ILZ,AWS)
			
     		  ILAST = ILAST - (AWS - KBOT )

		  IF (ILAST.LT.ILO) THEN
               GOTO 380
		  END IF
		  IF (ILAST.LT.IFIRST) THEN
			 GOTO 350
		  END IF

		
	   END IF
		 
		
*         DO 60 II = P, P+AWS-1, 2
*            IF ( ALPHAI(II).EQ.ZERO .AND. ALPHAI(II+1).NE.ZERO ) THEN
*               TEMP = ALPHAR(II)
*               ALPHAR(II) = ALPHAR(II+1)
*               ALPHAR(II+1) = ALPHAR(II+2)
*               ALPHAR(II+2) = TEMP
*               TEMP = ALPHAI(II)
*               ALPHAI(II) = ALPHAI(II+1)
*               ALPHAI(II+1) = ALPHAI(II+2)
*               ALPHAI(II+2) = TEMP
*               TEMP = BETA(II)
*               BETA(II) = BETA(II+1)
*               BETA(II+1) = BETA(II+2)
*               BETA(II+2) = TEMP
*            ELSE IF ( ALPHAI(II).NE.ZERO .AND. ALPHAI(II+1).EQ.ZERO )
*     $      THEN
*C              that shouldn't happen
*               PRINT*, 'PROBLEMS WITH EIGENVALUE SORTING'
*            END IF
*   60    CONTINUE



5	   CONTINUE	

	   TMPT2 = MPI_WTIME()

	   T1 = T1 + (TMPT2 - TMPT1)

*	   Try QZ early again if unreducded block is too small	
*	   IF ((ILAST-IFIRST+1) .LE. NB) GOTO 3
			 

4	   CONTINUE

	   TMPT1 = MPI_WTIME()
*
*        Search for converged elements at the bottom of the matrix
*
*        Two tests:
*           1: A(j,j-1)=0  or  j=ILO
*           2: B(j,j)=0
*
         IF ((DEBUG4.EQ.1)) THEN
            WRITE(*,*)'PQZ step1, JITER=', JITER, 
     $           ' IFIRST=', IFIRST, ' ILAST=', ILAST,
     $           ' IFRSTM=', IFRSTM, ' ILASTM=', ILASTM,
     $		' MYROW=', MYROW, 'MYCOL=',MYCOL
         ENDIF  
         IF( ILAST.EQ.ILO ) THEN
*
*           Special case: j=ILAST
*
            GO TO 80
         ELSE
*	      ..
*	      .. Test for A(ILAST, ILAST-1) = 0
*	      ..
*            CALL PDLACP3(2, ILAST-1, A, DESCA, WORK(IPTMP), LDTMP, 
*     $           -1, -1, 0)
*		  IF (ABS( WORK(IPTMP+1)).LT.ATOL) THEN
*			 WORK(IPTMP + 1) = ZERO
*              IACOL = INDXG2P(ILAST-1, NB, 0, 0, NPCOL)
*              IAROW = INDXG2P(ILAST-1, NB, 0, 0, NPROW)               			 
*			CALL PDLACP3(2, ILAST-1, A, DESCA,  WORK(IPTMP), LDTMP, 
*     $			IAROW, IACOL, -1)
*			GOTO 80
*		  END IF

*	      ..
*		  .. Test for B(ILAST, ILAST) = 0
*	      ..
*            CALL PDLACP3(1, ILAST, B, DESCB,  WORK(IPTMP), LDTMP, 
*     $           -1, -1, 0)
*            IF (ABS( WORK(IPTMP)).LE.BTOL) THEN
*               IACOL = INDXG2P(ILAST, NB, 0, 0, NPCOL)
*               IAROW = INDXG2P(ILAST, NB, 0, 0, NPROW)               
*               IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
*                  LC = INDXG2L(ILAST, NB, 0, 0, NPCOL)
*                  LR = INDXG2L(ILAST, NB, 0, 0, NPROW)
*                  B(LR+(LC-1)*LDB) = ZERO
*               ENDIF
*               GOTO 70
*            ENDIF
         END IF
	
		   
	   TMPT2 = MPI_WTIME()
	   
	   T2 = T2 + (TMPT2 - TMPT1)
*	   ..
*        .. Find a suitable start for the iteration
*	   .. This procedure sets a value for IFIRST, but might
*	   .. also have effect on ILAST
*	   .. GOTOL tells us what kind of action we should take
*	   .. If GOTOL = 110 we should start an iteration, if GOTOL=70
*	   .. a diagonal element of B is zero and needs to be dealt with
*	   .. If GOTOL = 80 then a subdiagonal element of A is zero and need to
*	   .. be dealt with
*	   ..
	   TMPT1 = MPI_WTIME()
	   IF (ILAST.LT.IFIRST) THEN		
			CALL FINDSTART(A,
     $			DESCA, LDA, 
     $			IFIRST, ILAST,
     $			ILO, IHI, GOTOL,ATOL)
			GOTO 350
	   ELSE
			CALL FINDSTART(A,
     $			DESCA, LDA, 
     $			IFIRST, ILAST,
     $			ILO,IHI, GOTOL,ATOL)


	   END IF
	   TMPT2 = MPI_WTIME()

	   T3 = T3 + (TMPT2 - TMPT1)

         IF ((DEBUG1.EQ.1)) THEN
            WRITE(*,*)'PQZ step2, JITER=', JITER, 
     $           ' IFIRST=', IFIRST, ' ILAST=', ILAST,
     $           ' IFRSTM=', IFRSTM, ' ILASTM=', ILASTM,
     $           ' GOTOL=', GOTOL,
     $		   ' MYROW=', MYROW, 'MYCOL=',MYCOL
         ENDIF  

         IF (GOTOL.EQ.70) THEN
            GOTO 70
         END IF
         IF (GOTOL.EQ.80) THEN
            GOTO 80
         END IF
         IF (GOTOL.EQ.110) THEN
            GOTO 110
         END IF
*
*        (Drop-through is "impossible")
*
         INFO = N + 1
         GO TO 420
*
*        B(ILAST,ILAST)=0 -- clear A(ILAST,ILAST-1) to split off a
*        1x1 block.
*
   70    CONTINUE
         CALL PDLACP3(2, ILAST-1, A, DESCA,  WORK(IPTMP), 
     $        LDTMP, -1, -1, 0)
                  
*         TEMP = TMPW( 2, 2 )
	   TEMP =  WORK(IPTMP + NB + 1)
         CALL DLARTG( TEMP,  WORK(IPTMP+1), C, S,
     $        WORK(IPTMP) )
         WORK(IPTMP+1) = ZERO
         
         IAROW = INDXG2P(ILAST-1, NB, 0, 0, NPROW)
         IACOL = INDXG2P(ILAST-1, NB, 0, 0, NPCOL)
         
         CALL PDLACP3(2, ILAST-1, A, DESCA,  WORK(IPTMP), LDTMP, 
     $        IAROW, IACOL, -1)


         CS(1) = C
         CS(2) = S

         CALL PDMYROT('X', A, DESCA,
     $        LDA, NB, IFRSTM, ILAST-1, 1, 1,
     $        IFRSTM, ILAST, 
     $        CS,  WORK(IPTMP),LDTMP)

         CALL PDMYROT('X', B, DESCB, 
     $        LDB, NB, IFRSTM, ILAST-1, 1, 1,
     $        IFRSTM, ILAST,
     $        CS,  WORK(IPTMP),LDTMP)

         IF (ILZ) THEN
            CALL PDMYROT('X', Z, DESCZ,
     $           LDZ, NB, IFRSTM, ILAST-1, 1, 1,
     $           1, N,
     $           CS,  WORK(IPTMP),LDTMP)
         END IF

*
*        A(ILAST,ILAST-1)=0 -- Standardize B, set ALPHAR, ALPHAI,
*                              and BETA
*
   80    CONTINUE

         CALL PDLACP3(1, ILAST, B, DESCB,  WORK(IPTMP), LDTMP, 
     $        -1, -1, 0)

         IF ( WORK(IPTMP).LT.ZERO) THEN         

            LC = INDXG2L(ILAST, NB, 0, 0, NPCOL)
            IACOL = INDXG2P(ILAST, NB, 0, 0, NPCOL)
            IF( ILSCHR ) THEN
               DO 90 J = IFRSTM, ILAST
                  IAROW = INDXG2P(J, NB, 0, 0, NPROW)
                  IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
                     LR = INDXG2L(J, NB, 0, 0, NPROW)
                     A( LR + (LC-1)*LDA ) = -A( LR + (LC-1)*LDA )
                     B( LR + (LC-1)*LDB ) = -B( LR + (LC-1)*LDB )
                  END IF
   90          CONTINUE
            ELSE
               IAROW = INDXG2P(ILAST, NB, 0, 0, NPROW)
               IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
                  LR = INDXG2L(ILAST, NB, 0, 0, NPROW)
                  A( LR + (LC-1)*LDA ) = -A( LR + (LC-1)*LDA )
                  B( LR + (LC-1)*LDB ) = -B( LR + (LC-1)*LDB )
               END IF
            END IF
            IF( ILZ ) THEN
               DO 100 J = 1, N
                  IAROW = INDXG2P(J, NB, 0, 0, NPROW)
                  IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
                     LR = INDXG2L(J, NB, 0, 0, NPROW)
                     Z( LR + (LC-1)*LDZ ) = -Z( LR + (LC-1)*LDZ )
                  END IF
  100          CONTINUE
            END IF
         END IF

         IACOL = INDXG2P(ILAST, NB, 0, 0, NPCOL)
         IAROW = INDXG2P(ILAST, NB, 0, 0, NPROW)
         IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
	      LC = INDXG2L(ILAST, NB, 0, 0, NPCOL)
		  LR = INDXG2L(ILAST, NB, 0, 0, NPROW)

            ALPHAR( ILAST ) = A( LR + (LC-1)*LDA )
            ALPHAI( ILAST ) = ZERO
            BETA( ILAST ) = B( LR + (LC-1)*LDB )
	      WORK(IPTMP) = ALPHAR( ILAST )
		  WORK(IPTMP+1) = ALPHAI( ILAST )
	      WORK(IPTMP+2) = BETA( ILAST )

		  Call BSend(3, 1, ALL,  WORK(IPTMP), LDTMP, ICTXT)
	   ELSE
              Call BRecv(3,1,ALL,WORK(IPTMP),LDTMP,IAROW,IACOL,ICTXT)
	      ALPHAR( ILAST ) =  WORK(IPTMP)
	      ALPHAI( ILAST ) =  WORK(IPTMP+1)
		  BETA( ILAST ) =  WORK(IPTMP+2)
         END IF
*
*        Go to next block -- exit if finished.
*
         ILAST = ILAST - 1

         IF( ILAST.LT.ILO )
     $      GO TO 380
*
*        Reset counters
*
         IITER = 1
         ESHIFT = ZERO
         IF( .NOT.ILSCHR ) THEN
            ILASTM = ILAST
            IF( IFRSTM.GT.ILAST )
     $         IFRSTM = ILO
         END IF
         GO TO 350
*
*        QZ step
*
*        This iteration only involves rows/columns IFIRST:ILAST. We
*        assume IFIRST < ILAST, and that the diagonal of B is non-zero.
*
  110    CONTINUE
         IITER = IITER + 1
         DT = DT + 1
         IF( .NOT.ILSCHR ) THEN
            IFRSTM = IFIRST
         END IF

  200    CONTINUE

         IF ((DEBUG1.EQ.1)) THEN
            WRITE(*,*)'PQZ step3, JITER=', JITER, 
     $           ' IFIRST=', IFIRST, ' ILAST=', ILAST,
     $           ' ILASTM=', ILASTM, ' GOTOL=', GOTOL,
     $		   ' MYROW=', MYROW, 'MYCOL=',MYCOL
         ENDIF

*	   ..
*	   .. Introduce the shifts in A and B at IFIRST
*	   .. We use as much as MAXSHIFTS/2 double implicit shifts
*	   .. or if IFIRST+1=ILAST just one single shift applied to the rows 
*	   .. IFIRST:ILAST of A and B.
*	   ..
	   TMPT1 = MPI_WTIME()
		 
c	   CALL INSHIFTS(A, B, Q, Z,
c     $        ALPHAI, ALPHAR, BETA,
c     $        LDA, LDB, LDQ, LDZ,
c     $         WORK(IPTMP),  WORK(IPTMP2),
c     $		LDTMP,
c     $        DESCA, DESCB, DESCQ, DESCZ,
c     $        IFRSTM, IFIRST, ILAST, ILASTM,
c     $        ILO, NB, N, MAXSHIFTS, SDIST, LDIST, 
c     $        NBULGES, IITER, ESHIFT, MAXIT, 
c     $        ILQ, ILZ, ILSCHR, SAFMIN,
c    $		ALPHAR(P), ALPHAI(P), BETA(P), 
c     $		WORK(IPWORK), JITER,ATOL,WORK(IPU), WORK(IPV), LDTMP,
c     $		AWS)

	   
	   IF (DEBUG0.EQ.1.AND.IAM.EQ.0) THEN
			WRITE(*,*)'Using ', NBULGES, ' bulges'
	   END IF
	   TMPT2 = MPI_WTIME()

	   T4 = T4 + (TMPT2 - TMPT1)
		
         IF ((DEBUG1.EQ.1)) THEN
            WRITE(*,*)'PQZ step4, JITER=', JITER, 
     $           ' IFIRST=', IFIRST, ' ILAST=', ILAST,
     $           ' ILASTM=', ILASTM, ' GOTOL=', GOTOL,
     $           ' NBULGES=', NBULGES,
     $		   ' MYROW=', MYROW, 'MYCOL=',MYCOL
         ENDIF
	
		
         IF (ILAST.LT.ILO) THEN
            GOTO 380
         END IF
         IF (ILAST.LT.IFIRST) THEN
            GOTO 350
         END IF

         IF ((DEBUG1.EQ.1)) THEN
            WRITE(*,*)'PQZ step5, JITER=', JITER, 
     $           ' IFIRST=', IFIRST, ' ILAST=', ILAST,
     $           ' ILASTM=', ILASTM, ' GOTOL=', GOTOL,
     $           ' NBULGES=', NBULGES,
     $		   ' MYROW=', MYROW, 'MYCOL=',MYCOL
         ENDIF
*	   ..
*	   .. Chase the created bulge(s) down the diagonal of 
*	   .. A and B. 
*	   .. We do this when we use double implicit shifts.
*	   ..
		

	   TMPT1 = MPI_WTIME()
		
c         CALL CHASESHIFTS(A, B, Q, Z,
c     $        DESCA, DESCB, DESCQ, DESCZ,
c     $        LDA, LDB, LDQ, LDZ,
c     $        WORK(IPTMP),  WORK(IPTMP2),
c     $		LDTMP,
c     $        NBULGES, NBULGES,
c     $        IFIRST, ILAST, ILASTM, IFRSTM,
c     $        SAFMIN, ILO, NB, N,
c     $        ILQ, ILZ, SDIST, 
c     $	    ALPHAR(P), ALPHAI(P), BETA(P),  
c     $		WORK(IPWORK), JITER,ATOL,
c     $		WORK(IPU), WORK(IPV), LDTMP,
c     $		CHT1, CHT2, CHT3, CHT4, CHT5,
c     $		AWS)				
		
	
		
	   TMPT2 = MPI_WTIME()

	   T5 = T5 + (TMPT2 - TMPT1)

         IF ((DEBUG1.EQ.1)) THEN
            WRITE(*,*)'PQZ step6, JITER=', JITER, 
     $           ' IFIRST=', IFIRST, ' ILAST=', ILAST,
     $           ' ILASTM=', ILASTM, ' GOTOL=', GOTOL,
     $           ' NBULGES=', NBULGES,
     $		   ' MYROW=', MYROW, 'MYCOL=',MYCOL
         ENDIF
		
*
*        End of iteration loop
*
  350    CONTINUE

  360 CONTINUE
*
*     Drop-through = non-convergence
*
  370 CONTINUE
      INFO = ILAST
      GO TO 420
*
*     Successful completion of all QZ steps
*
  380 CONTINUE
*
*     Set Eigenvalues 1:ILO-1
*
      
      DO 410 J = 1, ILO - 1
         IACOL = INDXG2P(J, NB, 0, 0, NPCOL)
	   LC = INDXG2L(J, NB, 0, 0, NPCOL)
         CALL PDLACP3(1, J, B, DESCB, WORK(IPTMP), LDTMP, 
     $      -1, -1, 0)
	   
         IF(WORK(IPTMP).LT.ZERO ) THEN
            IF( ILSCHR ) THEN
               DO 390 JR = 1, J
                  IAROW = INDXG2P(JR, NB, 0, 0, NPROW)
                  IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
                     LR = INDXG2L(JR, NB, 0, 0, NPROW)
                     A( LR + (LC-1)*LDA ) = -A( LR + (LC-1)*LDA )
                     B( LR + (LC-1)*LDB ) = -B( LR + (LC-1)*LDB )
                  END IF
 390           CONTINUE
            ELSE
	         IAROW = INDXG2P(J, NB, 0, 0, NPROW)
			 LR = INDXG2L(J, NB, 0, 0, NPROW)
               IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
                  A( LR + (LC-1)*LDA ) = -A( LR + (LC-1)*LDA )
                  B( LR + (LC-1)*LDB ) = -B( LR + (LC-1)*LDB )
               END IF
            END IF
            IF( ILZ ) THEN
               DO 400 JR = 1, N
                  IAROW = INDXG2P(JR, NB, 0, 0, NPROW)
                  IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN   
                     LR = INDXG2L(JR, NB, 0, 0, NPROW)
                     Z( LR + (LC-1)*LDZ ) = -Z( LR + (LC-1)*LDZ )
                  END IF
 400           CONTINUE
            END IF
         END IF

         IACOL = INDXG2P(J, NB, 0, 0, NPCOL)
         IAROW = INDXG2P(J, NB, 0, 0, NPROW)
         IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
            LR = INDXG2L(J, NB, 0, 0, NPROW)
            LC = INDXG2L(J, NB, 0, 0, NPCOL)
            ALPHAR( J ) = A( LR + (LC-1)*LDA )
            ALPHAI( J ) = ZERO
            BETA( J ) = B( LR + (LC-1)*LDB )
         END IF
 410  CONTINUE
*
*     Normal Termination
*
      INFO = 0
*
*     Exit (other than argument error) -- return optimal workspace size
*
 420  CONTINUE
	IF (DEBUG0.EQ.1.AND.IAM.EQ.0) THEN
		WRITE(*,*)'Number of itterations:',JITER, MAXIT
		WRITE(*,*)'T1=',T1,MYROW, MYCOL
		WRITE(*,*)'T2=',T2,MYROW, MYCOL
		WRITE(*,*)'T3=',T3,MYROW, MYCOL
		WRITE(*,*)'T4=',T4,MYROW, MYCOL
		WRITE(*,*)'T5=',T5,MYROW, MYCOL
		WRITE(*,*)'CHT1=',CHT1,MYROW, MYCOL
		WRITE(*,*)'CHT2=',CHT2,MYROW, MYCOL
		WRITE(*,*)'CHT3=',CHT3,MYROW, MYCOL
		WRITE(*,*)'CHT4=',CHT4,MYROW, MYCOL
		WRITE(*,*)'CHT5=',CHT5,MYROW, MYCOL

	END IF

      WORK( 1 ) = DBLE( N )
 800  FORMAT(20F7.2)
      RETURN
*
*     End of PDHGEQZ
*
      END

	





      SUBROUTINE PDHINFO( SUMMRY, NOUT, NMAT, NVAL, 
     $                      LDNVAL, NNB, NBVAL, LDNBVAL, NGRIDS, PVAL,
     $                      LDPVAL, QVAL, LDQVAL, THRESH, WORK, IAM,
     $                      NPROCS, MAXNOJ, NBULGES, BVAL, LDBVAL )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine (version 1.2) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 10, 1996
*
*     .. Scalar Arguments ..
      INTEGER               IAM, LDNBVAL, LDNVAL, LDPVAL, LDQVAL,
     $                      NGRIDS, NMAT, NNB, NOUT, NPROCS, NBULGES,
     $                      LDBVAL
      REAL                  THRESH
*     ..
*     .. Array Arguments ..
      CHARACTER*( * )       SUMMRY
      INTEGER               NBVAL(LDNBVAL),NVAL(LDNVAL),BVAL( LDBVAL),
     $                      PVAL( LDPVAL ), QVAL( LDQVAL ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  PDHRDINFO get the needed startup information for the Hessenberg
*  reduction tests and transmits it to all processes.
*
*  Arguments
*  =========
*
*  SUMMRY   (global output) CHARACTER*(*)
*           Name of output (summary) file (if any). Only defined for
*           process 0.
*
*  NOUT     (global output) INTEGER
*           The unit number for output file. NOUT = 6, output to screen,
*           NOUT = 0, output to stderr.  Only defined for process 0.
*
*  NMAT     (global output) INTEGER
*           The number of different values that can be used for
*           N, IHI & ILO.
*
*  NVAL     (global output) INTEGER array, dimension (LDNVAL)
*           The values of N (number of rows & columns in matrix).
*
*  LDNVAL   (global input) INTEGER
*           The maximum number of different values that can be used for
*           N, ILO and IHI. LDNVAL >=  NMAT.
*
*  NNB      (global output) INTEGER
*           The number of different values that can be used for NB.
*
*  NBVAL    (global output) INTEGER array, dimension (LDNBVAL)
*           The values of NB (blocksize) to run the code with.
*
*  LDNBVAL  (global input) INTEGER
*           The maximum number of different values that can be used for
*           NB, LDNBVAL >= NNB.
*
*  NGRIDS   (global output) INTEGER
*           The number of different values that can be used for P & Q.
*
*  PVAL     (global output) INTEGER array, dimension (LDPVAL)
*           The values of P (number of process rows) to run the code
*           with.
*
*  LDPVAL   (global input) INTEGER
*           The maximum number of different values that can be used for
*           P, LDPVAL >= NGRIDS.
*
*  QVAL     (global output) INTEGER array, dimension (LDQVAL)
*           The values of Q (number of process columns) to run the code
*           with.
*
*  LDQVAL   (global input) INTEGER
*           The maximum number of different values that can be used for
*           Q, LDQVAL >= NGRIDS.
*
*  THRESH   (global output) REAL
*           Indicates what error checks shall be run and printed out:
*           = 0 : Perform no error checking
*           > 0 : report all residuals greater than THRESH.
*
*  WORK     (local workspace) INTEGER array, dimension >=
*           3*LDNVAL+LDNBVAL+2*LDPVAL.  Used to pack all input arrays
*           in order to send info in one message.
*
*  IAM      (local input) INTEGER
*           My process number.
*
*  NPROCS   (global input) INTEGER
*           The total number of processes.
*
*  Note
*  ====
*
*  For packing the information we assumed that the length in bytes of an
*  integer is equal to the length in bytes of a real single precision.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            NIN
      PARAMETER          ( NIN = 11 )
*     ..
*     .. Local Scalars ..
      CHARACTER*79       USRINFO
      INTEGER            I, ICTXT, MAXNOJ
      DOUBLE PRECISION   EPS
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_ABORT, BLACS_GET, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINIT, BLACS_SETUP, ICOPY, IGEBR2D,
     $                   IGEBS2D, SGEBR2D, SGEBS2D
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   PDLAMCH
      EXTERNAL           PDLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Process 0 reads the input data, broadcasts to other processes and
*     writes needed information to NOUT
*
      IF( IAM.EQ.0 ) THEN
*
*        Open file and skip data file header
*
         OPEN( UNIT = NIN, FILE = 'H.dat', STATUS = 'OLD' )
         READ( NIN, FMT = * )SUMMRY
         SUMMRY = ' '
*
*        Read in user-supplied info about machine type, compiler, etc.
*
         READ( NIN, FMT = * ) USRINFO
*
*        Read name and unit number for summary output file
*
         READ( NIN, FMT = * ) SUMMRY
         READ( NIN, FMT = * ) MAXNOJ
         READ( NIN, FMT = * ) NOUT
         IF( NOUT.NE.0 .AND. NOUT.NE.6 )
     $      OPEN( UNIT = NOUT, FILE = SUMMRY, STATUS = 'UNKNOWN' )
*
*        Read and check the parameter values for the tests.
*
*        Get number of matrices
*
         READ( NIN, FMT = * ) NMAT
         IF( NMAT.LT.1. .OR. NMAT.GT.LDNVAL ) THEN
            WRITE( NOUT, FMT = 9997 ) 'N', LDNVAL
            GO TO 20
         END IF
*
*        Get values of N, ILO, IHI
*
         READ( NIN, FMT = * ) ( NVAL( I ), I = 1, NMAT )
*
*        Get values of NB
*
         READ( NIN, FMT = * ) NNB
         IF( NNB.LT.1 .OR. NNB.GT.LDNBVAL ) THEN
            WRITE( NOUT, FMT = 9997 ) 'NB', LDNBVAL
            GO TO 20
         END IF
         READ( NIN, FMT = * ) ( NBVAL( I ), I = 1, NNB )
*
*        Get number of grids
*
         READ( NIN, FMT = * ) NGRIDS
         IF( NGRIDS.LT.1 .OR. NGRIDS.GT.LDPVAL ) THEN
            WRITE( NOUT, FMT = 9997 ) 'Grids', LDPVAL
            GO TO 20
         ELSE IF( NGRIDS.GT.LDQVAL ) THEN
            WRITE( NOUT, FMT = 9997 ) 'Grids', LDQVAL
            GO TO 20
         END IF
*
*        Get values of P and Q
*
         READ( NIN, FMT = * ) ( PVAL( I ), I = 1, NGRIDS )
         READ( NIN, FMT = * ) ( QVAL( I ), I = 1, NGRIDS )
*
*	   Get values of Bulges
*
         READ( NIN, FMT = * ) NBULGES
         IF( NBULGES.LT.1 .OR. NBULGES.GT.LDBVAL ) THEN
            WRITE( NOUT, FMT = 9997 ) 'NB', LDBVAL
            GO TO 20
         END IF
         READ( NIN, FMT = * ) ( BVAL( I ), I = 1, NBULGES )


*
*        Get level of checking
*
         READ( NIN, FMT = * ) THRESH
*
*        Close input file
*
         CLOSE( NIN )
*
*        For pvm only: if virtual machine not set up, allocate it and
*        spawn the correct number of processes.
*
         IF( NPROCS.LT.1 ) THEN
            NPROCS = 0
            DO 10 I = 1, NGRIDS
               NPROCS = MAX( NPROCS, PVAL( I )*QVAL( I ) )
   10       CONTINUE
            CALL BLACS_SETUP( IAM, NPROCS )
         END IF
*
*        Temporarily define blacs grid to include all processes so
*        information can be broadcast to all processes
*
         CALL BLACS_GET( -1, 0, ICTXT )
         CALL BLACS_GRIDINIT( ICTXT, 'Row-major', 1, NPROCS )
*
*        Compute machine epsilon
*
         EPS = PDLAMCH( ICTXT, 'eps' )
*
*        Pack information arrays and broadcast
*
         IF ( NPROCS .GT. 1) THEN
            CALL SGEBS2D( ICTXT, 'All', ' ', 1, 1, THRESH, 1 )
         ENDIF
*
         WORK( 1 ) = NMAT
         WORK( 2 ) = NNB
         WORK( 3 ) = NGRIDS
         WORK( 4 ) = MAXNOJ
	   WORK( 5 ) = NBULGES

         IF ( NPROCS .GT. 1) THEN
            CALL IGEBS2D( ICTXT, 'All', ' ', 1, 5, WORK, 1 )
         ENDIF
*
         I = 1
         CALL ICOPY( NMAT, NVAL, 1, WORK( I ), 1 )
         I = I + NMAT
         CALL ICOPY( NNB, NBVAL, 1, WORK( I ), 1 )
         I = I + NNB
         CALL ICOPY( NGRIDS, PVAL, 1, WORK( I ), 1 )
         I = I + NGRIDS
         CALL ICOPY( NGRIDS, QVAL, 1, WORK( I ), 1 )
         I = I + NGRIDS
         CALL ICOPY( NBULGES, BVAL, 1, WORK( I ), 1 )
         I = I + NBULGES -1

         IF ( NPROCS .GT. 1) THEN
            CALL IGEBS2D( ICTXT, 'All', ' ', 1, I, WORK, 1 )
         ENDIF
*
*        regurgitate input
*
         WRITE( NOUT, FMT = 9999 )
     $      'ScaLAPACK Reduction routine to Block H-T form.'
         WRITE( NOUT, FMT = 9999 ) USRINFO
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 )
     $               'Tests of the parallel '//
     $               'real double precision block H-T'
         WRITE( NOUT, FMT = 9999 ) 'reduction routines.'
         WRITE( NOUT, FMT = 9999 )
     $               'The following scaled residual '//
     $               'checks will be computed:'
         WRITE( NOUT, FMT = 9999 )
     $               ' ||??????|| / (||???* eps * N)'
         WRITE( NOUT, FMT = 9999 )
     $               'The matrises A and B are randomly '//
     $               'generated for each test.'
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 )
     $               'An explanation of the input/output '//
     $               'parameters follows:'
         WRITE( NOUT, FMT = 9999 )
     $               'TIME     : Indicates whether WALL or '//
     $               'CPU time was used.'
         WRITE( NOUT, FMT = 9999 )
     $               'N        : The number of rows and columns '//
     $               'of the matrix A.'
         WRITE( NOUT, FMT = 9999 )
     $               'NB       : The size of the square blocks'//
     $               ' the matrix A is split into.'
         WRITE( NOUT, FMT = 9999 )
     $      '         on to the next column of processes.'
         WRITE( NOUT, FMT = 9999 )
     $               'P        : The number of process rows.'
         WRITE( NOUT, FMT = 9999 )
     $               'Q        : The number of process columns.'
         WRITE( NOUT, FMT = 9999 )
     $               'HRD time : Time in seconds to compute HRD '
         WRITE( NOUT, FMT = 9999 )
     $               'MFLOPS   : Rate of execution for HRD ' //
     $               'reduction.'
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9999 )
     $               'The following parameter values will be used:'
         WRITE( NOUT, FMT = 9995 )
     $               'N    ', ( NVAL( I ), I = 1, MIN( NMAT, 10 ) )
         IF( NMAT.GT.10 )
     $      WRITE( NOUT, FMT = 9994 ) ( NVAL( I ), I = 11, NMAT )
         WRITE( NOUT, FMT = 9995 )
     $               'NB   ', ( NBVAL( I ), I = 1, MIN( NNB, 10 ) )
         IF( NNB.GT.10 )
     $      WRITE( NOUT, FMT = 9994 ) ( NBVAL( I ), I = 11, NNB )
         WRITE( NOUT, FMT = 9995 )
     $               'P    ', ( PVAL( I ), I = 1, MIN( NGRIDS, 10 ) )
         IF( NGRIDS.GT.10 )
     $      WRITE( NOUT, FMT = 9994 ) ( PVAL( I ), I = 11, NGRIDS )
         WRITE( NOUT, FMT = 9995 )
     $               'Q    ', ( QVAL( I ), I = 1, MIN( NGRIDS, 10 ) )
         IF( NGRIDS.GT.10 )
     $      WRITE( NOUT, FMT = 9994 ) ( QVAL( I ), I = 11, NGRIDS )
         WRITE( NOUT, FMT = 9995 )
     $               'MAXBULGES    ',( BVAL( I ),I = 1,MIN(NBULGES,10))

         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9996 ) EPS
         WRITE( NOUT, FMT = 9993 ) THRESH
*
      ELSE
*
*        If in pvm, must participate setting up virtual machine
*
         IF( NPROCS.LT.1 )
     $      CALL BLACS_SETUP( IAM, NPROCS )
*
*        Temporarily define blacs grid to include all processes so
*        all processes have needed startup information
*
         CALL BLACS_GET( -1, 0, ICTXT )
         CALL BLACS_GRIDINIT( ICTXT, 'Row-major', 1, NPROCS )
*
*        Compute machine epsilon
*
         EPS = PDLAMCH( ICTXT, 'eps' )
*
         CALL SGEBR2D( ICTXT, 'All', ' ', 1, 1, THRESH, 1, 0, 0 )
         CALL IGEBR2D( ICTXT, 'All', ' ', 1, 5, WORK, 1, 0, 0 )
         NMAT   = WORK( 1 )
         NNB    = WORK( 2 )
         NGRIDS = WORK( 3 )
         MAXNOJ = WORK( 4 )
	   NBULGES = WORK( 5 )
*
         I = 1 + NMAT + NNB + 2*NGRIDS + NBULGES - 1
         CALL IGEBR2D( ICTXT, 'All', ' ', 1, I, WORK, 1, 0, 0 )
*
         I = 1
         CALL ICOPY( NMAT, WORK( I ), 1, NVAL, 1 )
         I = I + NMAT
         CALL ICOPY( NNB, WORK( I ), 1, NBVAL, 1 )
         I = I + NNB
         CALL ICOPY( NGRIDS, WORK( I ), 1, PVAL, 1 )
         I = I + NGRIDS
         CALL ICOPY( NGRIDS, WORK( I ), 1, QVAL, 1 )
         I = I + NGRIDS
         CALL ICOPY( NBULGES, WORK( I ), 1, BVAL, 1 )

*
      END IF
*
	  CALL BLACS_BARRIER(ICTXT, 'A')
	  CALL BLACS_BARRIER(ICTXT, 'C')
	  CALL BLACS_BARRIER(ICTXT, 'R')
	  CALL BLACS_BARRIER(ICTXT, 'ALL')
      CALL BLACS_GRIDEXIT( ICTXT )
*
      RETURN
*
   20 CONTINUE
      WRITE( NOUT, FMT = 9998 )
      CLOSE( NIN )
      IF( NOUT.NE.6 .AND. NOUT.NE.0 )
     $   CLOSE( NOUT )
      CALL BLACS_ABORT( ICTXT, 1 )
*
      STOP
*
 9999 FORMAT( A )
 9998 FORMAT( ' ILLEGAL INPUT IN FILE ', 40A, '.  ABORTING RUN.' )
 9997 FORMAT( ' NUMBER OF VALUES OF ', 5A,
     $      ' IS LESS THAN 1 OR GREATER ', 'THAN ', I2 )
 9996 FORMAT( 'Relative machine precision (eps) is taken to be ',
     $      E18.6 )
 9995 FORMAT( 2X, A5, ':        ', 10I6 )
 9994 FORMAT( '              ', 10I6 )
 9993 FORMAT( 'Routines pass computational tests if scaled residual is',
     $      ' less than ', G14.7 )
*
*     End of PDHRDINFO
*
      END
      SUBROUTINE PDLAREF( TYPE, A, LDA, N, NB, DESCA, BLOCK,
     $     IROW1, ICOL1, ISTART, ISTOP, ITMP1, ITMP2,
     $     VECS, V2, V3, T1, T2, T3,
     $     WORK, RET, LDTMP)

      IMPLICIT NONE
*
*
*     .. Scalar Arguments ..
      INTEGER            LDA, N, NB
      INTEGER            IROW1, ICOL1, ISTART, ISTOP
      INTEGER            ITMP1, ITMP2
      INTEGER            RET, LDTMP
      DOUBLE PRECISION   V2, V3, T1, T2, T3
      LOGICAL            BLOCK
      CHARACTER          TYPE

*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), VECS( * )
      INTEGER            DESCA(*)
      DOUBLE PRECISION   WORK(LDTMP, *)
*     ..
*

*  Purpose
*  =======
*
*  DLAREF applies one or several Householder reflectors of size 3
*     to one or two matrices (if column is specified) on either their
*     rows or columns.
*
*  Arguments
*  =========
*
*  TYPE    (global input) CHARACTER*1
*          If 'R': Apply reflectors to the rows of the matrix
*              (apply from left)
*          Otherwise: Apply reflectors to the columns of the matrix
*          Unchanged on exit.
*
*  A       (global input/output) DOUBLE PRECISION array, (LDA,*)
*          On entry, the matrix to receive the reflections.
*          The updated matrix on exit.
*
*  LDA     (local input) INTEGER
*          On entry, the leading dimension of the local A.  Unchanged on exit.
*
* 
*  N       (local input) INTEGER
*          On entry, the size of the global A. Unchanged on exit.
*     
*
*  DESCA   ...
*
*
*  NB      On entry, the block size of the global A. Unchanged on exit.
*
*
*  BLOCK   (global input) LOGICAL
*          If .TRUE., then apply several reflectors at once and read
*             their data from the VECS array.
*          If .FALSE., apply the single reflector given by V2, V3,
*             T1, T2, and T3.
*
*  IROW1   (local input/output) INTEGER
*          On entry, the local row element of A.
*          Undefined on output.
*
*
*  ICOL1   (local input/output) INTEGER
*          On entry, the local column element of A.
*          Undefined on output.
*
*  ISTART  (global input) INTEGER
*          Specifies the "number" of the first reflector.  This is
*              used as an index into VECS if BLOCK is set.
*              ISTART is ignored if BLOCK is .FALSE..
*
*  ISTOP   (global input) INTEGER
*          Specifies the "number" of the last reflector.  This is
*              used as an index into VECS if BLOCK is set.
*              ISTOP is ignored if BLOCK is .FALSE..
*
*  ITMP1   (local input) INTEGER
*          Starting range into A.  For rows, this is the local
*              first column.  For columns, this is the local first row.
*
*  ITMP2   (local input) INTEGER
*          Ending range into A.  For rows, this is the local last
*              column.  For columns, this is the local last row.
*
*  VECS    (global input) DOUBLE PRECISION array of size 3*N (matrix
*                                                             size)
*          This holds the size 3 reflectors one after another and this
*              is only accessed when BLOCK is .TRUE.
*
*  V2
*  V3
*  T1
*  T2
*  T3      (global input/output) DOUBLE PRECISION
*          This holds information on a single size 3 Householder
*              reflector and is read when BLOCK is .FALSE., and
*              overwritten when BLOCK is .TRUE.
*
*  Implemented by:  B. Adlerborn, August 22, 2002
*
*  =====================================================================
*
*     ..
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      INTEGER            DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
      INTEGER            ROWWISE, COLUMNWISE, ALL
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2, ALL=3)

      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, K
      INTEGER            IAROW1, IACOL1
	INTEGER			   LROW1, LCOL1
	INTEGER			   LROW2, LCOL2
	INTEGER			   LROW3, LCOL3

      INTEGER            IAROW2, IACOL2
      INTEGER            RCNT, DCNT, IR1, IC1, SCNT, LJ
	LOGICAL			   B1
      INTEGER            ICTXT, MYROW, MYCOL, NPROW, NPCOL 
      DOUBLE PRECISION   SUM
*     ..
*     .. External Subroutines ..
*     ..      
      INTEGER            INDXG2L, INDXG2P
      EXTERNAL           INDXG2L, INDXG2P

*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*     

	ICTXT = DESCA(CTXT_)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )  

      IF ((ISTART.GT.ISTOP ).OR.(ITMP1.GT.ITMP2 )) THEN
         RET = -1
         RETURN
      END IF


      IR1 = IROW1
      IC1 = ICOL1
      IF(TYPE.EQ.'R'  ) THEN
            DO 40 K = ISTART, ISTOP
               V2 = VECS( ( K-1 )*3+1 )
               V3 = VECS( ( K-1 )*3+2 )
               T1 = VECS( ( K-1 )*3+3 )
               T2 = T1*V2
               T3 = T1*V3

               IAROW1 = INDXG2P(IR1, NB, 0, 0, NPROW)
               IAROW2 = INDXG2P(IR1+2, NB, 0, 0, NPROW)
			 IF ((MYROW.NE.IAROW1).AND.(MYROW.NE.IAROW2)) GOTO 34
	         IF (IAROW1.EQ.IAROW2) THEN
                  DCNT = 0
                  LROW1 = INDXG2L(IR1, NB, 0, 0, NPROW)
			    LROW2 = LROW1 + 1
			    LROW3 = LROW1 + 2
               ELSE IF (IAROW1.EQ.INDXG2P(IR1+1, NB, 0, 0, NPROW)) THEN
                  DCNT = 1
                  LROW1 = INDXG2L(IR1, NB, 0, 0, NPROW)
			    LROW2 = LROW1 + 1
			    LROW3 = INDXG2L(IR1+2, NB, 0, 0, NPROW)
               ELSE
				DCNT = 2
                  LROW1 = INDXG2L(IR1, NB, 0, 0, NPROW)
			    LROW2 = INDXG2L(IR1+1, NB, 0, 0, NPROW)
			    LROW3 = LROW2 + 1
               END IF
			 IF (IAROW1.EQ.IAROW2) THEN
				DO 28 J = ITMP1, ITMP2
					B1 = ((J.EQ.ITMP1).OR.(MOD(J-1,NB).EQ.0) )			 
					IF (B1) THEN
						LCOL1 = INDXG2L(J, NB, 0, 0, NPCOL)
						IACOL1 = INDXG2P(J, NB, 0, 0, NPCOL)
					ELSE
						LCOL1 = LCOL1 + 1                      
					END IF
					IF (MYCOL.EQ.IACOL1) THEN
                        SUM = A( LROW1, LCOL1 ) + 
     $                       V2*A( LROW2, LCOL1 ) +
     $                       V3*A( LROW3, LCOL1 )
                        A( LROW1, LCOL1 ) = A( LROW1, LCOL1 ) - SUM*T1
                        A( LROW2, LCOL1 ) = A( LROW2, LCOL1 ) - SUM*T2
                        A( LROW3, LCOL1 ) = A( LROW3, LCOL1 ) - SUM*T3
                      END IF
 28			    CONTINUE
				GOTO 34
			 END IF
			 IF (MYROW.EQ.IAROW2) THEN
				LROW1 = INDXG2L(IR1+(3-DCNT), NB, 0, 0, NPROW)
			 END IF
			 IF (DCNT.EQ.2) THEN 
			   DO 30 J = ITMP1, ITMP2
			    B1 = ((J.EQ.ITMP1).OR.(MOD(J-1,NB).EQ.0) )
				IF (B1) THEN
					LJ = MOD(J-1, NB)+1
					LCOL1 = INDXG2L(J, NB, 0, 0, NPCOL)
					IACOL1 = INDXG2P(J, NB, 0, 0, NPCOL)               
				ELSE		
					LCOL1 = LCOL1 + 1
					LJ = LJ + 1
				END IF
				IF (MYCOL.NE.IACOL1) GOTO 30
*                    ..
*                    .. Copy elements and apply rotations ..
*                    ..
                     IF (MYROW.EQ.IAROW1) THEN	  
					  IF (B1) THEN
				           SCNT = MIN(NB - MOD(J-1, NB), ITMP2-J+1)
						   Call Send(1, SCNT, DOWN, 
     $							A(LROW1, LCOL1), LDA,
     $							MYROW, MYCOL, NPCOL, NPROW, ICTXT)
						   Call Recv(2, SCNT,DOWN, 
     $							WORK(1, LJ), LDTMP, 
     $							MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				      END IF

					  SUM = A( LROW1, LCOL1 ) + V2*WORK( 1, LJ ) +
     $                         V3*WORK( 2, LJ )
                        A( LROW1, LCOL1 ) = A( LROW1, LCOL1 )-SUM*T1
                     ELSE                       
					  LROW2 = LROW1 + 1					  
					  IF (B1) THEN
							SCNT = MIN(NB - MOD(J-1, NB), ITMP2-J+1)
							Call Send(2, SCNT, UP, 
     $							A(LROW1, LCOL1), LDA,
     $							MYROW, MYCOL, NPCOL, NPROW, ICTXT)
							Call Recv(1, SCNT, UP, 
     $							WORK(1, LJ), LDTMP, 
     $							MYROW, MYCOL, NPCOL, NPROW, ICTXT)			
				      END IF

					  SUM = WORK( 1, LJ ) + 
     $                         V2*A( LROW1, LCOL1 ) +
     $                         V3*A( LROW2, LCOL1 )
                        A( LROW1, LCOL1 ) = A( LROW1, LCOL1 )-SUM*T2 
                        A( LROW2, LCOL1 ) = A( LROW2, LCOL1 )-SUM*T3
                     END IF 
 30              CONTINUE
			 ELSE

			   DO 31 J = ITMP1, ITMP2
				B1 = ((J.EQ.ITMP1).OR.(MOD(J-1,NB).EQ.0) )
				IF (B1) THEN
					LJ = MOD(J-1, NB)+1
					LCOL1 = INDXG2L(J, NB, 0, 0, NPCOL)
					IACOL1 = INDXG2P(J, NB, 0, 0, NPCOL)               
				ELSE
					LJ = LJ + 1
					LCOL1 = LCOL1 + 1
				END IF
				IF (MYCOL.NE.IACOL1) GOTO 31
*                    ..
*                    .. Copy elements and apply rotations ..
*                    ..

                     IF (MYROW.EQ.IAROW1) THEN	  
					  IF (B1) THEN
				           SCNT = MIN(NB - MOD(J-1, NB), ITMP2-J+1)
						   Call Send(2, SCNT, DOWN, 
     $							A(LROW1, LCOL1), LDA,
     $							MYROW, MYCOL, NPCOL, NPROW, ICTXT)
						   Call Recv(1, SCNT,DOWN, 
     $							WORK(1, LJ), LDTMP, 
     $							MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				      END IF


                        SUM = A( LROW1, LCOL1 ) +
     $                       V2*A( LROW2, LCOL1) +
     $                       V3*WORK( 1, LJ )
                        A( LROW1, LCOL1 ) = A( LROW1, LCOL1 )-SUM*T1
                        A( LROW2, LCOL1 ) = A( LROW2, LCOL1 )-SUM*T2
                     ELSE                        				  
					  IF (B1) THEN
							SCNT = MIN(NB - MOD(J-1, NB), ITMP2-J+1)
							Call Send(1, SCNT, UP, 
     $							A(LROW1, LCOL1), LDA,
     $							MYROW, MYCOL, NPCOL, NPROW, ICTXT)
							Call Recv(2, SCNT, UP, 
     $							WORK(1, LJ), LDTMP, 
     $							MYROW, MYCOL, NPCOL, NPROW, ICTXT)			
				      END IF

                        SUM = WORK( 1, LJ ) + V2*WORK( 2, LJ ) +
     $                      V3*A( LROW1, LCOL1 )
					  A( LROW1, LCOL1 ) = A( LROW1, LCOL1 ) - SUM*T3	
                     END IF
 
 31              CONTINUE
			 END IF
 34		     CONTINUE
               IR1 = IR1 + 1

 40         CONTINUE
      ELSE
*     
*     Do column transforms
*
            DO 110 K = ISTART, ISTOP
               V2 = VECS( ( K-1 )*3+1 )
               V3 = VECS( ( K-1 )*3+2 )
               T1 = VECS( ( K-1 )*3+3 )
               T2 = T1*V2
               T3 = T1*V3     
			              
               IACOL1 = INDXG2P(IC1, NB, 0, 0, NPCOL)
               IACOL2 = INDXG2P(IC1+2, NB, 0, 0, NPCOL)
			 IF ((MYCOL.NE.IACOL1).AND.(MYCOL.NE.IACOL2)) GOTO 94
               IF (IACOL1.EQ.IACOL2) THEN
                  RCNT = 0
				LCOL1 = INDXG2L(IC1, NB, 0, 0, NPCOL)
				LCOL2 = LCOL1 + 1
				LCOL3 = LCOL1 + 2
               ELSE IF (IACOL1.EQ.INDXG2P(IC1+1, NB, 0, 0, NPCOL)) 
     $                 THEN
                  RCNT = 1
				LCOL1 = INDXG2L(IC1, NB, 0, 0, NPCOL)
				LCOL2 = LCOL1 + 1
				LCOL3 = INDXG2L(IC1+2, NB, 0, 0, NPCOL)
               ELSE
                  RCNT = 2
				LCOL1 = INDXG2L(IC1, NB, 0, 0, NPCOL)
				LCOL2 = INDXG2L(IC1+1, NB, 0, 0, NPCOL)
				LCOL3 = LCOL2 + 1
               END IF

			 IF (IACOL1.EQ.IACOL2) THEN				
				DO 88 J = ITMP1, ITMP2
					B1 = ((J.EQ.ITMP1).OR.(MOD(J-1,NB).EQ.0))
					IF (B1) THEN
						IAROW1 = INDXG2P(J, NB, 0, 0, NPROW)
						LROW1 = INDXG2L(J, NB, 0, 0, NPROW)				
                      ELSE
						LROW1 = LROW1 + 1				
					END IF
					IF (MYROW.EQ.IAROW1) THEN
                        SUM = A( LROW1, LCOL1 ) + V2*A( LROW1, LCOL2 )
     $                       + V3*A( LROW1, LCOL3 )
                        A( LROW1, LCOL1 ) = A( LROW1, LCOL1 ) - SUM*T1
                        A( LROW1, LCOL2 ) = A( LROW1, LCOL2 ) - SUM*T2
                        A( LROW1, LCOL3 ) = A( LROW1, LCOL3 ) - SUM*T3
                      END IF
88				CONTINUE
				GOTO 94
			 END IF
			 IF (MYCOL.EQ.IACOL2) THEN
				LCOL1 = INDXG2L(IC1+(3-RCNT), NB, 0, 0, NPCOL)
			 END IF
			 IF (RCNT.EQ.2) THEN
			   
                 DO 90 J = ITMP1, ITMP2
				B1 = ((J.EQ.ITMP1).OR.(MOD(J-1,NB).EQ.0))
				IF (B1) THEN
					LJ = MOD(J-1, NB)+1
					IAROW1 = INDXG2P(J, NB, 0, 0, NPROW)
					LROW1 = INDXG2L(J, NB, 0, 0, NPROW)
				ELSE
					LROW1 = LROW1 + 1
					LJ = LJ + 1
				END IF
				IF (MYROW.NE.IAROW1) GOTO 90
*                    ..
*                    .. Copy elements and apply rotations..
*                    ..
                     IF (MYCOL.EQ.IACOL1) THEN					  
					  IF (B1) THEN
							SCNT = MIN(NB - MOD(J-1, NB), ITMP2-J+1)
							Call Send(SCNT, 1, RIGHT, 
     $							A(LROW1, LCOL1), LDA,
     $							MYROW, MYCOL, NPCOL, NPROW, ICTXT)
							Call Recv(SCNT, 2, RIGHT, 
     $							WORK(LJ,1), LDTMP, 
     $							MYROW, MYCOL, NPCOL, NPROW, ICTXT)		
				      END IF
                        SUM = A( LROW1, LCOL1 ) + V2*WORK( LJ, 1 ) +
     $                        V3*WORK( LJ, 2 )
					  A( LROW1, LCOL1 ) = A( LROW1, LCOL1 ) - SUM*T1
                     ELSE                        						
                        LCOL2 = LCOL1 + 1					  
					  IF (B1) THEN
							SCNT = MIN(NB - MOD(J-1, NB), ITMP2-J+1)
							Call Send(SCNT, 2, LEFT, 
     $							A(LROW1, LCOL1), LDA,
     $							MYROW, MYCOL, NPCOL, NPROW, ICTXT)
							Call Recv(SCNT, 1, LEFT, 
     $							WORK(LJ,1), LDTMP, 
     $							MYROW, MYCOL, NPCOL, NPROW, ICTXT)
					  END IF
                        SUM = WORK(LJ, 1) + V2*A(LROW1, LCOL1) + 
     $                       V3*A(LROW1, LCOL2)
                        A( LROW1, LCOL1 ) = A( LROW1, LCOL1 )-SUM*T2
                        A( LROW1, LCOL2 ) = A( LROW1, LCOL2 )-SUM*T3
                     END IF                     
 90              CONTINUE
			 ELSE
                 DO 91 J = ITMP1, ITMP2
				
	            B1 = ((J.EQ.ITMP1).OR.(MOD(J-1,NB).EQ.0))
				IF (B1) THEN
					LJ = MOD(J-1, NB)+1
		            IAROW1 = INDXG2P(J, NB, 0, 0, NPROW)
			        LROW1 = INDXG2L(J, NB, 0, 0, NPROW)
				ELSE
					LJ = LJ +1
			        LROW1 = LROW1 + 1				
				END IF
				IF (MYROW.NE.IAROW1) GOTO 91
*                    ..
*                    .. Copy elements and apply rotations..
*                    ..
                     IF (MYCOL.EQ.IACOL1) THEN					  
					  IF (B1) THEN
							SCNT = MIN(NB - MOD(J-1, NB), ITMP2-J+1)
							Call Send(SCNT, 2, RIGHT, 
     $							A(LROW1, LCOL1), LDA,
     $							MYROW, MYCOL, NPCOL, NPROW, ICTXT)
							Call Recv(SCNT, 1, RIGHT, 
     $							WORK(LJ,1), LDTMP, 
     $							MYROW, MYCOL, NPCOL, NPROW, ICTXT)		
				      END IF
                        SUM = A( LROW1, LCOL1 ) + 
     $                       V2*A( LROW1, LCOL2 ) +
     $                       V3*WORK( LJ, 1 )
                        A( LROW1, LCOL1 ) = A( LROW1, LCOL1 )-SUM*T1
                        A( LROW1, LCOL2 ) = A( LROW1, LCOL2 )-SUM*T2 
                     ELSE                        			  
					  IF (B1) THEN
							SCNT = MIN(NB - MOD(J-1, NB), ITMP2-J+1)
							Call Send(SCNT, 1, LEFT, 
     $							A(LROW1, LCOL1), LDA,
     $							MYROW, MYCOL, NPCOL, NPROW, ICTXT)
							Call Recv(SCNT, 2, LEFT, 
     $							WORK(LJ,1), LDTMP, 
     $							MYROW, MYCOL, NPCOL, NPROW, ICTXT)
					  END IF
                        SUM = WORK(LJ, 1) + V2*WORK(LJ, 2) + 
     $                       V3*A(LROW1, LCOL1)
                        A( LROW1, LCOL1 ) = A( LROW1, LCOL1 )-SUM*T3

                     END IF                     
 91              CONTINUE
		     END IF
 94			 CONTINUE
               IC1 = IC1 + 1
 110        CONTINUE

      END IF

      RETURN
*
*     End of PDLAREF
*
      END
      SUBROUTINE PDMATGEN( ICTXT, AFORM, DIAG, M, N, MB, NB, A, LDA,
     $                     IAROW, IACOL, ISEED, IROFF, IRNUM, ICOFF,
     $                     ICNUM, MYROW, MYCOL, NPROW, NPCOL )
*
*  -- ScaLAPACK routine (version 1.4) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 17, 1996
*
*     .. Scalar Arguments ..
      CHARACTER*1        AFORM, DIAG
      INTEGER            IACOL, IAROW, ICNUM, ICOFF, ICTXT, IRNUM,
     $                   IROFF, ISEED, LDA, M, MB, MYCOL, MYROW, N,
     $                   NB, NPCOL, NPROW
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  PDMATGEN : Parallel Real Double precision MATrix GENerator.
*  Generate (or regenerate) a distributed matrix A (or sub-matrix of A).
*
*  Arguments
*  =========
*
*  ICTXT   (global input) INTEGER
*          The BLACS context handle, indicating the global context of
*          the operation. The context itself is global.
*
*  AFORM   (global input) CHARACTER*1
*          if AFORM = 'S' : A is returned is a symmetric matrix.
*          if AFORM = 'H' : A is returned is a Hermitian matrix.
*          if AFORM = 'T' : A is overwritten with the transpose of
*                           what would normally be generated.
*          if AFORM = 'C' : A is overwritten with the conjugate trans-
*                           pose of what would normally be generated.
*          otherwise a random matrix is generated.
*
*  DIAG    (global input) CHARACTER*1
*          if DIAG = 'D' : A is diagonally dominant.
*
*  M       (global input) INTEGER
*          The number of rows in the generated distributed matrix.
*
*  N       (global input) INTEGER
*          The number of columns in the generated distributed
*          matrix.
*
*  MB      (global input) INTEGER
*          The row blocking factor of the distributed matrix A.
*
*  NB      (global input) INTEGER
*          The column blocking factor of the distributed matrix A.
*
*  A       (local output) DOUBLE PRECISION, pointer into the local
*          memory to an array of dimension ( LDA, * ) containing the
*          local pieces of the distributed matrix.
*
*  LDA     (local input) INTEGER
*          The leading dimension of the array containing the local
*          pieces of the distributed matrix A.
*
*  IAROW   (global input) INTEGER
*          The row processor coordinate which holds the first block
*          of the distributed matrix A.
*
*  IACOL   (global input) INTEGER
*          The column processor coordinate which holds the first
*          block of the distributed matrix A.
*
*  ISEED   (global input) INTEGER
*          The seed number to generate the distributed matrix A.
*
*  IROFF   (local input) INTEGER
*          The number of local rows of A that have already been
*          generated.  It should be a multiple of MB.
*
*  IRNUM   (local input) INTEGER
*          The number of local rows to be generated.
*
*  ICOFF   (local input) INTEGER
*          The number of local columns of A that have already been
*          generated.  It should be a multiple of NB.
*
*  ICNUM   (local input) INTEGER
*          The number of local columns to be generated.
*
*  MYROW   (local input) INTEGER
*          The row process coordinate of the calling process.
*
*  MYCOL   (local input) INTEGER
*          The column process coordinate of the calling process.
*
*  NPROW   (global input) INTEGER
*          The number of process rows in the grid.
*
*  NPCOL   (global input) INTEGER
*          The number of process columns in the grid.
*
*  Notes
*  =====
*
*  The code is originally developed by David Walker, ORNL,
*  and modified by Jaeyoung Choi, ORNL.
*
*  Reference: G. Fox et al.
*  Section 12.3 of "Solving problems on concurrent processors Vol. I"
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            MULT0, MULT1, IADD0, IADD1
      PARAMETER        ( MULT0=20077, MULT1=16838, IADD0=12345,
     $                   IADD1=0 )
      DOUBLE PRECISION   ONE, TWO
      PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            SYMM, HERM, TRAN
      INTEGER            I, IC, IK, INFO, IOFFC, IOFFR, IR, J, JK,
     $                   JUMP1, JUMP2, JUMP3, JUMP4, JUMP5, JUMP6,
     $                   JUMP7, MAXMN, MEND, MOFF, MP, MRCOL, MRROW,
     $                   NEND, NOFF, NPMB, NQ, NQNB
*     ..
*     .. Local Arrays ..
      INTEGER            IADD(2), IA1(2), IA2(2), IA3(2), IA4(2),
     $                   IA5(2), IB1(2), IB2(2), IB3(2), IC1(2), IC2(2),
     $                   IC3(2), IC4(2), IC5(2), IRAN1(2), IRAN2(2),
     $                   IRAN3(2), IRAN4(2), ITMP1(2), ITMP2(2),
     $                   ITMP3(2), JSEED(2), MULT(2)
*     ..
*     .. External Subroutines ..
      EXTERNAL           JUMPIT, PXERBLA, SETRAN, XJUMPM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MOD
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICEIL, NUMROC
      DOUBLE PRECISION   PDRAND
      EXTERNAL           ICEIL, NUMROC, LSAME, PDRAND
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      MP   = NUMROC( M, MB, MYROW, IAROW, NPROW )
      NQ   = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
      SYMM = LSAME( AFORM, 'S' )
      HERM = LSAME( AFORM, 'H' )
      TRAN = LSAME( AFORM, 'T' )
*
      INFO = 0
      IF( .NOT.LSAME( DIAG, 'D' ) .AND.
     $         .NOT.LSAME( DIAG, 'N' )        ) THEN
         INFO = 3
      ELSE IF( SYMM.OR.HERM ) THEN
         IF( M.NE.N ) THEN
            INFO = 5
         ELSE IF( MB.NE.NB ) THEN
            INFO = 7
         END IF
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( MB.LT.1 ) THEN
         INFO = 6
      ELSE IF( NB.LT.1 ) THEN
         INFO = 7
      ELSE IF( LDA.LT.0 ) THEN
         INFO = 9
      ELSE IF( ( IAROW.LT.0 ).OR.( IAROW.GE.NPROW ) ) THEN
         INFO = 10
      ELSE IF( ( IACOL.LT.0 ).OR.( IACOL.GE.NPCOL ) ) THEN
         INFO = 11
      ELSE IF( MOD(IROFF,MB).GT.0 ) THEN
         INFO = 13
      ELSE IF( IRNUM.GT.(MP-IROFF) ) THEN
         INFO = 14
      ELSE IF( MOD(ICOFF,NB).GT.0 ) THEN
         INFO = 15
      ELSE IF( ICNUM.GT.(NQ-ICOFF) ) THEN
         INFO = 16
      ELSE IF( ( MYROW.LT.0 ).OR.( MYROW.GE.NPROW ) ) THEN
         INFO = 17
      ELSE IF( ( MYCOL.LT.0 ).OR.( MYCOL.GE.NPCOL ) ) THEN
         INFO = 18
      END IF
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDMATGEN', INFO )
         RETURN
      END IF
*
      MRROW = MOD( NPROW+MYROW-IAROW, NPROW )
      MRCOL = MOD( NPCOL+MYCOL-IACOL, NPCOL )
      NPMB  = NPROW * MB
      NQNB  = NPCOL * NB
      MOFF  = IROFF / MB
      NOFF  = ICOFF / NB
      MEND  = ICEIL(IRNUM, MB) + MOFF
      NEND  = ICEIL(ICNUM, NB) + NOFF
*
      MULT(1)  = MULT0
      MULT(2)  = MULT1
      IADD(1)  = IADD0
      IADD(2)  = IADD1
      JSEED(1) = ISEED
      JSEED(2) = 0
*
*     Symmetric or Hermitian matrix will be generated.
*
      IF( SYMM.OR.HERM ) THEN
*
*        First, generate the lower triangular part (with diagonal block)
*
         JUMP1 = 1
         JUMP2 = NPMB
         JUMP3 = M
         JUMP4 = NQNB
         JUMP5 = NB
         JUMP6 = MRCOL
         JUMP7 = MB*MRROW
*
         CALL XJUMPM( JUMP1, MULT, IADD, JSEED, IRAN1, IA1,   IC1 )
         CALL XJUMPM( JUMP2, MULT, IADD, IRAN1, ITMP1, IA2,   IC2 )
         CALL XJUMPM( JUMP3, MULT, IADD, IRAN1, ITMP1, IA3,   IC3 )
         CALL XJUMPM( JUMP4, IA3,  IC3,  IRAN1, ITMP1, IA4,   IC4 )
         CALL XJUMPM( JUMP5, IA3,  IC3,  IRAN1, ITMP1, IA5,   IC5 )
         CALL XJUMPM( JUMP6, IA5,  IC5,  IRAN1, ITMP3, ITMP1, ITMP2 )
         CALL XJUMPM( JUMP7, MULT, IADD, ITMP3, IRAN1, ITMP1, ITMP2 )
         CALL XJUMPM( NOFF,  IA4,  IC4,  IRAN1, ITMP1, ITMP2, ITMP3 )
         CALL XJUMPM( MOFF,  IA2,  IC2,  ITMP1, IRAN1, ITMP2, ITMP3 )
         CALL SETRAN( IRAN1, IA1,  IC1 )
*
         DO 10 I = 1, 2
            IB1(I) = IRAN1(I)
            IB2(I) = IRAN1(I)
            IB3(I) = IRAN1(I)
   10    CONTINUE
*
         JK = 1
         DO 80 IC = NOFF+1, NEND
            IOFFC = ((IC-1)*NPCOL+MRCOL) * NB
            DO 70 I = 1, NB
               IF( JK .GT. ICNUM ) GO TO 90
*
               IK = 1
               DO 50 IR = MOFF+1, MEND
                  IOFFR = ((IR-1)*NPROW+MRROW) * MB
*
                  IF( IOFFR .GT. IOFFC ) THEN
                     DO 20 J = 1, MB
                        IF( IK .GT. IRNUM ) GO TO 60
                        A(IK,JK) = ONE - TWO*PDRAND(0)
                        IK = IK + 1
   20                CONTINUE
*
                  ELSE IF( IOFFC .EQ. IOFFR ) THEN
                     IK = IK + I - 1
                     IF( IK .GT. IRNUM ) GO TO 60
                     DO 30 J = 1, I-1
                        A(IK,JK) = ONE - TWO*PDRAND(0)
   30                CONTINUE
                     A(IK,JK) = ONE - TWO*PDRAND(0)
                     DO 40 J = 1, MB-I
                        IF( IK+J .GT. IRNUM ) GO TO 60
                        A(IK+J,JK) = ONE - TWO*PDRAND(0)
                        A(IK,JK+J) = A(IK+J,JK)
   40                CONTINUE
                     IK = IK + MB - I + 1
                  ELSE
                     IK = IK + MB
                  END IF
*
                  CALL JUMPIT( IA2, IC2, IB1, IRAN2 )
                  IB1(1) = IRAN2(1)
                  IB1(2) = IRAN2(2)
   50          CONTINUE
*
   60          CONTINUE
               JK = JK + 1
               CALL JUMPIT( IA3, IC3, IB2, IRAN3 )
               IB1(1) = IRAN3(1)
               IB1(2) = IRAN3(2)
               IB2(1) = IRAN3(1)
               IB2(2) = IRAN3(2)
   70       CONTINUE
*
            CALL JUMPIT( IA4, IC4, IB3, IRAN4 )
            IB1(1) = IRAN4(1)
            IB1(2) = IRAN4(2)
            IB2(1) = IRAN4(1)
            IB2(2) = IRAN4(2)
            IB3(1) = IRAN4(1)
            IB3(2) = IRAN4(2)
   80    CONTINUE
*
*        Next, generate the upper triangular part.
*
   90    CONTINUE
         MULT(1)  = MULT0
         MULT(2)  = MULT1
         IADD(1)  = IADD0
         IADD(2)  = IADD1
         JSEED(1) = ISEED
         JSEED(2) = 0
*
         JUMP1 = 1
         JUMP2 = NQNB
         JUMP3 = N
         JUMP4 = NPMB
         JUMP5 = MB
         JUMP6 = MRROW
         JUMP7 = NB*MRCOL
*
         CALL XJUMPM( JUMP1, MULT, IADD, JSEED, IRAN1, IA1,   IC1 )
         CALL XJUMPM( JUMP2, MULT, IADD, IRAN1, ITMP1, IA2,   IC2 )
         CALL XJUMPM( JUMP3, MULT, IADD, IRAN1, ITMP1, IA3,   IC3 )
         CALL XJUMPM( JUMP4, IA3,  IC3,  IRAN1, ITMP1, IA4,   IC4 )
         CALL XJUMPM( JUMP5, IA3,  IC3,  IRAN1, ITMP1, IA5,   IC5 )
         CALL XJUMPM( JUMP6, IA5,  IC5,  IRAN1, ITMP3, ITMP1, ITMP2 )
         CALL XJUMPM( JUMP7, MULT, IADD, ITMP3, IRAN1, ITMP1, ITMP2 )
         CALL XJUMPM( MOFF,  IA4,  IC4,  IRAN1, ITMP1, ITMP2, ITMP3 )
         CALL XJUMPM( NOFF,  IA2,  IC2,  ITMP1, IRAN1, ITMP2, ITMP3 )
         CALL SETRAN( IRAN1, IA1,  IC1 )
*
         DO 100 I = 1, 2
            IB1(I) = IRAN1(I)
            IB2(I) = IRAN1(I)
            IB3(I) = IRAN1(I)
  100    CONTINUE
*
         IK = 1
         DO 150 IR = MOFF+1, MEND
            IOFFR = ((IR-1)*NPROW+MRROW) * MB
            DO 140 J = 1, MB
               IF( IK .GT. IRNUM ) GO TO 160
               JK = 1
               DO 120 IC = NOFF+1, NEND
                  IOFFC = ((IC-1)*NPCOL+MRCOL) * NB
                  IF( IOFFC .GT. IOFFR ) THEN
                     DO 110 I = 1, NB
                        IF( JK .GT. ICNUM ) GO TO 130
                        A(IK,JK) = ONE - TWO*PDRAND(0)
                        JK = JK + 1
  110                CONTINUE
                  ELSE
                     JK = JK + NB
                  END IF
                  CALL JUMPIT( IA2, IC2, IB1, IRAN2 )
                  IB1(1) = IRAN2(1)
                  IB1(2) = IRAN2(2)
  120          CONTINUE
*
  130          CONTINUE
               IK = IK + 1
               CALL JUMPIT( IA3, IC3, IB2, IRAN3 )
               IB1(1) = IRAN3(1)
               IB1(2) = IRAN3(2)
               IB2(1) = IRAN3(1)
               IB2(2) = IRAN3(2)
  140       CONTINUE
*
            CALL JUMPIT( IA4, IC4, IB3, IRAN4 )
            IB1(1) = IRAN4(1)
            IB1(2) = IRAN4(2)
            IB2(1) = IRAN4(1)
            IB2(2) = IRAN4(2)
            IB3(1) = IRAN4(1)
            IB3(2) = IRAN4(2)
  150    CONTINUE
  160    CONTINUE
*
*     (Conjugate) Transposed matrix A will be generated.
*
      ELSE IF( TRAN .OR. LSAME( AFORM, 'C' ) ) THEN
*
         JUMP1 = 1
         JUMP2 = NQNB
         JUMP3 = N
         JUMP4 = NPMB
         JUMP5 = MB
         JUMP6 = MRROW
         JUMP7 = NB*MRCOL
*
         CALL XJUMPM( JUMP1, MULT, IADD, JSEED, IRAN1, IA1,   IC1 )
         CALL XJUMPM( JUMP2, MULT, IADD, IRAN1, ITMP1, IA2,   IC2 )
         CALL XJUMPM( JUMP3, MULT, IADD, IRAN1, ITMP1, IA3,   IC3 )
         CALL XJUMPM( JUMP4, IA3,  IC3,  IRAN1, ITMP1, IA4,   IC4 )
         CALL XJUMPM( JUMP5, IA3,  IC3,  IRAN1, ITMP1, IA5,   IC5 )
         CALL XJUMPM( JUMP6, IA5,  IC5,  IRAN1, ITMP3, ITMP1, ITMP2 )
         CALL XJUMPM( JUMP7, MULT, IADD, ITMP3, IRAN1, ITMP1, ITMP2 )
         CALL XJUMPM( MOFF,  IA4,  IC4,  IRAN1, ITMP1, ITMP2, ITMP3 )
         CALL XJUMPM( NOFF,  IA2,  IC2,  ITMP1, IRAN1, ITMP2, ITMP3 )
         CALL SETRAN( IRAN1, IA1,  IC1 )
*
         DO 170 I = 1, 2
            IB1(I) = IRAN1(I)
            IB2(I) = IRAN1(I)
            IB3(I) = IRAN1(I)
  170    CONTINUE
*
         IK = 1
         DO 220 IR = MOFF+1, MEND
            IOFFR = ((IR-1)*NPROW+MRROW) * MB
            DO 210 J = 1, MB
               IF( IK .GT. IRNUM ) GO TO 230
               JK = 1
               DO 190 IC = NOFF+1, NEND
                  IOFFC = ((IC-1)*NPCOL+MRCOL) * NB
                  DO 180 I = 1, NB
                     IF( JK .GT. ICNUM ) GO TO 200
                     A(IK,JK) = ONE - TWO*PDRAND(0)
                     JK = JK + 1
  180             CONTINUE
                  CALL JUMPIT( IA2, IC2, IB1, IRAN2 )
                  IB1(1) = IRAN2(1)
                  IB1(2) = IRAN2(2)
  190          CONTINUE
*
  200          CONTINUE
               IK = IK + 1
               CALL JUMPIT( IA3, IC3, IB2, IRAN3 )
               IB1(1) = IRAN3(1)
               IB1(2) = IRAN3(2)
               IB2(1) = IRAN3(1)
               IB2(2) = IRAN3(2)
  210       CONTINUE
*
            CALL JUMPIT( IA4, IC4, IB3, IRAN4 )
            IB1(1) = IRAN4(1)
            IB1(2) = IRAN4(2)
            IB2(1) = IRAN4(1)
            IB2(2) = IRAN4(2)
            IB3(1) = IRAN4(1)
            IB3(2) = IRAN4(2)
  220    CONTINUE
  230    CONTINUE
*
*     A random matrix is generated.
*
      ELSE
*
         JUMP1 = 1
         JUMP2 = NPMB
         JUMP3 = M
         JUMP4 = NQNB
         JUMP5 = NB
         JUMP6 = MRCOL
         JUMP7 = MB*MRROW
*
         CALL XJUMPM( JUMP1, MULT, IADD, JSEED, IRAN1, IA1,   IC1 )
         CALL XJUMPM( JUMP2, MULT, IADD, IRAN1, ITMP1, IA2,   IC2 )
         CALL XJUMPM( JUMP3, MULT, IADD, IRAN1, ITMP1, IA3,   IC3 )
         CALL XJUMPM( JUMP4, IA3,  IC3,  IRAN1, ITMP1, IA4,   IC4 )
         CALL XJUMPM( JUMP5, IA3,  IC3,  IRAN1, ITMP1, IA5,   IC5 )
         CALL XJUMPM( JUMP6, IA5,  IC5,  IRAN1, ITMP3, ITMP1, ITMP2 )
         CALL XJUMPM( JUMP7, MULT, IADD, ITMP3, IRAN1, ITMP1, ITMP2 )
         CALL XJUMPM( NOFF,  IA4,  IC4,  IRAN1, ITMP1, ITMP2, ITMP3 )
         CALL XJUMPM( MOFF,  IA2,  IC2,  ITMP1, IRAN1, ITMP2, ITMP3 )
         CALL SETRAN( IRAN1, IA1,  IC1 )
*
         DO 240 I = 1, 2
            IB1(I) = IRAN1(I)
            IB2(I) = IRAN1(I)
            IB3(I) = IRAN1(I)
  240    CONTINUE
*
         JK = 1
         DO 290 IC = NOFF+1, NEND
            IOFFC = ((IC-1)*NPCOL+MRCOL) * NB
            DO 280 I = 1, NB
               IF( JK .GT. ICNUM ) GO TO 300
               IK = 1
               DO 260 IR = MOFF+1, MEND
                  IOFFR = ((IR-1)*NPROW+MRROW) * MB
                  DO 250 J = 1, MB
                     IF( IK .GT. IRNUM ) GO TO 270
                     A(IK,JK) = ONE - TWO*PDRAND(0)
                     IK = IK + 1
  250             CONTINUE
                  CALL JUMPIT( IA2, IC2, IB1, IRAN2 )
                  IB1(1) = IRAN2(1)
                  IB1(2) = IRAN2(2)
  260          CONTINUE
*
  270          CONTINUE
               JK = JK + 1
               CALL JUMPIT( IA3, IC3, IB2, IRAN3 )
               IB1(1) = IRAN3(1)
               IB1(2) = IRAN3(2)
               IB2(1) = IRAN3(1)
               IB2(2) = IRAN3(2)
  280       CONTINUE
*
            CALL JUMPIT( IA4, IC4, IB3, IRAN4 )
            IB1(1) = IRAN4(1)
            IB1(2) = IRAN4(2)
            IB2(1) = IRAN4(1)
            IB2(2) = IRAN4(2)
            IB3(1) = IRAN4(1)
            IB3(2) = IRAN4(2)
  290    CONTINUE
  300    CONTINUE
      END IF
*
*     Diagonally dominant matrix will be generated.
*
      IF( LSAME( DIAG, 'D' ) ) THEN
         IF( MB.NE.NB ) THEN
            WRITE(*,*) 'Diagonally dominant matrices with rowNB not'//
     $                 ' equal colNB is not supported!'
            RETURN
         END IF
*
         MAXMN = MAX(M, N)
         JK    = 1
         DO 340 IC = NOFF+1, NEND
            IOFFC = ((IC-1)*NPCOL+MRCOL) * NB
            IK    = 1
            DO 320 IR = MOFF+1, MEND
               IOFFR = ((IR-1)*NPROW+MRROW) * MB
               IF( IOFFC.EQ.IOFFR ) THEN
                  DO 310 J = 0, MB-1
                     IF( IK .GT. IRNUM ) GO TO 330
                     A(IK,JK+J) = ABS(A(IK,JK+J)) + MAXMN
                     IK = IK + 1
  310             CONTINUE
               ELSE
                  IK = IK + MB
               END IF
  320       CONTINUE
  330       CONTINUE
            JK = JK + NB
  340    CONTINUE
      END IF
*
      RETURN
*
*     End of PDMATGEN
*
      END
      SUBROUTINE PDMYROT( TYPE, A, DESCA, LDA, NB, IROW1, ICOL1, ISTART,
     $                 ISTOP, ITMP1, ITMP2,
     $                 CS, TMPA, LDTMP)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            NB, LDTMP
      INTEGER            ICOL1, IROW1, ISTART, ISTOP, ITMP1, ITMP2, LDA

*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), CS( * )
      DOUBLE PRECISION   TMPA(LDTMP, LDTMP)
	INTEGER			   DESCA( * ) 
*     ..
*
*  Purpose
*  =======
*
*  PDMYROT applies one or several Givens rotations to a matrix
*     on either their rows or columns.
*
*  Arguments
*  =========
*
*  TYPE    (global input) CHARACTER*1
*          If 'R': Apply reflectors to the rows of the matrix
*              (apply from left)
*          Otherwise: Apply reflectors to the columns of the matrix
*          Unchanged on exit.
*
*  A       (global input/output) DOUBLE PRECISION array, (LDA,*)
*          On entry, the matrix to receive the reflections.
*          The updated matrix on exit.
*
*  LDA     (local input) INTEGER
*          On entry, the leading dimension of A.  Unchanged on exit.
*
*  IROW1   (local input/output) INTEGER
*          On entry, the local row element of A.
*          Undefined on output.
*
*
*  ICOL1   (local input/output) INTEGER
*          On entry, the local column element of A.
*          Undefined on output.
*
*  ISTART  (global input) INTEGER
*          Specifies the "number" of the first cspair.  This is
*              used as an index into VECS.
*
*  ISTOP   (global input) INTEGER
*          Specifies the "number" of the last reflector.
*
*  ITMP1   (local input) INTEGER
*          Starting range into A.  For rows, this is the local
*              first column.  For columns, this is the local first row.
*
*  ITMP2   (local input) INTEGER
*          Ending range into A.  For rows, this is the local last
*              column.  For columns, this is the local last row.
*
*  CS    (global input) DOUBLE PRECISION array of size 2*N (matrix
*                                                             size)
*          This holds the C and S one after another and this
*
*  =====================================================================
*     ..
*     .. Parameters ..
      INTEGER            DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
      INTEGER            ROWWISE, COLUMNWISE, ALL
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2, ALL=3)
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )


*     ..
*     .. Local Arrays


*     ..
*     .. Local Scalars ..
      INTEGER            J, K, JC, JR, CNT, JJC, JJR
      INTEGER            IAROW1, IACOL1, LROW1, LCOL1
      INTEGER            IAROW2, IACOL2, LROW2, LCOL2
      INTEGER            MYROW, MYCOL, NPROW, NPCOL, ICTXT

      DOUBLE PRECISION   TEMP, C, S
*     ..
*     .. External Functions ..
      INTEGER            INDXG2P, INDXG2L
      EXTERNAL           INDXG2P, INDXG2L

*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     
	ICTXT = DESCA(CTXT_)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL ) 

      IF (( ISTART .GT. ISTOP ) .OR. (ITMP1 .GT. ITMP2 )) THEN
         RETURN
      END IF
      IF( TYPE.EQ.'R' )  THEN
	   DO 6 K = ISTART, ISTOP
	     J = IROW1
           IAROW1 = INDXG2P(J, NB, 0, 0, NPROW)
           IAROW2 = INDXG2P(J+1, NB, 0, 0, NPROW)
           LROW1 = INDXG2L(J, NB, 0, 0, NPROW)
           LROW2 = INDXG2L(J+1, NB, 0, 0, NPROW)
           C = CS( ( K-1 )*2+1 )
           S = CS( ( K-1 )*2+2 )
		 DO 5 JC = ITMP1, ITMP1+(NB-MOD(ITMP1-1, NB))-1, NB
               CNT = MIN(NB-MOD(ITMP1-1, NB), ITMP2-JC+1)
			 LCOL1 = INDXG2L(JC, NB, 0, 0, NPCOL)
			 IACOL1 = INDXG2P(JC, NB, 0, 0, NPCOL)
            
               IF (IAROW1.NE.IAROW2) THEN
                  IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL1)) THEN
                     Call Send(1, CNT, DOWN, 
     $                    A(LROW1, LCOL1), LDA,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                     Call Recv(1, CNT, DOWN, 
     $                    TMPA, LDTMP,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                  ELSE IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.IACOL1)) 
     $                    THEN
                     Call Send(1, CNT, UP, 
     $                    A(LROW2, LCOL1), LDA,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                     Call Recv(1, CNT, UP, 
     $                    TMPA, LDTMP,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                  END IF                  
               END IF
               DO 7 JJC = JC, JC+CNT-1
                  LCOL1 = INDXG2L(JJC, NB, 0, 0, NPCOL)
                  IF (IAROW1.EQ.IAROW2) THEN
                     IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL1)) THEN
                        TEMP = C*A( LROW1, LCOL1 ) + S*A( LROW2, LCOL1 )
                        A( LROW2, LCOL1 ) = -S*A( LROW1, LCOL1 ) + 
     $                       C*A( LROW2, LCOL1 )
                        A( LROW1, LCOL1 ) = TEMP
                     END IF
                  ELSE
                     IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL1)) THEN
                        A( LROW1, LCOL1 ) = C*A( LROW1, LCOL1 ) + 
     $                       S*TMPA( 1, JJC-JC+1)
                     ELSE IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.IACOL1))
     $                       THEN
                        A( LROW2, LCOL1 ) = -S*TMPA( 1, JJC-JC+1) + 
     $                       C*A( LROW2, LCOL1 )
                     END IF
                  END IF
 7             CONTINUE
               
 5		  CONTINUE

		  DO 10 JC = ITMP1+(NB-MOD(ITMP1-1, NB)), ITMP2, NB
			 CNT = MIN(NB, ITMP2-JC+1)
               IACOL1 = INDXG2P(JC, NB, 0, 0, NPCOL)
               LCOL1 = INDXG2L(JC, NB, 0, 0, NPCOL)               

               IF (IAROW1.NE.IAROW2) THEN
                  IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL1)) THEN
                     Call Send(1, CNT, DOWN, 
     $                    A(LROW1, LCOL1), LDA,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                     Call Recv(1, CNT, DOWN, 
     $                    TMPA, LDTMP,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                  ELSE IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.IACOL1)) 
     $                    THEN
                     Call Send(1, CNT, UP, 
     $                    A(LROW2, LCOL1), LDA,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                     Call Recv(1, CNT, UP, 
     $                    TMPA, LDTMP,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                  END IF                  
               END IF

               DO 12 JJC = JC, JC+CNT-1
                  LCOL1 = INDXG2L(JJC, NB, 0, 0, NPCOL)
                  IF (IAROW1.EQ.IAROW2) THEN
                     IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL1)) THEN
                        TEMP = C*A( LROW1, LCOL1 ) + S*A( LROW2, LCOL1 )
                        A( LROW2, LCOL1 ) = -S*A( LROW1, LCOL1 ) + 
     $                       C*A( LROW2, LCOL1 )
                        A( LROW1, LCOL1 ) = TEMP
                     END IF
                  ELSE
                     IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL1)) THEN
                        A( LROW1, LCOL1 ) = C*A( LROW1, LCOL1 ) + 
     $                       S*TMPA( 1, JJC-JC+1)
                     ELSE IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.IACOL1))
     $                       THEN
                        A( LROW2, LCOL1 ) = -S*TMPA( 1, JJC-JC+1) + 
     $                       C*A( LROW2, LCOL1 )
                     END IF
                  END IF
 12            CONTINUE
 10         CONTINUE
		  J = J + 1
 6       CONTINUE
      ELSE IF( TYPE.EQ. 'C'  ) THEN
	   J = ICOL1
         LCOL1 = INDXG2L(J, NB, 0, 0, NPCOL)
         LCOL2 = INDXG2L(J+1, NB, 0, 0, NPCOL)
         IACOL1 = INDXG2P(J, NB, 0, 0, NPCOL)
         IACOL2 = INDXG2P(J+1, NB, 0, 0, NPCOL)

         DO 40 K = ISTART, ISTOP
            C = CS( ( K-1 )*2+1 )
            S = CS( ( K-1 )*2+2 )
            


            DO 41 JR = ITMP1,ITMP1+(NB-MOD(ITMP1-1, NB))-1, NB
         
			 CNT = MIN(NB-MOD(ITMP1-1, NB), ITMP2-JR+1)
	         
               IAROW1 = INDXG2P(JR, NB, 0, 0, NPROW)
               LROW1 = INDXG2L(JR, NB, 0, 0, NPROW)
               IF (IACOL1.NE.IACOL2) THEN
                  IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL1)) THEN
                     Call Send(CNT, 1, RIGHT, 
     $                    A(LROW1, LCOL1), LDA,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                     Call Recv(CNT, 1, RIGHT, 
     $                    TMPA, LDTMP,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                  ELSE IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL2)) 
     $                    THEN
                     Call Send(CNT, 1, LEFT, 
     $                    A(LROW1, LCOL2), LDA,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                     Call Recv(CNT, 1, LEFT, 
     $                    TMPA, LDTMP,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                  END IF                  
               END IF                  

               DO 42 JJR = JR, JR+CNT-1
                  LROW1 = INDXG2L(JJR, NB, 0, 0, NPROW)
                  IF (IACOL1.EQ.IACOL2) THEN
				   IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL1)) THEN 	
	                   TEMP = C*A( LROW1, LCOL2 ) + S*A( LROW1,LCOL1)
		               A( LROW1, LCOL1 ) = -S*A( LROW1, LCOL2 ) + 
     $			           C*A( LROW1, LCOL1 )
					   A( LROW1, LCOL2 ) = TEMP
	               END IF
                  ELSE
                     IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL1)) THEN
                        A( LROW1, LCOL1 ) = -S*TMPA( JJR-JR+1, 1 ) + 
     $                       C*A( LROW1, LCOL1 )
                     ELSE IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL2)) 
     $                       THEN
                        A( LROW1, LCOL2 ) = C*A( LROW1, LCOL2 ) + 
     $                       S*TMPA( JJR-JR+1, 1)
                     END IF
                  END IF
 42            CONTINUE
 41         CONTINUE


            DO 43 JR = ITMP1+(NB-MOD(ITMP1-1, NB)), ITMP2, NB
               CNT = MIN(NB, ITMP2-JR+1)
               IAROW1 = INDXG2P(JR, NB, 0, 0, NPROW)
               LROW1 = INDXG2L(JR, NB, 0, 0, NPROW)


               IF (IACOL1.NE.IACOL2) THEN
                  IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL1)) THEN
                     Call Send(CNT, 1, RIGHT, 
     $                    A(LROW1, LCOL1), LDA,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                     Call Recv(CNT, 1, RIGHT, 
     $                    TMPA, LDTMP,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                  ELSE IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL2)) THEN
                     Call Send(CNT, 1, LEFT, 
     $                    A(LROW1, LCOL2), LDA,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                     Call Recv(CNT, 1, LEFT, 
     $                    TMPA, LDTMP,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                  END IF                  
               END IF 

               DO 44 JJR = JR, JR+CNT-1
                  LROW1 = INDXG2L(JJR, NB, 0, 0, NPROW)
                  IF (IACOL1.EQ.IACOL2) THEN
					IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL1)) THEN 		
						TEMP = C*A( LROW1, LCOL2 ) + S*A( LROW1,LCOL1)
						A( LROW1, LCOL1 ) = -S*A( LROW1, LCOL2 ) + 
     $						C*A( LROW1, LCOL1 )
						A( LROW1, LCOL2 ) = TEMP
					END IF
                  ELSE
                     IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL1)) THEN
                        A( LROW1, LCOL1 ) = -S*TMPA( JJR-JR+1, 1 ) +
     $                       C*A( LROW1, LCOL1 )
                     ELSE IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL2)) 
     $                       THEN 
                        A( LROW1, LCOL2 ) = C*A( LROW1, LCOL2 ) + 
     $                       S*TMPA( JJR-JR+1, 1)
                     END IF
                  END IF

 44            CONTINUE
 43         CONTINUE
            J = J + 1
 40      CONTINUE
	ELSE
	   J = ICOL1
         LCOL1 = INDXG2L(J, NB, 0, 0, NPCOL)
         LCOL2 = INDXG2L(J+1, NB, 0, 0, NPCOL)
         IACOL1 = INDXG2P(J, NB, 0, 0, NPCOL)
         IACOL2 = INDXG2P(J+1, NB, 0, 0, NPCOL)
         DO 60 K = ISTART, ISTOP
            C = CS( ( K-1 )*2+1 )
            S = CS( ( K-1 )*2+2 )

            DO 61 JR = ITMP1, ITMP1+(NB-MOD(ITMP1-1, NB))-1, NB
               CNT = MIN(NB-MOD(ITMP1-1, NB), ITMP2-JR+1)

               IAROW1 = INDXG2P(JR, NB, 0, 0, NPROW)
               LROW1 = INDXG2L(JR, NB, 0, 0, NPROW)
               IF (IACOL1.NE.IACOL2) THEN
                  IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL1)) THEN
                     Call Send(CNT, 1, RIGHT, 
     $                    A(LROW1, LCOL1), LDA,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                     Call Recv(CNT, 1, RIGHT, 
     $                    TMPA, LDTMP,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                  ELSE IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL2))
     $                    THEN
                     Call Send(CNT, 1, LEFT, 
     $                    A(LROW1, LCOL2), LDA,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                     Call Recv(CNT, 1, LEFT, 
     $                    TMPA, LDTMP,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                  END IF                  
               END IF                  

               DO 62 JJR = JR, JR+CNT-1
                  LROW1 = INDXG2L(JJR, NB, 0, 0, NPROW)
                  IF (IACOL1.EQ.IACOL2) THEN
				   IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL1)) THEN 	 	
						TEMP = C*A( LROW1, LCOL1 ) + S*A( LROW1,LCOL2)
						A( LROW1, LCOL2 ) = -S*A( LROW1, LCOL1 ) + 
     $						C*A( LROW1, LCOL2 )
						A( LROW1, LCOL1 ) = TEMP
				   END IF
                  ELSE
                     IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL1)) THEN
                        A( LROW1, LCOL1 ) = C*A( LROW1, LCOL1 ) + 
     $                       S*TMPA( JJR-JR+1, 1)
                      ELSE IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL2)) 
     $                       THEN
                         A( LROW1, LCOL2 ) = -S*TMPA( JJR-JR+1, 1 ) + 
     $                        C*A( LROW1, LCOL2 )
                     END IF
                  END IF
 62            CONTINUE
 61         CONTINUE

            DO 63 JR = ITMP1+(NB-MOD(ITMP1-1, NB)), ITMP2, NB
               CNT = MIN(NB, ITMP2-JR+1)
               IAROW1 = INDXG2P(JR, NB, 0, 0, NPROW)
               LROW1 = INDXG2L(JR, NB, 0, 0, NPROW)


               IF (IACOL1.NE.IACOL2) THEN
                  IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL1)) THEN
                     Call Send(CNT, 1, RIGHT, 
     $                    A(LROW1, LCOL1), LDA,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                     Call Recv(CNT, 1, RIGHT, 
     $                    TMPA, LDTMP,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                  ELSE IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL2)) 
     $                    THEN
                     Call Send(CNT, 1, LEFT, 
     $                    A(LROW1, LCOL2), LDA,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                     Call Recv(CNT, 1, LEFT, 
     $                    TMPA, LDTMP,
     $                    MYROW, MYCOL, NPCOL, NPROW, ICTXT)
                  END IF                  
               END IF 
               DO 64 JJR = JR, JR+CNT-1
                  LROW1 = INDXG2L(JJR, NB, 0, 0, NPROW)
                  IF (IACOL1.EQ.IACOL2) THEN
				   IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL1)) THEN 	
						TEMP = C*A( LROW1, LCOL1 ) + S*A( LROW1,LCOL2)
						A( LROW1, LCOL2 ) = -S*A( LROW1, LCOL1 ) + 
     $						C*A( LROW1, LCOL2 )
						A( LROW1, LCOL1 ) = TEMP
	               END IF
                  ELSE
                     IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL1)) THEN
                        A( LROW1, LCOL1 ) = C*A( LROW1, LCOL1 ) + 
     $                       S*TMPA( JJR-JR+1, 1 )
                     ELSE IF ((MYROW.EQ.IAROW1).AND.(MYCOL.EQ.IACOL2)) 
     $                       THEN
                        A( LROW1, LCOL2 ) = -S*TMPA( JJR-JR+1, 1 ) + 
     $                       C*A( LROW1, LCOL2 )
                     END IF
                  END IF

 64            CONTINUE
 63         CONTINUE
            J = J + 1
 60      CONTINUE
      END IF
      RETURN
*     ..
*     .. End of PDMYROT
*     ..
      END
      SUBROUTINE PQZ2( N, NB, A, LDA, DESCA, IFRST,
     $     B, LDB, DESCB, VL, VR, NCOL, 
     $	 TALPHAI, TALPHAR, TBETA,SAFMIN,
     $	 ESHIFT, TMPW, TMPW2, LDTMP,JITER,ATOL)

      IMPLICIT NONE
*
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, N, NCOL, NB, IFRST
      INTEGER            ESHIFT, LDTMP,JITER
	DOUBLE PRECISION   SAFMIN, ATOL
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), VL( * ), VR( * )
      DOUBLE PRECISION   TMPW(LDTMP, *), TMPW2(LDTMP, *)
	DOUBLE PRECISION   TALPHAI(*), TALPHAR(*), TBETA(*)
      INTEGER            DESCA(*), DESCB(*)

*  Purpose
*  =======
*
*     PQZ2 chases bulges N steps downwards the diagonal of A and B

*     ..
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      INTEGER            DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
      INTEGER            ROWWISE, COLUMNWISE, ALL
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2, ALL=3)

      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )

*     ..
*     .. Local Scalars ..
      LOGICAL            ILPIVT
      INTEGER            IACOL, IAROW, IACOL2, IAROW2, IAROW3
      INTEGER            LROW1, LROW2, LROW3
      INTEGER            LCOL1, LCOL2, LCOL3
      INTEGER            DCNT, RCNT
      INTEGER            J, JC, JR, JMP, LDV
	INTEGER			   ICTXT, MYROW, MYCOL, NPROW, NPCOL
      DOUBLE PRECISION   SCALE, T, VS, U2
      DOUBLE PRECISION   TAU, TEMP, TEMP2, U1, W11, W12, W21, W22
	DOUBLE PRECISION   D1, D2, D3, D4, D5
      INTEGER            RET, MAXJ, NUM1, NUM2
	INTEGER			   INFO

*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   V(4)
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLANHS, DLAPY2, DLAPY3
      EXTERNAL           DLANHS, DLAPY2, DLAPY3
      EXTERNAL           INDXG2P, INDXG2L
      INTEGER            INDXG2P, INDXG2L
      
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFG, PDLACP3

*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Machine Constants
*

	ICTXT = DESCA(CTXT_)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )  
	LDV = 4
      JMP = 0

	
	MAXJ = IFRST+1 + (N - NCOL) -2
      DO 90 J = IFRST+1, MAXJ
*
*        Use 3x3 Householder transforms.
*
*        Zero (j-1)st row of A
*
	   NUM1 = MIN(IFRST+N-1,J+2)	
	   NUM2 = MIN(IFRST+N-1,J+3)	
         IAROW = INDXG2P(J, NB, 0, 0, NPROW)
         IACOL = INDXG2P(J, NB, 0, 0, NPCOL)
	   IACOL2 = INDXG2P(J-1, NB, 0, 0, NPCOL)
	   IF (INDXG2P(J+2, NB, 0, 0, NPROW).EQ.IAROW) THEN
		DCNT = 0
	   ELSE IF (INDXG2P(J+1, NB, 0, 0, NPROW).EQ.IAROW) THEN
		DCNT = 1
	   ELSE 
		DCNT = 2
	   END IF

	   

	   CALL PDLACP3(4, J-1, A, DESCA, TMPW, LDTMP, IAROW, IACOL, 0)
	   CALL PDLACP3(3, J, B, DESCB, TMPW2, LDTMP, IAROW, IACOL, 0)
*        Calculate V
*	   ..
*	   .. Might have to restart bump if there was a vigilant deflation
*	   .. No check here    
         IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN   
			IF ((ABS(TMPW(2,1)).LT.ATOL).AND.
     $			(ABS(TMPW(3,1)).LT.ATOL).AND.
     $			(ABS(TMPW(4,1)).LT.ATOL)) THEN
  			

			   TMPW(2,1) = ZERO
			   TMPW(3,1) = ZERO
			   TMPW(4,1) = ZERO	  
	
	  		   CALL CREATEINITBUMP(TMPW(2,2), TMPW2, 
     $		  		ESHIFT, NB, 
     $				TALPHAR, TALPHAI, TBETA, V,LDTMP)
			   TAU = V(4)		    
			ELSE
		   V( 1 ) = TMPW( 2, 1 )
		   V( 2 ) = TMPW( 3, 1 )
		   V( 3 ) = TMPW( 4, 1 )
		   CALL DLARFG( 3, TMPW( 2, 1 ), V( 2 ), 1, TAU )
		   V( 1 ) = ONE           
		   V( 4 ) = TAU
		   TMPW(3,1) = ZERO
		   TMPW(4,1) = ZERO
	END IF
	   END IF		 


	   CALL PDLACP3(4, J-1, A, DESCA, TMPW, LDTMP, IAROW, IACOL, -1)


         

         IACOL2 = INDXG2P(IFRST+N-1, NB, 0, 0, NPCOL)
         IAROW2 = INDXG2P(J+2, NB, 0, 0, NPROW)

		
*        Communicate V to all processors needing it
	   IF (DCNT.GT.0) THEN
		 IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
	        V(1) = V(2)
			V(2) = V(3)
	        V(3) = TAU
			CALL Send(3, 1, DOWN, 
     $			V, LDV, MYROW, MYCOL, NPCOL, NPROW, ICTXT)
		 ELSE IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.IACOL)) THEN
			CALL Recv(3, 1, UP, 
     $	        V, LDV, MYROW, MYCOL, NPCOL, NPROW, ICTXT)
		 END IF
	   ELSE
		 IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
	        V(1) = V(2)
			V(2) = V(3)
	        V(3) = TAU
		 END IF
	   END IF
	   IF (IACOL.NE.IACOL2) THEN
		 IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
			CALL Send(3, 1, RIGHT, 
     $			V, LDV, MYROW, MYCOL, NPCOL, NPROW, ICTXT)
		 ELSE IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL2)) THEN
			CALL Recv(3, 1, LEFT, 
     $	        V, LDV, MYROW, MYCOL, NPCOL, NPROW, ICTXT)
		 END IF
		 IF (DCNT.GT.0) THEN
			IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.IACOL)) THEN
				CALL Send(3, 1, RIGHT, 
     $		       V, LDV, MYROW, MYCOL, NPCOL, NPROW, ICTXT)
			ELSE IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.IACOL2)) THEN
				CALL Recv(3, 1, LEFT, 
     $	           V, LDV, MYROW, MYCOL, NPCOL, NPROW, ICTXT)
			END IF
		 END IF

	   END IF

21	   CONTINUE         
*
*        Save Householder transforms.
*   
         VL( JMP + 1) = V( 1 )
         VL( JMP + 2) = V( 2 )
         VL( JMP + 3) = V( 3 )


*	   Update A and B from the left
	   CALL PDLAREF('Rows', A, LDA,
     $        N, NB, DESCA, 
     $       .TRUE., J, J, 1, 1, 
     $        J, IFRST+N-1,
     $        V, D1, D2, D3, D4, D5,
     $        TMPW, RET,LDTMP)

         CALL PDLAREF('Rows', B, LDB,
     $        N, NB, DESCB, 
     $       .TRUE., J, J, 1, 1, 
     $        J, IFRST+N-1,
     $        V, D1, D2, D3, D4, D5,
     $        TMPW, RET,LDTMP)

	  
*
*       Zero j-th column of B (see DLAGBC for details)
*


        IF (INDXG2P(J+2, NB, 0, 0, NPCOL).EQ.IACOL) THEN
           RCNT = 0
        ELSE IF (INDXG2P(J+1, NB, 0, 0, NPCOL).EQ.IACOL) THEN
           RCNT = 1
        ELSE 
           RCNT = 2
        END IF

        CALL PDLACP3(3, J, B, DESCB, TMPW, LDTMP, IAROW, IACOL, 0)


*       Calculate V
        IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
		 CALL INVHSE(3, TMPW, TMPW, LDTMP, TAU,V,INFO)
           V( 4 ) = TAU
        END IF



        IACOL2 = INDXG2P(J+2, NB, 0, 0, NPCOL)
	  IAROW2 = INDXG2P(IFRST, NB, 0, 0, NPROW)
	  IAROW3 = INDXG2P(NUM2, NB, 0, 0, NPROW)
	  


*        Communicate V to all processors needing it
	   IF (RCNT.GT.0) THEN
		 IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
      	    V(1) = V(2)
	        V(2) = V(3)
	        V(3) = TAU
			CALL Send(3, 1, RIGHT, 
     $			V, LDV, MYROW, MYCOL, NPCOL, NPROW, ICTXT)
		 ELSE IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL2)) THEN
			CALL Recv(3, 1, LEFT, 
     $	        V, LDV, MYROW, MYCOL, NPCOL, NPROW, ICTXT)
		 END IF
	   ELSE
		 IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
      	    V(1) = V(2)
	        V(2) = V(3)
	        V(3) = TAU
		END IF
	   END IF

	   IF ((IAROW.NE.IAROW2).OR.(IAROW.NE.IAROW3)) THEN
		 IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
			IF (MYROW.NE.IAROW2) THEN
				CALL Send(3, 1, UP, 
     $				V, LDV, MYROW, MYCOL, 
     $				NPCOL, NPROW, ICTXT)
			ELSE
				CALL Send(3, 1, DOWN, 
     $				V, LDV, MYROW, MYCOL, 
     $				NPCOL, NPROW, ICTXT)
			END IF
		 ELSE IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.IACOL)) THEN
			CALL Recv(3, 1, DOWN, 
     $	        V, LDV, MYROW, MYCOL, 
     $			NPCOL, NPROW, ICTXT)
		 ELSE IF ((MYROW.EQ.IAROW3).AND.(MYCOL.EQ.IACOL)) THEN
			CALL Recv(3, 1, UP, 
     $	        V, LDV, MYROW, MYCOL, 
     $			NPCOL, NPROW, ICTXT)
		 END IF
	     IF (RCNT.GT.0) THEN 
			IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL2)) THEN
				IF (MYROW.NE.IAROW2) THEN
					CALL Send(3, 1, UP, 
     $					V, LDV, MYROW, MYCOL,	
     $					NPCOL, NPROW, ICTXT)
				ELSE
					CALL Send(3, 1, DOWN, 
     $					V, LDV, MYROW, MYCOL, 
     $					NPCOL, NPROW, ICTXT)
				END IF
			ELSE IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.IACOL2)) THEN
				CALL Recv(3, 1, DOWN, 
     $	            V, LDV, MYROW, MYCOL, 
     $				NPCOL, NPROW, ICTXT)
			ELSE IF ((MYROW.EQ.IAROW3).AND.(MYCOL.EQ.IACOL2)) THEN
				CALL Recv(3, 1, UP, 
     $				V, LDV, MYROW, MYCOL, 
     $				NPCOL, NPROW, ICTXT)			
			END IF
		 END IF
	   END IF

 22	   CONTINUE
*	   ..
*        .. Save Householder transforms.
*	   ..  
         VR( JMP + 1) = V( 1 )
         VR( JMP + 2) = V( 2 )
         VR( JMP + 3) = V( 3 )

*	   Update A and B from the right
         CALL PDLAREF('Cols', A, LDA,
     $        N, NB, DESCA, 
     $       .TRUE., J, J, 1, 1, 
     $        IFRST, NUM2,
     $        V, D1, D2, D3, D4, D5,
     $        TMPW, RET,LDTMP)

         CALL PDLAREF('Cols', B, LDB,
     $        N, NB, DESCB, 
     $       .TRUE., J, J, 1, 1, 
     $        IFRST, NUM1,
     $        V, D1, D2, D3, D4, D5,
     $        TMPW, RET,LDTMP)

        IAROW = INDXG2P(J+1, NB, 0, 0, NPROW)
        IAROW2 = INDXG2P(J+2, NB, 0, 0, NPROW)

	  LCOL1 = INDXG2L(J, NB, 0, 0, NPCOL)
*       Zero out elements of B, column j
        IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
           LROW1 = INDXG2L(J+1, NB, 0, 0, NPROW)
           B(LROW1, LCOL1) = ZERO
        END IF
        IF ((MYROW.EQ.IAROW2).AND.(MYCOL.EQ.IACOL)) THEN
           LROW1 = INDXG2L(J+2, NB, 0, 0, NPROW)
           B(LROW1, LCOL1) = ZERO
        END IF
        JMP = JMP + 3
 90   CONTINUE
 100	CONTINUE
      RETURN
*
*     End of PQZ2
*
800   FORMAT(20F7.4)
      END

*	 HAVE TO ADD Support for vagilenta deflationer

      SUBROUTINE QZ2Local( N, A, LDA, B, LDB, VL, VR, NCOL,GJ,WSIZE,
     $			ATOL,NB,ESHIFT,TALPHAR,TALPHAI, TBETA)
	  IMPLICIT NONE
*
*
*     .. Scalar Arguments ..
      INTEGER             LDA, LDB, N, NCOL,GJ,WSIZE,ESHIFT,NB
	DOUBLE PRECISION	ATOL
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), VL( * ), VR( * )
	DOUBLE PRECISION   TALPHAI(*), TALPHAR(*), TBETA(*)
*     ..
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILPIVT
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   V( 4 )
*     ..
*     .. External Functions ..
 
      DOUBLE PRECISION   DLAMCH, DLANHS, DLAPY2, DLAPY3
	  DOUBLE PRECISION   TIME1, TIME2, TIME3, TM0, TM1, TM2, TM3, TMDUM
      EXTERNAL           DLAMCH, DLANHS, DLAPY2, DLAPY3
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAG2, DLARFG, DLARTG, DLASET, DLASV2, DROT,
     $                   XERBLA
	  INTEGER            J, JC, JR, JMP, INFO,STARTJ, NMAX
      DOUBLE PRECISION   SCALE, T, VS, U2, 
     $                   SAFMIN,SAFMAX, TAU,
     $                   TEMP, TEMP2, U1, W11, W12, W21, W22

*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Machine Constants
*
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
	  JMP = 0
      DO 290 J = 2, N - NCOL 
*
*        Use 3x3 Householder transforms.
*
*        Zero (j-1)st column of A
*
			IF ((ABS(A(J,J-1)).LT.ATOL).AND.
     $			(ABS(A(J+1,J-1)).LT.ATOL).AND.
     $			(ABS(A(J+2,J-1)).LT.ATOL)) THEN
  			

			   A(J,J-1) = ZERO
			   A(J+1,J-1) = ZERO
			   A(J+2,J-1) = ZERO	  
	
	  		   CALL CREATEINITBUMP(A(J,J), B(J,J), 
     $		  		ESHIFT, NB, 
     $				TALPHAR, TALPHAI, TBETA, V,LDA)			
			   TAU = V(4)		 
			ELSE

             V( 1 ) = A( J, J-1 )
             V( 2 ) = A( J+1, J-1 )
             V( 3 ) = A( J+2, J-1 )
*
             CALL DLARFG( 3, A( J, J-1 ), V( 2 ), 1, TAU )
             V( 1 ) = ONE
             A( J+1, J-1 ) = ZERO
             A( J+2, J-1 ) = ZERO
			ENDIF
*
*            Save Householder transforms.
*
			 VL( JMP + 1) = V( 2 )
			 VL( JMP + 2) = V( 3 )
			 VL( JMP + 3) = TAU 
*
		   NMAX = MAX(N,WSIZE)
             DO 230 JC = J, NMAX
                TEMP = TAU*( A( J, JC )+V( 2 )*A( J+1, JC )+V( 3 )*
     $                   A( J+2, JC ) )
                A( J, JC ) = A( J, JC ) - TEMP
                A( J+1, JC ) = A( J+1, JC ) - TEMP*V( 2 )
                A( J+2, JC ) = A( J+2, JC ) - TEMP*V( 3 )
                TEMP2 = TAU*( B( J, JC )+V( 2 )*B( J+1, JC )+V( 3 )*
     $                  B( J+2, JC ) )
                B( J, JC ) = B( J, JC ) - TEMP2
                B( J+1, JC ) = B( J+1, JC ) - TEMP2*V( 2 )
                B( J+2, JC ) = B( J+2, JC ) - TEMP2*V( 3 )
  230        CONTINUE
*
*            Zero j-th column of B
*
		 CALL INVHSE(3, B(J,J), B(J,J), LDB, TAU,V,INFO)


*
*             Save Householder transforms.
*
			  VR( JMP + 1) = V( 2 )
			  VR( JMP + 2) = V( 3 )
			  VR( JMP + 3) = TAU 
*
*             Apply transformations from the right.
*
			STARTJ = 1
*			IF (GJ.EQ.1) STARTJ = 2
              DO 260 JR = STARTJ, MIN( J+3, N )
                  TEMP = TAU*( A( JR, J )+V( 2 )*A( JR, J+1 )+V( 3 )*
     $                   A( JR, J+2 ) )
                  A( JR, J ) = A( JR, J ) - TEMP
                  A( JR, J+1 ) = A( JR, J+1 ) - TEMP*V( 2 )
                  A( JR, J+2 ) = A( JR, J+2 ) - TEMP*V( 3 )
  260         CONTINUE
              DO 270 JR = STARTJ, J + 2
                  TEMP = TAU*( B( JR, J )+V( 2 )*B( JR, J+1 )+V( 3 )*
     $                   B( JR, J+2 ) )
                 B( JR, J ) = B( JR, J ) - TEMP
                 B( JR, J+1 ) = B( JR, J+1 ) - TEMP*V( 2 )
                 B( JR, J+2 ) = B( JR, J+2 ) - TEMP*V( 3 )
  270         CONTINUE
               B( J+1, J ) = ZERO
               B( J+2, J ) = ZERO
			   JMP = JMP + 3
  290 CONTINUE
      RETURN

	END
      SUBROUTINE QZEARLY( N, S, PARA, H, LDH, T, LDT, KBOT, SR, SI,
     $                    SBETA, Q, LDQ, Z, LDZ, DWORK, LDWORK, INFO )
      IMPLICIT NONE
C
C     PURPOSE
C
C     To perform aggressive early deflation on the n-by-n
C     Hessenberg-triangular pencil (H,T). This routine computes
C     a possibly reduced Hessenberg-triangular pencil (HR,TR) and
C     orthogonal matrices Q, Z so that
C
C               (1) H*Z = Q*HR,  T*Z = Q*TR,
C
C               (2) the trailing principal subpencil
C                      ( HR(KBOT:N,KBOT:N), TR(KBOT:N,KBOT:N) )
C                   is in generalized Schur form, normalized as
C                   returned by the LAPACK routine DHGEQZ,
C
C               (3) the components of S*Q(1,KBOT:N,1) are small enough
C                   to be neglected (that also depends on the parameters
C                   in PARA, see TTQZ), and if KBOT is returned with a
C                   value of zero, then S*Q(1,1) is also negligible.
C
C     The procedure is designed for deflating eigenvalues during the QZ
C     iteration.
C
C     Things NOT done:
C     - Implementation of different deflation strategies
C     - deflation of infinite eigenvalues (necessary?)
C     - use of modified QZ for small-sized problem
C     - proper LDWORK checking
C
C     ARGUMENTS
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the matrices H and T.  N >= 0.
C
C     S       (input) DOUBLE PRECISION
C             Scaling factor.
C
C     PARA    (input) DOUBLE PRECISION array, dimension (8)
C             This array contains all parameters of the multishift
C             pipelined QZ algorithm with aggressive deflation.
C             For details, see the routine TTQZ.
C
C     H       (input/output) DOUBLE PRECISION array, dimension (LDH,N)
C             On entry, the leading N-by-N part of this array must
C             contain the upper Hessenberg matrix H.
C             On exit, the leading N-by-N part of this array contains
C             the (possibly reduced) upper Hessenberg matrix HR.
C
C     LDH     (input) INTEGER
C             The leading dimension of the array H. LDH >= max(1,N).
C
C     T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)
C             On entry, the leading N-by-N part of this array must
C             contain the upper triangular matrix T.
C             On exit, the leading N-by-N part of this array contains
C             the upper triangular matrix TR.
C
C     LDT     (input) INTEGER
C             The leading dimension of the array T. LDT >= max(1,N).
C
C     KBOT    (output) INTEGER
C             On exit, KBOT specifies that the pencil
C             (HR(KBOT:N,KBOT:N),TR(KBOT:N,KBOT:N)) is in generalized
C             Schur form. If the entire vector S*Q(1,:) is negligible,
C             then KBOT is returned with a value of zero.
C
C     SR      (output) DOUBLE PRECISION array, dimension (N)
C     SI      (output) DOUBLE PRECISION array, dimension (N)
C     SBETA   (output) DOUBLE PRECISION array, dimension (N)
C             On exit, the leading N elements of the arrays SR, SI and
C             SBETA contain the eigenvalues of the pencil (HR,TR) as
C             returned by the LAPACK routine DHGEQZ. The last N-KBOT+1
C             elements correspond to the eigenvalues of the subpencil
C             (HR(KBOT:N,KBOT:N),TR(KBOT:N,KBOT:N)).
C
C     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
C             On entry, the leading N-by-N part of this array must
C             contain an orthogonal matrix Q.
C             On exit, the leading N-by-N part of this array contains
C             the matrix Q post-multiplied by transpose of the
C             transformations which are applied to H and T on the left.
C
C     LDQ     (input) INTEGER
C             The leading dimension of the array Q.  LDQ >= MAX(1,N).
C
C     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
C             On entry, the leading N-by-N part of this array must
C             contain an orthogonal matrix Z.
C             On exit, the leading N-by-N part of this array contains
C             the matrix Z post-multiplied by the transformations which
C             are applied to H and T on the right.
C
C     LDZ     (input) INTEGER
C             The leading dimension of the array Z.  LDZ >= MAX(1,N).
C
C     Workspace
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) returns the optimal
C             value of LDWORK, ????.
C             On exit, if  INFO = -17, DWORK(1) returns the minimum
C             value of LDWORK.
C
C     LDWORK  INTEGER
C             The length of the array DWORK.  LDWORK >= ???.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal 
C                   value;
C             > 0:  DHGEQZ failed to compute the first INFO eigenvalues
C                   of (H,T). KBOT is set to N and the not converged
C                   eigenvalues are replaced by the diagonal elements
C                   of H and T.
C
C     CONTRIBUTOR
C
C     D. Kressner, Univ. Umea, Sweden, Sept. 2002.
C
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, KBOT, LDH, LDQ, LDT, LDWORK, LDZ, N
      DOUBLE PRECISION  S
C     .. Array Arguments ..
      DOUBLE PRECISION  DWORK(LDWORK), H(LDH,*), PARA(*), Q(LDQ,*),
     $                  SBETA(*), SI(*), SR(*), T(LDT,*), Z(LDZ,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, OVFL, RTDET, SMLNUM, TAU, TB1,
     $                  TB2, TI1, TI2, TR1, TR2, TST1, ULP, UNFL
      INTEGER           I, IERR, IFST, ILST, J, KNT, WRKMIN
      LOGICAL           AGGRRQ, BULGE, DFLATE
C     .. External Functions ..
      DOUBLE PRECISION  DLAMCH
      EXTERNAL          DLAMCH
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGGHRD, DHGEQZ, DLABAD, DLARF,
     $                  DLARFG, DLARTG, DLASET, DROT, DTGEXC, XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, INT, MAX, SQRT
C
C     .. Executable Statements ..
C
C     Decode and check the scalar input parameters.
C
      INFO = 0
      WRKMIN = 4*N + 16
      AGGRRQ = PARA(9).NE.ZERO
C
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSE IF ( LDH.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF ( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF ( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -13
      ELSE IF ( LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -15
      ELSE IF ( LDWORK.LT.WRKMIN ) THEN
         INFO = -17
         DWORK(1) = DBLE( WRKMIN )
      END IF
C
C     Return if there were illegal values.
C
      IF ( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGESQB', -INFO )
         RETURN
      END IF
C
C     Quick return if N = 0.
C
      IF ( N.EQ.0 ) THEN
         DWORK(1) = ONE
         RETURN
      END IF


C
C     Set machine constants for stopping criterion.
C
      UNFL   = DLAMCH('SAFE MINIMUM')
      OVFL   = ONE/UNFL
      CALL DLABAD( UNFL, OVFL )
      ULP    = DLAMCH( 'PRECISION' )
      SMLNUM = UNFL*( N / ULP )
C
C     Convert to (spike-)triangular form.
C
      CALL DHGEQZ( 'Schur form', 'V', 'V', N, 1, N, H, LDH, T, LDT, SR,
     $             SI, SBETA, Q, LDQ, Z, LDZ, DWORK, LDWORK, INFO )
C
C     Note that the documentation of the parameter INFO in DHGEQZ is
C     not consistent with its implementation.
C     Failure of shift calculation can never happen so that IERR is
C     always less than N+1.
C
      IF ( INFO.GT.1 ) THEN
C
C        Set not converged eigenvalues and leave.
C
         PRINT*, 'WARNING: DHGEQZ DID NOT CONVERGE'
         DO 10 I = 1, INFO
            SR(I) = H(I,I)
            SI(I) = ZERO
            SBETA(I) = T(I,I)
   10    CONTINUE
         KBOT = N
         RETURN
      END IF
C
      ILST = 1
      KBOT = N
      KNT = 0
C
C     Don't do aggressive stuff
C
      IF ( S.EQ.ZERO )  THEN
		KBOT = 0      
		RETURN	
	END IF
C      IF (N.LE.25) THEN
C         KBOT = N
C         RETURN
C      END IF
C     Check spike for deflated eigenvalues.
C     while ( knt < n )
C
   20 CONTINUE
      IF ( KNT.LT.N ) THEN
C
         IF ( KBOT.EQ.1 ) THEN
            BULGE = .FALSE.
         ELSE
            BULGE = H(KBOT, KBOT-1).NE.ZERO
         END IF
C
         IF ( .NOT.BULGE ) THEN
 
            KNT = KNT + 1
            TST1 = ABS( S ) + ABS( H(KBOT,KBOT) )
C		  TST1 innehåller H(-1, 0) + H_konv(N,N)
C		  UPL är maskinprecisionen
C		  TST1 motsvarar H^44 (del av Q^T*H33*Z)  I D kressners paper
C		  S*Q(1,QBOT) mot s4
C		  dvs om den förändring som sker på H(-1,0)
C		  pga reduceringen (S*Q(1,KBOT) för att egenvärdet på H_konv(KBOT,KBOT)
C		  ska framträda är mindre än det poteniella avrundningsfelet för
C		  den reducering som verkligen görs så antar man att denna är OK.  
            DFLATE = ABS( S*Q(1,KBOT)).LE.MAX( ULP*TST1, SMLNUM )
         ELSE
            KNT = KNT + 2
            IF ( ( H(KBOT-1,KBOT-1 ).EQ.ZERO )  .AND.
     $           ( H(KBOT,KBOT-1).EQ.ZERO ) ) THEN
               RTDET = ZERO
            ELSE IF ( ABS( H(KBOT-1,KBOT-1) ).GE.
     $                ABS( H(KBOT,KBOT-1) ) ) THEN
               RTDET = SQRT( ABS( H(KBOT-1,KBOT-1) ) )*
     $                 SQRT( ABS( H(KBOT,KBOT) - ( H(KBOT,KBOT-1)/
     $                       H(KBOT-1,KBOT-1) )*H(KBOT-1,KBOT) ) )
            ELSE
               RTDET = SQRT( ABS( H(KBOT,KBOT-1) ) )*
     $                 SQRT( ABS( H(KBOT-1,KBOT) - ( H(KBOT-1,KBOT-1)/
     $                       H(KBOT,KBOT-1) )*H(KBOT,KBOT) ) )
            END IF
            TST1 = ABS(S) + RTDET
            DFLATE = MAX( ABS( S*Q(1,KBOT) ), ABS( S*Q(1,KBOT-1) ) ).LE.
     $               MAX( ULP*TST1, SMLNUM )
         END IF
C
         IF ( DFLATE ) THEN
C
C           Deflation.
C
            IF ( BULGE ) THEN
               KBOT = KBOT - 2
            ELSE
               KBOT = KBOT - 1
            END IF
         ELSE
C
C           Move undeflatable eigenvalue up.
C
            IFST = KBOT
            CALL DTGEXC( .TRUE., .TRUE., N, H, LDH, T, LDT, Q, LDQ, Z,
     $                   LDZ, IFST, ILST, DWORK, LDWORK, IERR )
C
C           Go directly to reduction to Hessenberg-triangular form if
C           eigenvalue swapping failed.
C
            IF ( IERR.NE.0 ) THEN
                PRINT*, 'WARNING: SWAPPING FAILED with reason:',IERR
                GO TO 42
            END IF
C
            IF ( .NOT.BULGE ) THEN
               TR1 = SR(KBOT)
               TI1 = SI(KBOT)
               TB1 = SBETA(KBOT)
               DO 30 I = KBOT - 1, ILST, -1
                  SR(I+1) = SR(I)
                  SI(I+1) = SI(I)
                  SBETA(I+1) = SBETA(I)
   30          CONTINUE
               SR(ILST) = TR1
               SI(ILST) = TI1
               SBETA(ILST) = TB1
               ILST   = ILST + 1
            ELSE
               TR1 = SR(KBOT-1)
               TI1 = SI(KBOT-1)
               TB1 = SBETA(KBOT-1)
               TR2 = SR(KBOT)
               TI2 = SI(KBOT)
               TB2 = SBETA(KBOT)
               DO 40 I = KBOT - 2, ILST, -1
                  SR(I+2) = SR(I)
                  SI(I+2) = SI(I)
                  SBETA(I+2) = SBETA(I)
   40          CONTINUE
               SR(ILST) =  TR1
               SI(ILST) = -TI1
               SBETA(ILST) = TB1
               SR(ILST+1) =  TR2
               SI(ILST+1) =  TI2
               SBETA(ILST+1) = TB2
               ILST = ILST + 2
            END IF
         END IF
C
         GO TO 20
C
      END IF
C
C     Return to Hessenberg form.
C
   42 CONTINUE
C
      IF ( KBOT.GT.0 ) THEN
         CALL DCOPY( KBOT, Q, LDQ, DWORK, 1 )
         ALPHA = DWORK(1)
         CALL DLARFG( KBOT, ALPHA, DWORK(2), 1, TAU )
         DWORK(1) = ONE
         CALL DLARF( 'Left', KBOT, N, DWORK, 1, TAU, H, LDH,
     $               DWORK(N+1) )
         CALL DLARF( 'Right', N, KBOT, DWORK, 1, TAU, Q, LDQ,
     $               DWORK(N+1) )
C
         IF ( KBOT.GT.1 .AND. AGGRRQ ) THEN
            CALL DLARF( 'Left', KBOT, N, DWORK(1), 1, TAU, T, LDT,
     $                  DWORK(N+1) )
            CALL DGERQF( KBOT, KBOT, T, LDT, DWORK(1), DWORK(KBOT+1),
     $                   LDWORK-KBOT, IERR )
            CALL DORMRQ( 'Right', 'Transpose', KBOT, KBOT, KBOT, T, LDT,
     $                   DWORK(1), H, LDH, DWORK(KBOT+1), LDWORK-KBOT,
     $                   IERR )
            CALL DORMRQ( 'Right', 'Transpose', N, KBOT, KBOT, T, LDT,
     $                   DWORK(1), Z, LDZ, DWORK(KBOT+1), LDWORK-KBOT,
     $                   IERR )
            CALL DLASET( 'Lower', KBOT-1, KBOT-1, ZERO, ZERO, T(2,1),
     $                   LDT )
         ELSE IF ( KBOT.GT.1 ) THEN
C
C           The updated matrix (I-tau*u*u') T = T - tau*u*(T'*u)' is a
C           rank-1 perturbation of an upper triangular matrix. Thus,
C           an RQ-update consisting of 2*KBOT-2 Givens rotations can be
C           employed.
C
            IF ( N.GT.KBOT ) THEN
               CALL DLARF( 'Left', KBOT, N-KBOT, DWORK(1), 1, TAU,
     $                     T(1,KBOT+1), LDT, DWORK(N+1) )
            END IF
            CALL DGEMV( 'Transpose', KBOT, KBOT, -TAU, T, LDT, DWORK, 1,
     $                  ZERO, DWORK(KBOT+1), 1 )
C
            DO 17 I = 1, KBOT-1
               ALPHA = DWORK(KBOT+I+1)
               CALL DLARTG( ALPHA, DWORK(KBOT+I), DWORK(2*KBOT+I),
     $                      DWORK(3*KBOT+I), DWORK(KBOT+I+1) )
               DWORK(3*KBOT+I) = -DWORK(3*KBOT+I)
               CALL DROT( I+1, T(1,I), 1, T(1,I+1), 1, DWORK(2*KBOT+I),
     $                    DWORK(3*KBOT+I) )
   17       CONTINUE
            CALL DAXPY( KBOT, DWORK(2*KBOT), DWORK, 1, T(1,KBOT), 1 )
            DO 18 I = 1, KBOT-1
               CALL DROT( KBOT, H(1,I), 1, H(1,I+1), 1, DWORK(2*KBOT+I),
     $                    DWORK(3*KBOT+I) )
   18       CONTINUE
            DO 19 I = 1, KBOT-1
               CALL DROT( N, Z(1,I), 1, Z(1,I+1), 1, DWORK(2*KBOT+I),
     $                    DWORK(3*KBOT+I) )
   19       CONTINUE
C
            DO 28 I = KBOT-1, 1, -1
               ALPHA = T(I+1,I+1)
               CALL DLARTG( ALPHA, T(I+1,I), DWORK(I), DWORK(KBOT+I),
     $                      T(I+1,I+1) )
               DWORK(KBOT+I) = -DWORK(KBOT+I)
               T(I+1,I) = ZERO
               CALL DROT( I, T(1,I), 1, T(1,I+1), 1, DWORK(I),
     $                    DWORK(KBOT+I) )
   28       CONTINUE
            DO 21 I = KBOT-1, 1, -1
               CALL DROT( KBOT, H(1,I), 1, H(1,I+1), 1, DWORK(I),
     $                    DWORK(KBOT+I) )
   21       CONTINUE
            DO 22 I = KBOT-1, 1, -1
               CALL DROT( N, Z(1,I), 1, Z(1,I+1), 1, DWORK(I),
     $                    DWORK(KBOT+I) )
   22       CONTINUE
         END IF
C
C        Hessenberg-triangular reduction. 
C
         CALL DGGHRD( 'V', 'V', N, 1, KBOT, H, LDH, T, LDT, Q,
     $                LDQ, Z, LDZ, IERR ) 
      END IF
C
C     Last line of QZEARLY
C
      END

      SUBROUTINE QZFCOL( N, NSHIFT, PARA, SR, SI, SBETA, H, LDH, T, LDT,
     $                   V, INFO )
      IMPLICIT NONE
C
C     Forms a multiple of the first column of the shift polynomial for
C     the generalized eigenvalue problem.
C     To do: Avoidance of overflow while shifting. Documentation.
C     Check of inputs. Alt. implementation based on matrix-vector-mult.
C     (??)
C     It is assumed that NSHIFT is even and >= 2, <= MIN(12,N-1).
C
C     .. Parameters ..
      INTEGER           BULGMX
      PARAMETER         ( BULGMX = 12 )
      INTEGER           LDTMP
      PARAMETER         ( LDTMP = BULGMX+1 )
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
C     .. Scalar Arguments ..
      INTEGER           INFO, LDH, LDT, N, NSHIFT
C     .. Array Arguments ..
      DOUBLE PRECISION  H(LDH,*), PARA(*),
     $                  SBETA(*), SI(*), SR(*), T(LDT,*), V(*)
C     .. Local Scalars ..
      INTEGER           I, J, NGIV
      DOUBLE PRECISION  ALPHA, CSE, SCAL, SNE
C     .. Local Arrays ..
      DOUBLE PRECISION  TMP(LDTMP,LDTMP), CS1(BULGMX+1),
     $                  CS2(BULGMX+1), SN1(BULGMX+1), SN2(BULGMX+1)
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DLACPY, DLARTG, DROT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, INT, MAX, SQRT
C
C     .. Executable Statements ..
C
      INFO = 0
      NGIV = 0
      DO 100  I = NSHIFT, 2, -2
C
         IF ( NGIV.GT.0 ) THEN
            CALL DLACPY( 'All', NGIV+1, NGIV+1, T, LDT, TMP, LDTMP )
            DO 10  J = NGIV, 1, -1
               CALL DROT( NGIV-J+2, TMP(J,J), LDTMP, TMP(J+1,J), LDTMP,
     $                    CS1(J), SN1(J) )
               ALPHA = TMP(J+1,J+1)
               CALL DLARTG( ALPHA, TMP(J+1,J), CS2(J), SN2(J),
     $                      TMP(J+1,J+1) )
               SN2(J) = -SN2(J)
               CALL DROT( J, TMP(1,J), 1, TMP(1,J+1), 1, CS2(J),
     $                    SN2(J) )
   10       CONTINUE
            SCAL = TMP(1,1)
         ELSE 
            SCAL = T(1,1)
         END IF
C
         CALL DLACPY( 'All', NGIV+2, NGIV+2, H, LDH, TMP, LDTMP )
         DO 20  J = 1, NGIV+2
            CALL DSCAL( MIN( J+1, NGIV+2 ), SBETA(I), TMP(1,J), 1 )
            CALL DAXPY( J, -SR(I), T(1,J), 1, TMP(1,J), 1 )
   20    CONTINUE
         DO 30  J = NGIV, 1, -1
            CALL DROT( NGIV+2, TMP(1,J), 1, TMP(1,J+1), 1, CS2(J),
     $                 SN2(J) )
   30    CONTINUE
         DO 40  J = NGIV+1, 1, -1
            ALPHA = TMP(J,1)
            CALL DLARTG( ALPHA, TMP(J+1,1), CS2(J), SN2(J), TMP(J,1) )
   40    CONTINUE
         NGIV = NGIV + 1
C
         CALL DLARTG( TMP(1,1), SI(I)*SCAL, CSE, SNE, ALPHA )
         CALL DLACPY( 'All', NGIV+1, NGIV+1, T, LDT, TMP, LDTMP )
         DO 50  J = NGIV, 1, -1
            CALL DROT( NGIV-J+2, TMP(J,J), LDTMP, TMP(J+1,J), LDTMP,
     $                 CS2(J), SN2(J) )
            ALPHA = TMP(J+1,J+1)
            CALL DLARTG( ALPHA, TMP(J+1,J), CS2(J), SN2(J),
     $                   TMP(J+1,J+1) )
            SN2(J) = -SN2(J)
            CALL DROT( J, TMP(1,J), 1, TMP(1,J+1), 1, CS2(J), SN2(J) )
   50    CONTINUE
         CALL DLARTG( CSE, -SNE*TMP(1,1), CSE, SNE, ALPHA )
         SNE = -SNE
C
         CALL DLACPY( 'All', NGIV+2, NGIV+2, H, LDH, TMP, LDTMP )
         DO 60  J = 1, NGIV+2
            CALL DSCAL( MIN( J+1, NGIV+2 ), SBETA(I-1), TMP(1,J), 1 )
            CALL DAXPY( J, -SR(I-1), T(1,J), 1, TMP(1,J), 1 )
   60    CONTINUE
         DO 70  J = NGIV, 1, -1
            CALL DROT( NGIV+2, TMP(1,J), 1, TMP(1,J+1), 1, CS2(J),
     $                 SN2(J) )
   70    CONTINUE
         CALL DLASET( 'All', NGIV+2, NGIV+1, ZERO, -SI(I-1), TMP(1,2),
     $                LDTMP )
         DO 80  J = NGIV-1, 1, -1
            CALL DROT( NGIV, TMP(1,J+1), 1, TMP(1,J+2), 1, CS1(J),
     $                 SN1(J) )
   80    CONTINUE
         CALL DROT( NGIV+2, TMP(1,1), 1, TMP(1,2), 1, CSE, SNE )
         IF ( I.GT.2 ) THEN
            DO 90  J = NGIV+1, 1, -1
               ALPHA = TMP(J,1)
               CALL DLARTG( ALPHA, TMP(J+1,1), CS1(J), SN1(J),
     $                      TMP(J,1) )
   90       CONTINUE
            NGIV = NGIV + 1
         ELSE
            CALL DCOPY( NGIV+2, TMP(1,1), 1, V, 1 )
         END IF
  100 CONTINUE

C
C     Last line of QZFCOL
C
      END
      SUBROUTINE STEP1(COMPQ, COMPZ, NOBLK, N, ILO, IHI , INFO, 
     $     A, DESCA, B, DESCB,
     $     Q, DESCQ, Z, DESCZ, TAU, WORK, NB, TT, LWORK)
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            NOBLK, N, ILO, IHI, NB, LWORK, INFO
      LOGICAL            COMPQ, COMPZ

*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   TAU(*), TT(*)
      DOUBLE PRECISION   WORK(*), A(*), B(*), Z(*), Q(*)
      INTEGER            DESCA(*), DESCB(*)
      INTEGER		 DESCZ(*), DESCQ(*)
*     ..
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          (ONE = 1.D0, ZERO = 0.D0)
*     ..
*     .. Local Arrays ..

*     ..
*     .. Local Scalars ..

      INTEGER            I, J, LASTBLK, NNB
      INTEGER            ICTXT, NPROW, NPCOL, MYROW, MYCOL, OLDI, NOCOL
      INTEGER            K, II, NOIT, NOROW
      INTEGER            IX, JX, IS, IA, IIA, JJA
      INTEGER            IAROW, IACOL , IH, JH

*     ..
*     .. External Subroutines ..
      EXTERNAL           PDGEQR2, PDLARFT, PDLARFB, INFOG2L
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, PXERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, MOD
*     ..
*     .. External Functions

*     ..
*     .. Executable Statements ..
*     ..


*     ..
*     .. Get grid parameters
*     ..
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL)
*     ..
*     .. Test the input parameters
*     ..
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -(700+CTXT_)
*      ELSE
*         CALL CHK1MAT( N, 1, N, 2, 1, 1, DESCA, 4, INFO )
*         CALL CHK1MAT( N, 1, N, 2, 1, 1, DESCB, 6, INFO )
      END IF
      
      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'STEP1', -INFO )
         RETURN
      END IF

*     ..
*     .. Quick Return
*     ..
      IF (N .LE. NB +1) THEN
         RETURN
      ENDIF
      
*     ..
*     .. Annihilate Blocks of A
*     ..
      DO 10 J = ILO, IHI-NB, NB
         LASTBLK = MOD((IHI-ILO)+1, NB)
         
         IF (LASTBLK .EQ. 0) THEN
            LASTBLK = NB
         ENDIF
         
         I = MAX(N- ((NOBLK-1)*NB+LASTBLK)+1, J + NB)
         NOIT = (I- (J +NB)) / ((NOBLK-1)*NB) + 1
         
         IF (MOD((I- (J +NB)),  ((NOBLK-1)*NB)) .NE. 0) THEN
            NOIT = NOIT + 1
         ENDIF
         
         NOROW = MIN(NOBLK*NB, N - I + 1)
         
         DO 20 IA = 1, NOIT
            NNB = MIN(NB, NOROW)
            CALL PDGEQR2(NOROW, NNB, A, I, J, DESCA, TAU,
     $           WORK, LWORK, INFO )
            CALL PDLARFT( 'Forward', 'Columnwise', NOROW, NNB, A, I, J,
     $           DESCA, TAU, TT, WORK )
            CALL PDLARFB( 'Left', 'Transpose', 'Forward', 'Columnwise',
     $           NOROW, N-(J+NNB)+1, NNB, A, I, J, DESCA, TT,
     $           A, I, J+NNB, DESCA, WORK )
            CALL PDLARFB( 'Left', 'Transpose', 'Forward', 'Columnwise',
     $           NOROW, N-I+1, NNB, A, I, J, DESCA, TT, B,
     $           I, I, DESCB, WORK )

*	    ..
*           .. Update Q if desired
*	    ..
            IF (COMPQ) THEN
               CALL PDLARFB( 'Right', 'No Transpose','Forward', 
     $              'Columnwise',
     $              N, NOROW , NNB, A, I, J, DESCA, TT,
     $              Q, 1, I, DESCQ, WORK )
            ENDIF

*           ..
*           .. Set triangular part to zero
*           ..
            CALL INFOG2L( I, J, DESCA, NPROW, NPCOL,
     $           MYROW, MYCOL, IIA, JJA,IAROW, IACOL )
            IF( MYROW.EQ.IAROW .AND. MYCOL.EQ.IACOL ) THEN
               K = 0
               DO 30 JX = JJA, JJA+ NNB-1
                  K = K + 1
                  DO 31 IX = IIA+K, IIA + NNB - 1
                     A(IX + (JX-1)*DESCA( LLD_ ) ) = ZERO
 31               CONTINUE
 30            CONTINUE
            ENDIF
*           ..
*           .. Set rectangular part to zero
*           ..
            DO 60 IH = I + NB, I + NOROW-1, NB
               CALL INFOG2L( IH, J, DESCA, NPROW, NPCOL,
     $              MYROW, MYCOL, IIA, JJA,IAROW, IACOL )
               IF( MYROW.EQ.IAROW .AND. MYCOL.EQ.IACOL ) THEN
                  DO 70 JX = JJA, JJA + NB - 1
                     DO 71 IX = IIA, IIA + MIN(NB,  N - IH +1) - 1
                        A(IX + (JX-1)*DESCA( LLD_ ) ) = ZERO
 71                  CONTINUE
 70               CONTINUE
               ENDIF
 60         CONTINUE
            
*           ..
*           .. Delete B Fill in 
*           ..
            IF (I .EQ. J + NB) THEN
               IS = I
            ELSE
               IS = I + NB
            ENDIF
            NOCOL = NOROW
            DO 40 II = I + NOROW - LASTBLK, IS, -NB
               NNB = MIN(NB, IHI-II+1)
               CALL PDGERQ2(NNB,NOCOL,B, II, I, DESCB,
     $                TAU,WORK,LWORK,INFO)
               CALL PDLARFT('Backward','Rowwise', NOCOL, NNB, B, II, I,
     $              DESCB, TAU, TT, WORK )
               CALL PDLARFB( 'Right', 'No transpose', 'Backward',
     $             'Rowwise', N, NOCOL, NNB, B , II, I, DESCB, TT,
     $             A, 1, I, DESCA, WORK )
               CALL PDLARFB( 'Right', 'No transpose', 'Backward',
     $              'Rowwise', II-1, NOCOL, NNB, B, II, I, DESCB, TT,
     $              B, 1, I, DESCB, WORK )

*	       ..
*	       .. Update Z if desired
*	       ..
               IF (COMPZ) THEN
                  CALL PDLARFB( 'Right', 'No transpose', 'Backward',
     $                 'Rowwise', N, NOCOL, NNB, B, II, I, DESCB, TT,
     $                 Z, 1, I, DESCZ, WORK )
               ENDIF
			
*              ..
*              .. Set triangular part to zero
*              ..
               CALL INFOG2L( II, II, DESCA, NPROW, NPCOL,
     $                   MYROW, MYCOL, IIA, JJA,IAROW, IACOL )
               IF( MYROW.EQ.IAROW .AND. MYCOL.EQ.IACOL ) THEN
                   K = 0 
                   DO 50 JX = JJA, JJA + NNB - 1
                      K = K + 1
                      DO 51 IX = IIA + K, IIA + NNB - 1
                         B(IX + (JX-1)*DESCA( LLD_ ) ) = ZERO
 51                   CONTINUE
 50                CONTINUE
               ENDIF
*              ..
*              .. Set rectangular part to zero
*              ..
               DO 80 JH = I, II-LASTBLK, NB
                  CALL INFOG2L( II, JH, DESCA, NPROW, NPCOL,
     $                 MYROW, MYCOL, IIA, JJA,IAROW, IACOL )
                  IF( MYROW.EQ.IAROW .AND. MYCOL.EQ.IACOL ) THEN
                     DO 90 JX = JJA, JJA + NB - 1
                        DO 91 IX = IIA, IIA + NNB - 1
                           B(IX + (JX-1)*DESCA( LLD_ ) ) = ZERO
 91                     CONTINUE
 90                  CONTINUE
                  ENDIF
 80            CONTINUE
               NOCOL = NOCOL - LASTBLK
               LASTBLK = NB
 40         CONTINUE
            OLDI = I
            I = MAX(I - (NOBLK-1)*NB, J+ NB)
            NOROW = OLDI - I + NB
 20      CONTINUE
 10   CONTINUE
 99   RETURN
      END
      SUBROUTINE STEP2(COMPQ, COMPZ, N, ILO, IHI, A, DESCA, B, DESCB,
     $     NB, INFO, CSROW, CSCOL, TMPARR, MAXNOJ,
     $     Q, Z, DESCQ, DESCZ, LWORK)
      
      IMPLICIT NONE
*     
*     .. Scalar Arguments ..
      INTEGER            N, LWORK
      INTEGER            NB, INFO, ILO, IHI, MAXNOJ
      LOGICAL            COMPQ, COMPZ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A(*), B(*), TMPARR(*)      
      DOUBLE PRECISION   Q(*), Z(*)
      DOUBLE PRECISION   CSROW(*), CSCOL(*)
      INTEGER            DESCA(*), DESCB(*)
      INTEGER            DESCQ(*), DESCZ(*)

*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      

      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $     LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      

*     ..
*     .. Local Arrays
      DOUBLE PRECISION   AA(20,20), BB(20,20), WORK(20)

*     ..
*     .. Local Scalars ..
      INTEGER             I, J, K, L, JJ, IIJ
      INTEGER             NNB, NB1, LDA, LDB, LDC, LDTMP
      INTEGER		  LDQ, LDZ
      INTEGER             NOJ, IC, IX, IND, JC, IY, NN
      INTEGER             ICTXT, NPROW, NPCOL, MYROW, MYCOL
      INTEGER             IAM

*     ..
*     .. External Subroutines ..
      EXTERNAL            BLACS_GRIDINFO
      EXTERNAL            BLACS_PINFO
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC           MAX, MIN, MOD        

*     ..
*     .. External Functions
      EXTERNAL            INDXG2L, INDXG2P
      INTEGER             INDXG2L, INDXG2P

*     ..
*     .. Executable Statements ..
*     ..

*     ..
*     .. Test the input parameters
*     ..

      INFO = 0
      IF(NPROW .EQ. -1 .OR. NPCOL .EQ. -1) THEN
         INFO = -(700+CTXT_)
      END IF

*     ..
*     .. Get grid parameters
*     ..
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
      IAM = MYROW + MYCOL * NPCOL

      LDA = DESCA(LLD_)
      LDB = DESCB(LLD_)
      LDQ = DESCQ(LLD_)
      LDZ = DESCQ(LLD_)
      LDC = MAXNOJ * 2
      LDTMP = NB

      NOJ = 0
*      IHI = 10 
      L = MIN(ILO + NB - 1, IHI - 1)
      IC = ILO + 1
      JC = ILO
      NN = IHI - ILO + 1


      DO 10 J = ILO, IHI - 2
         NOJ = NOJ + 1
         IF (NOJ .GT. 1) THEN
*           .. Move all bulges one step ahead
*           .. Apply saved cosine and sinus values to be able 
*           .. to reduce futher            
            JC = J - NOJ + 1
            I = IC
            DO 20 IX = 1, NOJ - 1
               IND = (IX - 1) * 2 + 1
               NNB = MIN(NB, IHI - NB - I + 1)
               NB1 = MIN(NB, IHI - I + 1)
               
               IF (NB1 .GT. 0) THEN
*                 ..
*                 .. Apply row rotations to Diagonal block of B
*                 ..                  
                  IF ((NPROW * NPCOL) .EQ. 1) THEN
                     CALL MVBB(NB1, B(I + LDB * (I - 1)), LDB, 
     $                    CSROW(IND + LDC * (I - 1)),
     $                    CSCOL(IND + LDC * (I - 1)), LDC)
                  ELSE                   
                     CALL pMVBB(B, LDB, CSROW, CSCOL, LDC, I, JC, NB,
     $                    NB1, NN, NPROW, NPCOL, IND, TMPARR, LDTMP,
     $                    MYROW,MYCOL,ICTXT)
                  ENDIF

*                 ..
*                 .. Apply Rotations to block column  of B (above diagonal block)
*                 ..
                  IF ((NPROW * NPCOL) .EQ. 1) THEN
                     CALL UPD(I - 1, NB1, B(1 + LDB * (I - 1)), LDB, NB,
     $                    JC, CSROW(IND), CSCOL(IND + LDC * (I - 1)),
     $                    LDC)
                  ELSE
                     CALL pUPDB(B, LDB, CSROW, CSCOL, LDC, I, JC, I - 1,
     $                    NB, NB1, NPROW, NPCOL, IND, TMPARR, LDTMP,
     $                    MYROW, MYCOL, ICTXT)
                  ENDIF


*                 ..
*                 .. Apply Rotations to block column  of A (above diagonal block)
*                 ..
                  IF ((NPROW * NPCOL) .EQ. 1) THEN
                     CALL UPD(MIN(I + NB - 1, IHI), NB1, 
     $                    A(1 + LDA * (I - 1)), LDA, NB, JC, CSROW(IND),
     $                    CSCOL(IND + LDC * (I - 1)), LDC)
                  ELSE
                     CALL pUPDA(A, LDA, CSROW, CSCOL, LDC,
     $                    I, JC, MIN(I + NB - 1, IHI), NB, NB1, NPROW, 
     $                    NPCOL, IND, TMPARR, LDTMP, 
     $                    MYROW, MYCOL, ICTXT)
                  ENDIF
               ENDIF
*              ..
*              .. Apply col-roations to diagonal block of A 
*              .. 
               IF (NNB .GT. 0) THEN
                  IF ((NPROW * NPCOL) .EQ. 1) THEN
                     CALL MVBA(NNB, NB, A(I + NB + LDA * (I - 1)), LDA, 
     $                    CSROW(IND + LDC * (I + NB - 1)),
     $                    CSCOL(IND + LDC * (I - 1)), LDC)
                  ELSE
                     CALL pMVBA(A, LDA, CSROW, CSCOL, LDC, I, NB,
     $                    NNB, NN, NPROW, NPCOL, IND, TMPARR, LDTMP, 
     $                    COMPQ,MYROW, MYCOL, ICTXT)
                  ENDIF
               ENDIF
               I = I - (NB - 1)
               JC = JC + 1
 20         CONTINUE
            IC = IC + NB
         ENDIF

*        ..
*        .. Annihilate Column J
*        ..
         IND = (NOJ - 1) * 2 + 1

         IF ((NPROW * NPCOL) .EQ. 1) THEN
            CALL KILLCOL(L - J + 1, A(J + 1 + LDA *(J - 1)), LDA, 
     $           CSROW(IND + LDC * J), LDC)
         ELSE
            CALL pKILLCOL(A, LDA, CSROW, CSCOL, LDC,
     $           J, L, NB, NN, NPROW, NPCOL, IND,
     $           TMPARR, LDTMP, COMPQ,MYROW,MYCOL,ICTXT)
         ENDIF


         IF ((NOJ .EQ. MAXNOJ) .OR. (J .EQ. (IHI - 2))) THEN

*           ..
*           .. Chase all MAXNOJ bulges down the diagonals of A and B
*           ..
            DO 40 IY = IC - (NOJ - 1)*(NB - 1), IHI, NB
               JC = J - NOJ + 1
               I = IC
               DO 50 IX = 1, NOJ
                  IND = (IX - 1) * 2 + 1
                  NNB = MIN(NB, IHI - NB - I + 1)
                  NB1 = MIN(NB, IHI - I + 1)

                  IF (NB1 .GT. 0) THEN
*                 ..
*                 .. Apply row rotations to Diagonal block of B
*                 .. and delete fill-in
*                 ..
                     IF ((NPROW * NPCOL) .EQ. 1) THEN
                        CALL MVBB(NB1, B(I + LDB * (I-1)), LDB, 
     $                       CSROW(IND + LDC * (I-1)),
     $                       CSCOL(IND + LDC * (I-1)), LDC)
                     ELSE
                        CALL pMVBB(B, LDB, CSROW, CSCOL, LDC, I, JC,
     $                       NB, NB1, NN, NPROW, NPCOL, IND,
     $                       TMPARR, LDTMP, MYROW, MYCOL, ICTXT)
                     ENDIF
*                    ..
*                    .. Apply Rotations to block column  of B (above diagonal block)
*                    ..
                     IF ((NPROW * NPCOL) .EQ. 1) THEN
                        CALL UPD(I - 1, NB1, B(1+LDB*(I-1)), LDB, NB,
     $                       JC, CSROW(IND), CSCOL(IND + LDC * (I - 1)), 
     $                       LDC)
                     ELSE
                        CALL pUPDB(B, LDB, CSROW, CSCOL, LDC,
     $                       I, JC, I - 1, NB, NB1, NPROW, NPCOL, IND, 
     $                       TMPARR, LDTMP, MYROW, MYCOL, ICTXT)
                     ENDIF

*                    ..
*                    .. Apply Rotations to block column  of A (above diagonal block)
*                    ..
                     IF ((NPROW * NPCOL) .EQ. 1) THEN
                        CALL UPD(MIN(I + NB - 1, IHI), NB1, 
     $                       A(1 + LDA * (I - 1)), LDA, NB, JC, 
     $                       CSROW(IND), CSCOL(IND + LDC * (I - 1)),
     $                       LDC)
                     ELSE
                        CALL pUPDA(A, LDA, CSROW, CSCOL, LDC,
     $                       I, JC, MIN(I + NB - 1, IHI), NB, NB1, 
     $                       NPROW, NPCOL, IND, TMPARR, LDTMP, 
     $                       MYROW,MYCOL,ICTXT)
                     ENDIF
                  ENDIF

*                 ..
*                 .. Apply col-roations to diagonal block of A
*                 ..
                  IF (NNB .GT. 0) THEN
                     IF ((NPROW * NPCOL) .EQ. 1) THEN
                        CALL MVBA(NNB, NB, A(I + NB + LDA * (I - 1)),
     $                       LDA, CSROW(IND + LDC * (I + NB - 1)),
     $                       CSCOL(IND + LDC * (I - 1)), LDC)
                     ELSE
                        CALL pMVBA(A, LDA, CSROW, CSCOL, LDC, I,
     $                       NB, NNB, NN, NPROW, NPCOL, IND, TMPARR, 
     $                       LDTMP, COMPQ,MYROW, MYCOL, ICTXT)
                     ENDIF
                  ENDIF
                  I = I - (NB - 1)
                  JC = JC + 1
 50            CONTINUE
               IC = IC + NB
 40         CONTINUE
            IC = J + 2
            JC = J + 1
*           ..
*           .. Compute Q if asked to do so
*           .. 
            IF (COMPQ) THEN
               IF ((NPROW * NPCOL) .EQ. 1) THEN
                  CALL UPDQ(IHI, ILO, N, J, NB, NOJ, Q, LDQ, CSROW, LDC)
               ELSE
                  CALL pUPDQ(IHI, ILO, N, J, NB, NOJ, Q, LDQ, CSROW, 
     $                 LDC, NPROW, NPCOL, MYROW, MYCOL, TMPARR, LDTMP, 
     $                 ICTXT)
               ENDIF
            ENDIF
*           ..
*           .. Compute Z if asked to do so
*           ..
            IF (COMPZ) THEN
               IF ((NPROW * NPCOL) .EQ. 1) THEN
                  CALL UPDZ(IHI, ILO, N, J, NB, NOJ, Z, LDZ, CSCOL, LDC)
               ELSE
                  CALL pUPDZ(IHI, ILO, N, J, NB, NOJ, Z, LDZ, CSCOL,
     $                 LDC, NPROW, NPCOL, MYROW, MYCOL, TMPARR, LDTMP, 
     $                 ICTXT)
               ENDIF
            ENDIF
*           ..
*           .. If the matrix is wider than high update remaining cols
*           ..
            IF (IHI .LT. N) THEN
               JJ = J - NOJ + 1
               DO 60 K = 1, NOJ
                  IND = (K - 1) * 2 + 1
                  IF ((NPROW * NPCOL) .EQ. 1) THEN
                     CALL UPDH(N, IHI, ILO, A,
     $                    LDA, NB, JJ, CSROW(IND), LDC)
                     CALL UPDH(N, IHI, ILO, B,
     $                    LDB, NB, JJ, CSROW(IND), LDC)
                  ELSE
                     CALL pUPDH(N, IHI, ILO, JJ, NB, NOJ, A, LDA, CSROW,
     $                    LDC, NPROW, NPCOL, MYROW, MYCOL, 
     $                    TMPARR, LDTMP, ICTXT)
                     CALL pUPDH(N, IHI, ILO, JJ, NB, NOJ, B, LDB, CSROW,
     $                    LDC, NPROW, NPCOL, MYROW, MYCOL, 
     $                    TMPARR, LDTMP, ICTXT)
                  ENDIF
                  JJ = JJ + 1
 60            CONTINUE
            ENDIF
            NOJ = 0
         ENDIF
         L = MIN(L + 1, IHI - 1)
 10   CONTINUE
*
*     End of STEP2
*
800   FORMAT(24F7.2)
999   RETURN
      END
*
   
      SUBROUTINE UPDH(N, IHI, ILO, A, LDA, NB, J, CSROW, LDC)
      IMPLICIT NONE
      INTEGER IHI, ILO, LDA, NB, J, LAST, LDC
      DOUBLE PRECISION A(LDA, *), CSROW(LDC,*)
      INTEGER K, KK, N
*
*     Apply row rotations
      DO 20 KK = MIN(NB+J, IHI), IHI, NB
         DO 30 K = KK, MAX(KK-NB+2, J+2), -1
            CALL DROT( N-IHI, A( K - 1,IHI+1 ), LDA,
     $           A( K, IHI+1 ), LDA, CSROW(1,K-1), CSROW(2,K-1) )
 30      CONTINUE
         LAST = KK
 20   CONTINUE
      DO 40 K = IHI, LAST+2, -1
C     CALL DROT( N, A( K - 1,IHI+1 ), LDA,
C     $        A( K, IHI+1 ), LDA, CSROW(1,K-1), CSROW(2,K-1) )
C     Ban 010213 Can't rowrot more than n-ihi columns !?! 
         CALL DROT( N-IHI, A( K - 1,IHI+1 ), LDA,
     $        A( K, IHI+1 ), LDA, CSROW(1,K-1), CSROW(2,K-1) )
         
40    CONTINUE
      RETURN
      END

      SUBROUTINE pUPDH(N, IHI, ILO, J, NB, NOJ, A, LDA, CSROW, LDC,
     $     NPROW, NPCOL, MYROW, MYCOL, TMPARR, LDTMP, CTX)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER		IHI, ILO, N, J, NB, NOJ, LDA, LDC
      INTEGER		NPROW, NPCOL, MYROW, MYCOL, LDTMP, CTX

*     .. Array Arguments ..
      DOUBLE PRECISION	A(LDA, *), CSROW(LDC,*), TMPARR(LDTMP)

*     .. Parameters ..
      INTEGER	        DOWN, UP, RIGHT, LEFT
      PARAMETER		( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )

*     .. Local scalars ..
      INTEGER	        K, KK, LAST, I
      INTEGER		bl_owner_r_u, bl_owner_r_d
      INTEGER		bl_owner_c_l, bl_owner_c_r
      INTEGER		bl_index_r_u, bl_index_r_d
      INTEGER		bl_index_c_l, bl_index_c_r
      INTEGER		LCNT, RCNT

*     .. External Functions ..
      EXTERNAL		INDXG2L, INDXG2P
      INTEGER		INDXG2L, INDXG2P

*     ..
*     .. Start execution
*     ..
      

*     ..
*     .. Apply row rotations
*     ..
      DO 20 KK = MIN(NB + J, IHI), IHI, NB
         DO 30 K = KK, MAX(KK - NB + 2, J + 2), -1
            bl_owner_r_u = INDXG2P(K - 1, NB, 0, 0, NPROW)
            bl_owner_r_d = INDXG2P(K, NB, 0, 0, NPROW)
            bl_index_r_u = INDXG2L(K - 1, NB, 0, 0, NPROW)
            bl_index_r_d = INDXG2L(K, NB, 0, 0, NPROW)
            DO 40 I = IHI + 1, N, NB
               LCNT = NB - MOD(I, NB)
               RCNT = MOD(I, NB)
               IF ((I + LCNT) .GT. N) THEN
                  RCNT = 0
                  LCNT = N - I
               ELSE IF ((I + LCNT + RCNT) .GT. N) THEN
                  RCNT = N - I - LCNT
               ENDIF
               
               bl_index_c_l = INDXG2L(I, NB, 0, 0, NPCOL)
               bl_owner_c_l = INDXG2P(I, NB, 0, 0, NPCOL)
               bl_index_c_r = INDXG2L(I + LCNT, NB, 0, 0, NPCOL)
               bl_owner_c_r = INDXG2P(I + LCNT, NB, 0, 0, NPCOL)
               
               
               IF (bl_owner_r_u .EQ. bl_owner_r_d) THEN
                  IF ((MYROW .EQ. bl_owner_r_u) .AND. 
     $                 (MYCOL .EQ. bl_owner_c_l)) THEN
*                    ..
*                    .. rowrot LCNT elements of current rows
*                    ..
                     CALL DROT(LCNT, A(bl_index_r_u, bl_index_c_l), LDA,
     $                    A(bl_index_r_d, bl_index_c_l), LDA, 
     $                    CSROW(1, K - 1), CSROW(2, K - 1))
                  ENDIF
                  IF ((RCNT .GT. 0) .AND. (MYROW .EQ. bl_owner_r_u) 
     $                 .AND. (MYCOL .EQ. bl_owner_c_r)) THEN
*                    ..
*                    .. rowrot LCNT elements of current rows
*                    ..
                     CALL DROT(RCNT, A(bl_index_r_u, bl_index_c_r), LDA,				
     $                    A(bl_index_r_d, bl_index_c_r), LDA, 
     $                    CSROW(1, K - 1), CSROW(2, K - 1))
                  ENDIF
               ELSE
                  IF ((MYROW .EQ. bl_owner_r_d) .AND. 
     $                 (MYCOL .EQ. bl_owner_c_l)) THEN
*                    ..
*                    .. Exchange part of row with proc located above before preforming rotation
*                    ..
                     CALL SEND(1, LCNT, UP, 
     $                    A(bl_index_r_d, bl_index_c_l), LDA,	
     $                    MYROW, MYCOL, NPCOL, NPROW, CTX)
                     CALL RECVTMP(1, LCNT, UP, TMPARR, LDTMP, 
     $                    MYROW, MYCOL, NPCOL, NPROW, CTX)
                     CALL DROT(LCNT, TMPARR, LDTMP, 
     $                    A(bl_index_r_d, bl_index_c_l), LDA,
     $                    CSROW(1, K - 1), CSROW(2, K - 1))
                  ENDIF
                  IF ((RCNT .GT. 0) .AND. (MYROW .EQ. bl_owner_r_d)
     $                 .AND. (MYCOL .EQ. bl_owner_c_r)) THEN
*                    ..
*                    .. Exchange part of row with proc located above before preforming rotation
*                    ..
                     CALL SEND(1, RCNT, UP, 
     $                    A(bl_index_r_d, bl_index_c_r), LDA,	
     $                    MYROW, MYCOL, NPCOL, NPROW, CTX)
                     CALL RECVTMP(1, RCNT, UP, TMPARR, LDTMP, 
     $                    MYROW, MYCOL, NPCOL, NPROW, CTX)
                     CALL DROT(RCNT, TMPARR, LDTMP, 
     $                    A(bl_index_r_d, bl_index_c_l), LDA, 
     $                    CSROW(1, K - 1), CSROW(2, K - 1))
                  ENDIF
                  IF ((MYROW .EQ. bl_owner_r_u) .AND. 
     $                 (MYCOL .EQ. bl_owner_c_l)) THEN
*                    ..
*                    .. Exchange part of row with proc located below before preforming rotation
*                    ..
                     CALL RECVTMP(1, LCNT, DOWN, TMPARR, LDTMP, 
     $                    MYROW, MYCOL, NPCOL, NPROW, CTX)
                     CALL SEND(1, LCNT, DOWN, 
     $                    A(bl_index_r_u, bl_index_c_l), LDA,	
     $                    MYROW, MYCOL, NPCOL, NPROW, CTX)
                     CALL DROT(LCNT, A(bl_index_r_u, bl_index_c_l), LDA,	
     $                    TMPARR, LDTMP, 
     $                    CSROW(1, K - 1), CSROW(2, K - 1))
                  ENDIF
                  IF ((RCNT .GT. 0) .AND. (MYROW .EQ. bl_owner_r_u) 
     $                 .AND. (MYCOL .EQ. bl_owner_c_r)) THEN
*                    ..
*                    .. Exchange part of row with proc located below before preforming rotation
*                    ..
                     CALL RECVTMP(1, RCNT, DOWN, TMPARR, LDTMP, 
     $                    MYROW, MYCOL, NPCOL, NPROW, CTX)
                     CALL SEND(1, RCNT, DOWN, 
     $                    A(bl_index_r_u, bl_index_c_r), LDA,	
     $                    MYROW, MYCOL, NPCOL, NPROW, CTX)						
                     CALL DROT(RCNT, A(bl_index_r_u, bl_index_c_r), LDA,	
     $                    TMPARR, LDTMP, 
     $                    CSROW(1, K - 1), CSROW(2, K - 1))
                  ENDIF
               ENDIF
 40         CONTINUE
            
 30      CONTINUE
         LAST = KK
 20   CONTINUE
      DO 50 K = IHI, LAST+2, -1
         bl_owner_r_u = INDXG2P(K - 1, NB, 0, 0, NPROW)
         bl_owner_r_d = INDXG2P(K, NB, 0, 0, NPCOL)
         bl_index_r_u = INDXG2L(K - 1, NB, 0, 0, NPROW)
         bl_index_r_d = INDXG2L(K, NB, 0, 0, NPROW)
         DO 60 I = IHI + 1, N, NB
            LCNT = NB - MOD(I, NB)
            RCNT = MOD(I, NB)
            IF ((I + LCNT) .GT. N) THEN
               RCNT = 0
               LCNT = N - I
            ELSE IF ((I + LCNT + RCNT) .GT. N) THEN
               RCNT = N - I - LCNT
            ENDIF
            
            bl_index_c_l = INDXG2L(I, NB, 0, 0, NPCOL)
            bl_owner_c_l = INDXG2P(I, NB, 0, 0, NPCOL)
            bl_index_c_r = INDXG2L(I + LCNT, NB, 0, 0, NPCOL)
            bl_owner_c_r = INDXG2P(I + LCNT, NB, 0, 0, NPCOL)
            
            IF (bl_owner_r_u .EQ. bl_owner_r_d) THEN
               IF ((MYROW .EQ. bl_owner_r_u) .AND. 
     $              (MYCOL .EQ. bl_owner_c_l)) THEN
*                 ..
*                 .. rowrot LCNT elements of current rows
*                 ..
                  CALL DROT(LCNT, A(bl_index_r_u, bl_index_c_l), LDA,				
     $                 A(bl_index_r_d, bl_index_c_l), LDA, 
     $                 CSROW(1, K - 1), CSROW(2, K - 1))
               ENDIF
               
               IF ((RCNT .GT. 0) .AND. (MYROW .EQ. bl_owner_r_u) 
     $              .AND. (MYCOL .EQ. bl_owner_c_r)) THEN
*                 ..
*                 .. rowrot RCNT elements of current rows
*                 ..
                  CALL DROT(RCNT, A(bl_index_r_u, bl_index_c_r), LDA,				
     $                 A(bl_index_r_d, bl_index_c_r), LDA, 
     $                 CSROW(1, K - 1), CSROW(2, K - 1))
               ENDIF
            ELSE
               IF ((MYROW .EQ. bl_owner_r_d) .AND. 
     $              (MYCOL .EQ. bl_owner_c_l)) THEN
*                 ..
*                 .. Exchange part of row with proc located above before preforming rotation
*                 ..
                  CALL SEND(1, LCNT, UP, 
     $                 A(bl_index_r_d, bl_index_c_l), LDA,	
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)
                  CALL RECVTMP(1, LCNT, UP, TMPARR, LDTMP, 
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)
                  CALL DROT(LCNT, TMPARR, LDTMP, 
     $                 A(bl_index_r_d, bl_index_c_l), LDA, 
     $                 CSROW(1, K - 1), CSROW(2, K - 1))
               ENDIF
               IF ((RCNT .GT. 0) .AND. (MYROW .EQ. bl_owner_r_d) 
     $              .AND. (MYCOL .EQ. bl_owner_c_r)) THEN
*                 ..
*                 .. Exchange part of row with proc located above before preforming rotation
*                 ..
                  CALL SEND(1, RCNT, UP, 
     $                 A(bl_index_r_d, bl_index_c_r), LDA,	
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)
                  CALL RECVTMP(1, RCNT, UP, TMPARR, LDTMP, 
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)
                  CALL DROT(RCNT, TMPARR, LDTMP, 
     $                 A(bl_index_r_d, bl_index_c_l), LDA, 
     $                 CSROW(1, K - 1), CSROW(2, K - 1))
               ENDIF
               IF ((MYROW .EQ. bl_owner_r_u) .AND. 
     $              (MYCOL .EQ. bl_owner_c_l)) THEN
*                 ..
*                 .. Exchange part of row with proc located below before preforming rotation
*                 ..
                  CALL RECVTMP(1, LCNT, DOWN, TMPARR, LDTMP, 
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)
                  CALL SEND(1, LCNT, DOWN, 
     $                 A(bl_index_r_u, bl_index_c_l), LDA,	
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)
                  CALL DROT(LCNT, A(bl_index_r_u, bl_index_c_l), LDA,	
     $                 TMPARR, LDTMP, 
     $                 CSROW(1, K - 1), CSROW(2, K - 1))
               ENDIF
               IF ((RCNT .GT. 0) .AND. (MYROW .EQ. bl_owner_r_u) 
     $              .AND. (MYCOL .EQ. bl_owner_c_r)) THEN
*                 ..
*                 .. Exchange part of row with proc located below before preforming rotation
*                 ..
                  CALL RECVTMP(1, RCNT, DOWN, TMPARR, LDTMP,
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)
                  CALL SEND(1, RCNT, DOWN, 
     $                 A(bl_index_r_u, bl_index_c_r), LDA,
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)
                  CALL DROT(RCNT, A(bl_index_r_u, bl_index_c_r), LDA,	
     $                 TMPARR, LDTMP, 
     $                 CSROW(1, K - 1), CSROW(2, K - 1))
               ENDIF
            ENDIF
 60      CONTINUE
 50   CONTINUE
      RETURN
      END
      


      SUBROUTINE UPDQ(IHI, ILO, N, JC, NB, NOJ, Q, LDQ, CS, LDC)
      IMPLICIT NONE
      INTEGER N, LDQ, JC, LDC, NB, IHI, ILO
      DOUBLE PRECISION Q(LDQ, *), CS(LDC,*)
      INTEGER IND, I, K, NOJ, LJ, LNOJ
      INTEGER  GJ, IX, NOIT, II
*
*     Number of active bulges (j columns to annihilate)
*
      LNOJ = 0
      NOIT = (IHI-ILO+1)/NB + NOJ
C      NOIT = N/NB + NOJ
*
*     Global index if first column to annihilate
      GJ = JC - NOJ + 2
      DO 10 IX = 1, NOIT
*
*        Increase the number of active bulges (up to max)
         IF (LNOJ .LT. NOJ) THEN
            LNOJ = LNOJ + 1
         ENDIF
*
*        Local J
         DO 25 II = 1, N, N
            LJ = GJ
            DO 20 I = 1, LNOJ
               IF (LJ .LT. IHI) THEN
                  IND = (I-1)*2+1
                  DO 30 K = MIN(IHI, LJ + NB -1), LJ+1, -1
                     CALL DROT( MIN(N, N-II+1), Q( II, K-1  ), 1,
     $                    Q( II, K ), 1, CS(IND,K-1), CS(IND+1,K-1) )
30                CONTINUE
               ENDIF
*
*              Next bulge
               LJ = LJ - (NB - 1)
20          CONTINUE
25       CONTINUE
         GJ = GJ + NB
10    CONTINUE
      RETURN
      END

      SUBROUTINE pUPDQ(IHI, ILO, N, JC, NB, NOJ, Q, LDQ, CS, LDC,
     $     NPROW, NPCOL, MYROW, MYCOL, TMPARR, LDTMP, CTX)
      
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER		IHI, ILO, N, JC, NB, NOJ, LDQ, LDC
      INTEGER		NPROW, NPCOL, MYROW, MYCOL, LDTMP, CTX
      
*     .. Array Arguments ..
      DOUBLE PRECISION	Q(LDQ, *), CS(LDC,*), TMPARR(LDTMP)
      
*     .. Parameters ..
      INTEGER		DOWN, UP, RIGHT, LEFT
      PARAMETER		( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
      
*     .. Local scalars ..
      INTEGER		IND, I, J, LJ, LNOJ
      INTEGER		GJ, IX, NOIT, II, RIGHT_K, LEFT_K
      INTEGER		bl_owner_r, bl_owner_cl, bl_owner_cr
      INTEGER		I_c, bl_index_r, bl_index_c
      INTEGER		t_I_c, t_bl_index_c	
      INTEGER		NNB_L, NNB_R, NNB_D
      
*     .. External Functions ..
      EXTERNAL		INDXG2L, INDXG2P
      INTEGER		INDXG2L, INDXG2P
      

	LOGICAL		DEBUG1
*     ..
*     .. Start execution
*     ..
      
	DEBUG1 = .FALSE.
*	DEBUG1 = .TRUE.
*     ..
*     .. Number of active bulges (j columns to annihilate)
*     ..
      LNOJ = 0
      NOIT = (IHI-ILO+1)/NB + NOJ
*     ..
*     .. Global index if first column to annihilate
*     ..
      GJ = JC - NOJ + 2

*     ..
*     .. If NPCOL = 1 then NPROW > 1
*     ..
      IF (NPCOL .EQ. 1) THEN
         
         DO 110 IX = 1, NOIT
*           ..
*           .. Increase the number of active bulges (up to max)
*           ..
            IF (LNOJ .LT. NOJ) THEN
               LNOJ = LNOJ + 1
            ENDIF
*           ..
*           .. Local J
*           ..
            DO 125 II = 1, N, N
               LJ = GJ
               DO 120 I = 1, LNOJ
                  IF (LJ .LT. IHI) THEN
                     IND = (I - 1) * 2 + 1
*                    .. Start column
                     RIGHT_K = MIN(IHI, LJ + NB - 1)
*                    .. End column
                     LEFT_K = LJ
                     

                     NNB_L = MIN(NB - MOD(LEFT_K - 1, NB),IHI-LEFT_K+1)
                     NNB_R = MOD(RIGHT_K, NB)


                     IF (NNB_L.EQ.(IHI-LEFT_K+1)) THEN
                        NNB_R = 0
                     END IF
                                                                                        
*                    ..
*                    .. Update the NNB_L + NNB_R columns on row II to N
*                    ..
                     DO 130 J = II, MIN(N, N - II + 1), NB
                        bl_owner_r = INDXG2P(J, NB, 0, 0, NPROW)
                        bl_index_r = INDXG2L(J, NB, 0, 0, NPROW)			
                        NNB_D = MIN(NB, MIN(N, N - II + 1) - J + 1)
                        IF (MYROW .EQ. bl_owner_r) THEN
*                       ..
*                       .. Preform column rotation 
*                       ..
                           I_c = LEFT_K + NNB_L-1
                           bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
                           IF ((NNB_R .GT. 0) .AND. 
     $                          (NNB_L .NE. NB)) THEN
                              t_I_c = RIGHT_K - NNB_R + 1
                              t_bl_index_c = INDXG2L(t_I_c, NB, 
     $                             0, 0, NPCOL)
			      IF (DEBUG1) THEN
		   	         WRITE(*,*)'PUPDQ 1 NNB_R=',NNB_R,
     $                             'NNB_D=',NNB_D,'t_I_c=',t_I_c,
     $                             'IND=',IND,'I_c=',I_c,'J=',J,
     $                             MYROW,MYCOL
			      END IF
                              CALL cRot2(NNB_R, NNB_D, 
     $                             Q(bl_index_r, t_bl_index_c), LDQ,
     $                             CS(IND, t_I_c), LDC)
                              CALL DROT(NNB_D, 
     $                             Q(bl_index_r, bl_index_c), 1, 
     $                             Q(bl_index_r, t_bl_index_c), 1,
     $                             CS(IND, I_c), CS(IND + 1, I_c))
                           END IF
                           I_c = LEFT_K
                           bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
 			   IF (DEBUG1) THEN
				WRITE(*,*)'PUPDQ 2 NNB_L=',NNB_L,
     $                             'NNB_D=',NNB_D,'I_c=',I_c,
     $                             'IND=',IND,'J=',J,
     $                             MYROW,MYCOL
			   END IF
                           
                           CALL cRot2(NNB_L, NNB_D, 
     $                          Q(bl_index_r, bl_index_c), LDQ, 
     $                          CS(IND, I_c), LDC)
                        END IF
 130                 CONTINUE
                  ENDIF
*                 ..
*                 .. Next bulge
*                 ..
                  LJ = LJ - (NB - 1)
 120           CONTINUE
 125        CONTINUE
            GJ = GJ + NB
 110     CONTINUE
         
         RETURN
      END IF

*     ..
*     .. Else we have a grid of minimum size 1x2
*     ..
      
      DO 10 IX = 1, NOIT
*        ..
*        .. Increase the number of active bulges (up to max)
*        ..
         IF (LNOJ .LT. NOJ) THEN
            LNOJ = LNOJ + 1
         ENDIF
*        ..
*        .. Local J
*        ..
         DO 25 II = 1, N, N
            LJ = GJ
		  IF (DEBUG1) THEN
			WRITE(*,*)'Debug 1 LJ=',LJ, 'II=',II,MYROW,MYCOL
	      END IF
            DO 20 I = 1, LNOJ
               IF (LJ .LT. IHI) THEN
                  IND = (I - 1) * 2 + 1
*                 .. Start column
                  RIGHT_K = MIN(IHI, LJ + NB - 1)
*                 .. End column
                  LEFT_K = LJ
                  
                  bl_owner_cl = INDXG2P(LEFT_K, NB, 0, 0, NPCOL)
                  bl_owner_cr = INDXG2P(RIGHT_K, NB, 0, 0, NPCOL)
                  
                  NNB_L = MIN(NB - MOD(LEFT_K-1,NB),IHI-LEFT_K+1)
                  NNB_R = MOD(RIGHT_K, NB)

                  IF (bl_owner_cl .EQ. bl_owner_cr) THEN
                     NNB_R = 0
                  ENDIF

			IF (DEBUG1) THEN
				WRITE(*,*)'Debug 2 I=',I, IND, RIGHT_K, LEFT_K, 
     $				NNB_L, NNB_R,MYROW,MYCOL				 
			END IF
*                 ..
*                 .. Update the NNB_L + NNB_R columns on row II to N
*                 ..
                  DO 30 J = II, MIN(N, N - II + 1), NB
                     bl_owner_r = INDXG2P(J, NB, 0, 0, NPROW)
                     bl_index_r = INDXG2L(J, NB, 0, 0, NPROW)				
                     NNB_D = MIN(NB, MIN(N, N - II + 1) - J + 1)
                     IF ((MYROW .EQ. bl_owner_r) .AND. 
     $                    (MYCOL .EQ. bl_owner_cl)) THEN
*                       ..
*                       .. Exchange elements with proc to the right
*                       ..
                        I_c = LEFT_K + NNB_L-1

                        bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
                        IF ((NNB_R .GT. 0) .AND. 
     $                       (NNB_L .NE. NB)) THEN

			 IF (DEBUG1) THEN
				WRITE(*,*)'Debug 3 I=',I,J,I_c, 
     $				NNB_D, bl_index_c,MYROW,MYCOL
			 END IF
                           CALL Send(NNB_D, 1, RIGHT, 
     $                          Q(bl_index_r, bl_index_c),
     $                          LDQ, MYROW, MYCOL, NPCOL, NPROW, CTX)

                           CALL Recv(NNB_D, 1, RIGHT, TMPARR, 
     $                          LDTMP, MYROW, MYCOL, NPCOL, NPROW, CTX)

                           CALL DROT(NNB_D, Q(bl_index_r, bl_index_c), 
     $                          1,TMPARR, 1,
     $                          CS(IND, I_c), CS(IND + 1, I_c))
                        END IF
                        I_c = LEFT_K
                        bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

		  IF (DEBUG1) THEN
			WRITE(*,*)'Debug 4 I=',I,J,I_c, 
     $			NNB_D, bl_index_c,MYROW,MYCOL
		  END IF
*                       ..
*                       .. Preform column rotation 
*                       ..
                        CALL cRot2(NNB_L, NNB_D, 
     $                       Q(bl_index_r, bl_index_c), LDQ,
     $                       CS(IND, I_c), LDC)
                     END IF
                     IF ((MYROW .EQ. bl_owner_r) .AND. 
     $                    (MYCOL .EQ. bl_owner_cr) .AND. 
     $                    (MYCOL .NE. bl_owner_cl)) THEN
                        I_c = RIGHT_K -NNB_R + 1
                        bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
		  IF (DEBUG1) THEN
			WRITE(*,*)'Debug 5 I=',I,J,I_c, 
     $			NNB_D, bl_index_c,MYROW,MYCOL
		  END IF
*                       ..
*                       .. Preform column rotation 
*                       ..
                        IF (NNB_R .GT. 1) THEN
                           CALL cRot2(NNB_R, NNB_D, 
     $                          Q(bl_index_r, bl_index_c), LDQ, 
     $                          CS(IND, I_c), LDC)
                        END IF
*                       ..
*                       .. Exchange elements with left proc
*                       .. 
                        CALL Send(NNB_D, 1, LEFT, 
     $                       Q(bl_index_r, bl_index_c),
     $                       LDQ, MYROW, MYCOL, NPCOL, NPROW, CTX)
                        
                        CALL Recv(NNB_D, 1, LEFT, TMPARR, 
     $                       LDTMP, MYROW, MYCOL, NPCOL, NPROW, CTX)

                        CALL DROT(NNB_D,
     $                       TMPARR, 1, Q(bl_index_r, bl_index_c), 1, 
     $                       CS(IND, I_c - 1), CS(IND + 1, I_c - 1))
                        
                     END IF
 30               CONTINUE
               ENDIF
*              ..
*              .. Next bulge
*              ..
               LJ = LJ - (NB - 1)
 20         CONTINUE
 25      CONTINUE
         GJ = GJ + NB
 10   CONTINUE
      RETURN
      END



      SUBROUTINE UPDZ(IHI, ILO, N, JC, NB, NOJ, Z, LDZ, CS, LDC)

      IMPLICIT NONE     
      INTEGER N, LDZ, JC, LDC, NB, IHI, ILO
      DOUBLE PRECISION Z(LDZ, *), CS(LDC,*)
      INTEGER IND, I, K, NOJ, LJ, LNOJ
      INTEGER  GJ, IX, NOIT, II
*
*     Number of active bulges (j columns to annihilate)
*     
      LNOJ = 0
      NOIT = (IHI-ILO+1)/NB + NOJ
*
*     Global index if first column to annihilate
      GJ = JC - NOJ + 2
      DO 10 IX = 1, NOIT
*
*        Increase the number of active bulges (up to max)
         IF (LNOJ .LT. NOJ) THEN
            LNOJ = LNOJ + 1
         ENDIF
*
*        Local J
         DO 25 II = 1, N, N
            LJ = GJ
            DO 20 I = 1, LNOJ
               IF (LJ .LT. IHI) THEN
                  IND = (I-1)*2+1
                  DO 30 K = MIN(IHI, LJ + NB -1), LJ + 1, -1
                     CALL DROT( MIN(N, N-II+1), Z( II, K  ), 1,
     $                    Z( II, K - 1 ), 1, CS(IND,K-1), CS(IND+1,K-1))
 30               CONTINUE
               ENDIF
*
*           Next bulge
               LJ = LJ - (NB-1)
 20         CONTINUE
 25      CONTINUE
         GJ = GJ + NB
10    CONTINUE
      END

      SUBROUTINE pUPDZ(IHI, ILO, N, JC, NB, NOJ, Z, LDZ, CS, LDC,
     $     NPROW, NPCOL, MYROW, MYCOL, TMPARR, LDTMP, CTX)
      IMPLICIT NONE
      
*     .. Scalar Arguments ..
      INTEGER           N, LDZ, JC, LDC, NB, IHI, ILO, LDTMP
      INTEGER           NPROW, NOJ, NPCOL, MYROW, MYCOL, CTX
*     .. Array Arguments ..
      DOUBLE PRECISION  Z(LDZ, *), CS(LDC,*), TMPARR(LDTMP, *)
      
*     .. Parameters ..
      INTEGER           DOWN, UP, RIGHT, LEFT
      PARAMETER         ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
      
*     .. Local scalars ..
      INTEGER           IND, I, J, LJ, LNOJ
      INTEGER           GJ, IX, NOIT, II, RIGHT_K, LEFT_K
      INTEGER           bl_owner_r, bl_owner_cl, bl_owner_cr
      INTEGER           I_c, bl_index_r, bl_index_c
      INTEGER           t_I_c, t_bl_index_c
      INTEGER           NNB_L, NNB_R, NNB_D
      LOGICAL           DEBUG1
      
*     .. External Functions ..
      EXTERNAL          INDXG2L, INDXG2P
      INTEGER           INDXG2L, INDXG2P
      
*     ..
*     .. Begin execution
*     .. 
*     ..
*     .. Number of active bulges (j columns to annihilate)
*     ..
      DEBUG1 = .FALSE.
*      DEBUG1 = .TRUE.
      LNOJ = 0


      NOIT = (IHI-ILO+1)/NB + NOJ
*     ..
*     .. Global index if first column to annihilate
*     ..
      GJ = JC - NOJ + 2
*     ..
*     .. If NPCOL = 1 then NPROW > 1
*     ..
      IF (NPCOL .EQ. 1) THEN
         DO 110 IX = 1, NOIT
*           ..
*           .. Increase the number of active bulges (up to max)
*           ..
            IF (LNOJ .LT. NOJ) THEN
               LNOJ = LNOJ + 1
            ENDIF
*           ..
*           .. Local J
*           ..
            DO 125 II = 1, N, N
               LJ = GJ
               DO 120 I = 1, LNOJ
                  IF (LJ .LT. IHI) THEN
                     IND = (I - 1) * 2 + 1                     
*                    .. Start column
                     RIGHT_K = MIN(IHI, LJ + NB - 1)
*                    .. End column
                     LEFT_K = LJ
                     
                     NNB_L = MIN(NB - MOD(LEFT_K-1,NB),IHI-LEFT_K+1)
                     NNB_R = MOD(RIGHT_K, NB)

                     IF (NNB_L.EQ.(IHI-LEFT_K+1)) THEN
                        NNB_R = 0
                     ENDIF
*                    ..
*                    .. Update the NNB_L + NNB_R columns on row II to N
*                    ..
                     DO 130 J = II, MIN(N, N - II + 1), NB
                        bl_index_r = INDXG2L(J, NB, 0, 0, NPROW)
                        bl_owner_r = INDXG2P(J, NB, 0, 0, NPROW)
                        
                        NNB_D = MIN(NB, MIN(N, N - II + 1) - J + 1)

                        
                        IF (MYROW .EQ. bl_owner_r) THEN
*                          ..
*                          .. Preform column rotation
*                          ..
                           I_c = LEFT_K + NNB_L-1
                           bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
                           IF ((NNB_R .GT. 0) .AND. 
     $                          (NNB_L .NE. NB)) THEN
                              t_I_c = RIGHT_K-NNB_R+1

                              t_bl_index_c = INDXG2L(t_I_c, NB, 
     $                             0, 0, NPCOL)
                              IF (DEBUG1) THEN
                              WRITE(*,*)'PUPDZ 1 NNB_R=',NNB_R,
     $                               'NNB_D=',NNB_D, 'IND=',IND,
     $                               't_I_c=',t_I_c,
     $                               'I_c =',I_c,'J=',J,
     $                               MYROW,MYCOL
                              END IF
                              CALL cRot(NNB_R, NNB_D, 
     $                             Z(bl_index_r, t_bl_index_c), LDZ,
     $                             CS(IND, t_I_c), LDC)
                              CALL DROT(NNB_D, 
     $                             Z(bl_index_r, t_bl_index_c), 1,
     $                             Z(bl_index_r, bl_index_c), 1,
     $                             CS(IND, I_c), CS(IND + 1, I_c))
                           END IF
                           I_c = LEFT_K
                           bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
*                          ..
*                          .. Preform column rotation 
*                          ..

                           IF (DEBUG1) THEN
                               WRITE(*,*)'PUPDZ 2 NNB_L=',NNB_L,
     $                         'NNB_D=',NNB_D,'I_c=',I_c,'IND=',IND,
     $                         'J=',J,MYROW,MYCOL
                           END IF
                           CALL cRot(NNB_L, NNB_D, 
     $                          Z(bl_index_r, bl_index_c), LDZ,
     $                          CS(IND, I_c), LDC)
                           
                        ENDIF
 130                 CONTINUE
                  ENDIF
*                 ..
*                 .. Next bulge
*                 ..
                  LJ = LJ - (NB - 1)
 120           CONTINUE
 125        CONTINUE
            GJ = GJ + NB
 110     CONTINUE

         RETURN
      END IF

*     ..
*     .. Else we have a grid of minimum size 1x2
*     ..

      DO 10 IX = 1, NOIT
*        ..
*        .. Increase the number of active bulges (up to max)
*        .. 
         IF (LNOJ .LT. NOJ) THEN
            LNOJ = LNOJ + 1
         ENDIF
*        ..
*        .. Local J
*        ..
         DO 25 II = 1, N, N
            LJ = GJ
            DO 20 I = 1, LNOJ
               IF (LJ .LT. IHI) THEN
                  IND = (I - 1) * 2 + 1            
*                 .. Start column
                  RIGHT_K = MIN(IHI, LJ + NB - 1)
*                 .. End column
                  LEFT_K = LJ

                  bl_owner_cl = INDXG2P(LEFT_K, NB, 0, 0, NPCOL)
                  bl_owner_cr = INDXG2P(RIGHT_K, NB, 0, 0, NPCOL)
                  
                  NNB_L = MIN(NB-MOD(LEFT_K-1,NB),IHI-LEFT_K+1)
                  NNB_R = MOD(RIGHT_K, NB)

                  IF (bl_owner_cl .EQ. bl_owner_cr) THEN
                     NNB_R = 0
                  ENDIF
*                 ..
*                 .. Update the NNB_L + NNB_R columns on row II to N
*                 ..
                  DO 30 J = II, MIN(N, N - II + 1), NB
                     bl_index_r = INDXG2L(J, NB, 0, 0, NPROW)
                     bl_owner_r = INDXG2P(J, NB, 0, 0, NPROW)
                     
                     NNB_D = MIN(NB, MIN(N, N - II + 1) - J + 1)



                     IF ((MYROW .EQ. bl_owner_r) .AND. 
     $                    (MYCOL .EQ. bl_owner_cl)) THEN
*                       ..
*                       .. Exchange elements with proc to the right
*                       ..
                        I_c = LEFT_K + NNB_L-1

                        bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
                        IF ((NNB_R .GT. 0) .AND. (NNB_L .NE. NB)) THEN
                           CALL Send(NNB_D, 1, RIGHT, 
     $                          Z(bl_index_r, bl_index_c),
     $                          LDZ, MYROW, MYCOL, NPCOL, NPROW, CTX)

                           CALL Recv(NNB_D, 1, RIGHT, TMPARR, 
     $                          LDTMP, MYROW, MYCOL, NPCOL, NPROW, CTX)

                           CALL DROT(NNB_D, TMPARR, 1,
     $                          Z(bl_index_r, bl_index_c), 1,
     $                          CS(IND, I_c), CS(IND + 1, I_c))
                        END IF
                        I_c = LEFT_K
                        bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
*                       ..
*                       .. Preform column rotation 
*                       ..
                        CALL cRot(NNB_L, NNB_D, 
     $                       Z(bl_index_r, bl_index_c), LDZ,
     $                       CS(IND, I_c), LDC)
                        
                     ENDIF
                     IF ((MYROW .EQ. bl_owner_r) .AND. 
     $                    (MYCOL .EQ. bl_owner_cr) .AND. 
     $                    (MYCOL .NE. bl_owner_cl)) THEN
                        
                        I_c = RIGHT_K-NNB_R+1
                        bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

*                       ..
*                       .. Preform column rotation 
*                       ..
                        IF (NNB_R .GT. 1) THEN
                           CALL cRot(NNB_R, NNB_D, 
     $                          Z(bl_index_r, bl_index_c), LDZ,
     $                          CS(IND, I_c), LDC)
                        END IF
*                       ..
*                       .. Exchange elements with left proc
*                       .. 
                        CALL Send(NNB_D, 1, LEFT, 
     $                       Z(bl_index_r, bl_index_c),
     $                       LDZ, MYROW, MYCOL, NPCOL, NPROW, CTX)

                        CALL Recv(NNB_D, 1, LEFT, TMPARR, 
     $                       LDTMP, MYROW, MYCOL, NPCOL, NPROW, CTX)

                        CALL DROT(NNB_D,
     $                       Z(bl_index_r, bl_index_c), 1, TMPARR, 1, 
     $                       CS(IND, I_c - 1), CS(IND + 1, I_c - 1))
                     ENDIF
 30               CONTINUE
               ENDIF
*              ..
*              .. Next bulge
*              ..
               LJ = LJ - (NB - 1)
 20         CONTINUE
 25      CONTINUE
         GJ = GJ + NB
 10   CONTINUE
      RETURN
      END

      SUBROUTINE KILLCOL(N, A, LDA, CSROW, LDC)
      IMPLICIT NONE

*     .. Scalar Arguments ..
      INTEGER            N, LDA, LDC

*     .. Array Arguments ..
      DOUBLE PRECISION   A(*)
      DOUBLE PRECISION   CSROW( LDC, *)
      INTEGER K

*     .. Parameters ..
      DOUBLE PRECISION   ZERO, TEMP
      PARAMETER          ( ZERO = 0.0D+0 )

      DO 20 K = N - 1, 1, -1
         TEMP = A(K)
         CALL DLARTG(TEMP, A(K + 1), CSROW(1, K),
     $               CSROW(2, K), A(K))
         A(K + 1) = ZERO
20    CONTINUE
      RETURN
      END

      SUBROUTINE UPD(M, N, A, LDA, NB, J, CSROW, CSCOL, LDC)
      INTEGER M, N, LDA, NB, J, LAST, LDC
      DOUBLE PRECISION A(LDA, *), CSROW(LDC,*), CSCOL(LDC,*)
      INTEGER NNB, JJ, BZ, IND
      PARAMETER (NNB = 15)
C      NNB = NB - 1
      LAST = J
      DO 5 JJ = N, 2, -NNB
*
*     Apply column rotations
         DO 10 K = JJ, MAX(JJ-NNB+1, 2), -1
            CALL DROT( M, A( 1, K  ), 1,
     $           A( 1, K-1 ), 1, CSCOL(1,K-1), CSCOL(2,K-1) )

10       CONTINUE
*
*     Apply row rotations
         IND = JJ-NNB+1
         BZ  = NNB
         IF (JJ .LE. NNB+1) THEN
            IND = 1
            BZ = JJ
         ENDIF
         DO 20 KK = J+NB, M, NB
            DO 30 K = KK, KK-NB+2, -1
               CALL DROT( BZ, A( K - 1,IND ), LDA,
     $              A( K, IND ), LDA, CSROW(1,K-1), CSROW(2,K-1) )
30          CONTINUE
            LAST = KK
20       CONTINUE
5     CONTINUE
      IF (N.EQ.1) THEN
*
*        Apply row rotations
         DO 200 KK = J+NB, M, NB
            DO 300 K = KK, KK-NB+2, -1
               CALL DROT( N, A( K - 1,1 ), LDA,
     $           A( K, 1 ), LDA, CSROW(1,K-1), CSROW(2,K-1) )
300         CONTINUE
            LAST = KK
200      CONTINUE
      ENDIF
*     If it is the last block column update all rows
      DO 40 K = M, LAST+2, -1
         CALL DROT( N, A( K - 1,1 ), LDA,
     $        A( K, 1 ), LDA, CSROW(1,K-1), CSROW(2,K-1) )
 40   CONTINUE
      RETURN 
      END

      SUBROUTINE MVBB(N, A, LDA, CSROW, CSCOL, LDC)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            N, LDA, LDC

*     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA, *)
      DOUBLE PRECISION   CSROW( LDC, *), CSCOL(LDC, *)
      INTEGER K

*     .. Parameters ..
      DOUBLE PRECISION   ZERO, TEMP
      PARAMETER          ( ZERO = 0.0D+0 )

*     .. Apply row rotations to Diagonal block of A
      DO 10 K = N-1, 1, -1
         CALL DROT(N-K + 1 , A( K, K ), LDA,
     $             A( K+1, K), LDA, CSROW(1,K), CSROW(2, K) )
10    CONTINUE
*      ..
*     .. Annihilate B fill-in
      DO 20 K = N, 2, -1
         TEMP = A(K, K)
         CALL DLARTG(TEMP, A(K,K-1), CSCOL(1,K-1),
     $               CSCOL(2,K-1),A(K,K))
         A( K, K-1) = ZERO
         CALL DROT( K - 1, A( 1, K ), 1,
     $              A( 1, K-1 ), 1, CSCOL(1,K-1), CSCOL(2,K-1) )
20    CONTINUE
      RETURN
      END

      SUBROUTINE MVBA(M, N, A, LDA, CSROW, CSCOL, LDC)
      IMPLICIT NONE
      INTEGER            M, N, LDA, LDC
      DOUBLE PRECISION   A(LDA, *)
      DOUBLE PRECISION   CSROW( LDC, *), CSCOL(LDC, *)
      INTEGER K, KK

*     .. Parameters ..
      DOUBLE PRECISION   ZERO, TEMP
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Apply col-roations to diagonal block of A (the part of the
*     .. rotations that not result in fill-in)
      DO 20 K = N-1, M, -1
         CALL DROT(M, A( 1, K+1 ), 1,
     $       A( 1, K ), 1, CSCOL(1,K), CSCOL(2,K))
20    CONTINUE
*     ..
*     .. Apply col rotations to diag block of A and delete fill-in
      KK = 0
      DO 80 K = M - 1, 1, -1
         CALL DROT(M-KK, A( 1, K+1 ), 1,
     $            A( 1, K ), 1, CSCOL(1,K), CSCOL(2,K))
         KK = KK + 1
         TEMP = A(K, K)         
         CALL DLARTG(TEMP,A(K+1,K),CSROW(1,K),
     $        CSROW(2, K), A( K, K))
         A( K+1, K) = ZERO
         CALL DROT(N-K, A( K, K+1 ), LDA, A( K+1, K+1 ),
     $        LDA, CSROW(1,K), CSROW(2,K) )
80    CONTINUE
      RETURN
      END


      SUBROUTINE droRot(N, A, LDA, CSROW, LDC)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            N, LDA, LDC

*     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA, *)
      DOUBLE PRECISION   CSROW( LDC, *)

*     .. Local scalars ..
      INTEGER K

*     ..
*     .. Apply row rotations to Diagonal block of A
*     ..
      DO 10 K = N - 1, 1, -1
         CALL DROT(N-K + 1, A( K, K ), LDA,
     $        A( K + 1, K), LDA, CSROW(1, K), CSROW(2, K))
10    CONTINUE
      RETURN 
      END

      SUBROUTINE drRot(N, A, LDA, CSROW, LDC)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            N, LDA, LDC

*     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA, *)
      DOUBLE PRECISION   CSROW( LDC, *)

*     .. Local scalars ..
      INTEGER K

*     .. Parameters ..
      DOUBLE PRECISION   ZERO, TEMP
      PARAMETER          ( ZERO = 0.0D+0 )

*     ..
*     .. Apply row rotations to diag block of A and delete fill-in
*     ..
      DO 80 K = N - 1, 1, -1
         TEMP = A(K, K)
         CALL DLARTG(TEMP, A(K + 1, K), CSROW(1, K),
     $        CSROW(2, K), A( K, K))
         A(K + 1, K) = ZERO
         CALL DROT(N - K, A(K, K + 1), LDA, A(K + 1, K + 1),
     $        LDA, CSROW(1, K), CSROW(2, K))
80    CONTINUE
      RETURN
      END

      SUBROUTINE rRot(N, M, A, LDA, CSROW, LDC, SKIP)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            N, M, LDA, LDC, SKIP

*     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA, *)
      DOUBLE PRECISION   CSROW( LDC, *)

*     .. Local scalars ..
      INTEGER K

*     ..
*     .. Apply row rotations to Diagonal block of A
*     ..
      DO 10 K = N - 1, 1, -1
         IF (SKIP .NE. K) THEN
            CALL DROT(M, A(K, 1), LDA,
     $           A(K + 1, 1), LDA, CSROW(1, K), CSROW(2, K))
         ENDIF
10    CONTINUE
      RETURN 
      END

      SUBROUTINE dcRot(N, A, LDA, CSCOL, LDC)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            N, LDA, LDC

*     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA, *)
      DOUBLE PRECISION   CSCOL( LDC, *)

*     .. Local scalars ..
      INTEGER K

*     .. Parameters ..
      DOUBLE PRECISION   ZERO, TEMP
      PARAMETER          ( ZERO = 0.0D+0 )

*     ..
*     .. Annihilate fill-in and apply rotation to the rest
*     ..

      DO 20 K = N, 2, -1
         TEMP = A(K, K)
         CALL DLARTG(TEMP, A(K, K - 1), CSCOL(1, K - 1),
     $               CSCOL(2, K - 1), A(K, K))
         A( K, K - 1) = ZERO
         CALL DROT( K - 1, A(1, K ), 1,
     $              A(1, K - 1), 1, CSCOL(1, K - 1), CSCOL(2, K - 1))
20    CONTINUE

      RETURN 
      END


      SUBROUTINE cRot(N, M, A, LDA, CSCOL, LDC)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            M, N, LDA, LDC

*     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA, *)
      DOUBLE PRECISION   CSCOL( LDC, *)

*     .. Local scalars ..
      INTEGER K


*     ..
*     .. Apply colrotations
*     ..
      DO 20 K = N, 2, -1
         CALL DROT( M, A(1, K ), 1,
     $              A(1, K - 1), 1, CSCOL(1, K - 1), CSCOL(2, K - 1))
20    CONTINUE
      RETURN 
      END

      SUBROUTINE cRot2(N, M, A, LDA, CSCOL, LDC)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            M, N, LDA, LDC

*     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA, *)
      DOUBLE PRECISION   CSCOL( LDC, *)

*     .. Local scalars ..
      INTEGER K


*     ..
*     .. Apply colrotations
*     ..
      DO 20 K = N, 2, -1
         CALL DROT( M, A(1, K-1 ), 1,
     $              A(1, K ), 1, CSCOL(1, K - 1), CSCOL(2, K - 1))
20    CONTINUE
        
      RETURN 
      END
      SUBROUTINE dcrRot(N, M, A, LDA, CSROW, CSCOL, LDC)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            N, M, LDA, LDC

*     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA, *)
      DOUBLE PRECISION   CSCOL( LDC, *), CSROW( LDC, *)

*     .. Local scalars ..
      INTEGER K

*     .. Parameters ..
      DOUBLE PRECISION   ZERO, TEMP
      PARAMETER          ( ZERO = 0.0D+0 )

*     ..
*     .. apply colrotation
*     ..
      DO 10 K = N, 2, -1
         CALL DROT(K, A(1, K), 1,
     $              A(1, K - 1), 1, CSCOL(1, K - 1), CSCOL(2, K - 1))


         TEMP = A(K - 1, K - 1)         
         CALL DLARTG(TEMP, A(K, K - 1), CSROW(1, K - 1),
     $        CSROW(2, K - 1), A( K - 1, K - 1))
         A(K, K - 1) = ZERO

         CALL DROT(M - K + 1, A(K - 1, K), LDA, A(K, K),
     $        LDA, CSROW(1, K - 1), CSROW(2, K - 1))
10    CONTINUE
      RETURN 
      END


      SUBROUTINE pMVBB(B,LDB,CSROW,CSCOL,LDC,I,JC,NB,NB1,N,NPROW,NPCOL,
     $     IND,TMPARR,LDTMP,MYROW,MYCOL,CTX)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            I, JC, NB, NB1, N, NPROW, NPCOL, LDB, LDC, IND 
      INTEGER            LDTMP,MYROW, MYCOL, CTX

*     .. Array Arguments ..
      DOUBLE PRECISION   B(LDB,*)
      DOUBLE PRECISION   CSCOL(LDC,*), CSROW(LDC,*)
      DOUBLE PRECISION   TMPARR(LDTMP)

*     .. Parameters ..
      DOUBLE PRECISION   ZERO, TEMP
      PARAMETER          ( ZERO = 0.0D+0 )
      INTEGER            DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
      INTEGER            ROWWISE, COLUMNWISE
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2)

*     .. Local scalars ..
      INTEGER            bl_owner_r, bl_owner_c, bl_right_c, bl_bottom_r
      INTEGER            I_c, I_r, bl_index_r, bl_index_c, NNB1, NNB2
      INTEGER            NNB3
      INTEGER            t_I_c, t_I_r, t_bl_index_r, t_bl_index_c
*     .. External Functions ..
      EXTERNAL           INDXG2L, INDXG2P
      INTEGER            INDXG2L,INDXG2P

*     ..
*     .. Start execution
*     ..

      NNB1 = MOD(I - 1, NB) - (NB - NB1)
      NNB2 = NB - MOD(I - 1, NB)

      IF (NNB1 .LT. 0) THEN
		NNB2 = NNB2 + NNB1
      ENDIF


*     ..
*     .. If NPROW = 1 then NPCOL > 1 and vice versa
*     ..
      IF (NPROW .EQ. 1) THEN
         bl_owner_r = INDXG2P(I, NB, 0, 0, NPROW)
         bl_owner_c = INDXG2P(I, NB, 0, 0, NPCOL)
         bl_right_c = INDXG2P(I + NB1 - 1, NB, 0, 0, NPCOL)

*        ..
*        .. The processor who owns B(I,I)
*        ..
         IF ((MYROW .EQ. bl_owner_r) .AND. (MYCOL .EQ. bl_owner_c)) THEN
            I_r = I + NNB2 - 1 
            I_c = I_r
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            IF ((MOD(I - 1, NB) .NE. 0) .AND. (NNB1 .GT. 0)) THEN
               t_I_r = I + NNB2
               t_bl_index_r = INDXG2L(t_I_r, NB, 0, 0, NPROW)
*              ..
*              .. Preform rowrotation
*              ..
               CALL DROT(1, B(bl_index_r, bl_index_c), LDB,
     $              B(t_bl_index_r, bl_index_c), LDB,
     $              CSROW(IND, t_I_r - 1), CSROW(IND + 1, t_I_r - 1))

               CALL Send(1, 1, RIGHT, B(t_bl_index_r, bl_index_c),
     $              LDB, MYROW, MYCOL, NPCOL, NPROW, CTX)

               CALL Recv(1, 1, RIGHT, TMPARR, LDTMP, MYROW, MYCOL, 
     $              NPCOL, NPROW, CTX)

               TEMP = TMPARR(1)
               CALL DLARTG (TEMP, B(t_bl_index_r, bl_index_c), 
     $              CSCOL(IND, I_c),
     $              CSCOL(IND + 1, I_c), TMPARR)

               B(t_bl_index_r, bl_index_c) = ZERO
            ENDIF
*           ..
*           .. Preform rowrotation
*           ..
            I_r = I
            I_c = I
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            CALL droRot(NNB2, B(bl_index_r, bl_index_c),
     $           LDB, CSROW(IND, I_r), LDC)

*           ..
*           .. Exchange elemets with right to let that block be colrotated
*           ..
            I_c = I + NNB2 - 1
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            IF ((MOD(I - 1, NB) .NE. 0) .AND. (NNB1 .GT. 0)) THEN
               CALL Send(NNB2, 1, RIGHT, B(bl_index_r, bl_index_c),
     $              LDB, MYROW, MYCOL, NPCOL, NPROW, CTX)

               CALL Recv(NNB2, 1, RIGHT, TMPARR, LDTMP,
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)

               CALL DROT(NNB2, TMPARR, 1,
     $              B(bl_index_r, bl_index_c), 1,
     $              CSCOL(IND, I_c), CSCOL(IND + 1, I_c))
            ENDIF
*           ..
*           .. Calculate elimination elements and preform colrotation
*           ..
            I_c = I
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

*           ..
*           .. diagonal column rot with elimination
*           ..
            CALL dcRot(NNB2, B(bl_index_r, bl_index_c),
     $           LDB, CSCOL(IND, I_c), LDC)
            
         ENDIF

*        ..
*        .. The processor to the right of the one who owns B(I,I)
*        ..         
         IF ((MYROW .EQ. bl_owner_r) .AND. (MYCOL .EQ. bl_right_c).AND.
     $        (MYCOL .NE. bl_owner_c)) THEN
            I_c = I + NNB2
            I_r = I_c - 1
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)                  
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)                  
*           ..
*           .. Preform rowrotation
*           ..
            t_I_r = I + NNB2
            t_bl_index_r = INDXG2L(t_I_r, NB, 0, 0, NPROW)
                  
            CALL droRot(NNB1, B(t_bl_index_r, bl_index_c), LDB,
     $           CSROW(IND, t_I_r), LDC)

            CALL DROT(NNB1, B(bl_index_r, bl_index_c), LDB,
     $           B(t_bl_index_r, bl_index_c), LDB,
     $           CSROW(IND, t_I_r - 1), CSROW(IND + 1, t_I_r - 1))
            I_r = I
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
         
            CALL rRot(NNB2, NNB1, B(bl_index_r, bl_index_c),
     $           LDB, CSROW(IND, I_r), LDC, -1)

*           ..
*           .. Preform colrotation
*           ..
            CALL dcRot(NNB1, B(t_bl_index_r, bl_index_c), LDB,
     $           CSCOL(IND, I_c), LDC)
       
            CALL Send(1, 1, LEFT, B(t_bl_index_r, bl_index_c), LDB,
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)

            CALL Recv(1, 1, LEFT, TMPARR, LDTMP, MYROW, MYCOL, 
     $           NPCOL, NPROW, CTX)

            TEMP = B(t_bl_index_r, bl_index_c)
            CALL DLARTG (TEMP, TMPARR, CSCOL(IND, I_c - 1),
     $           CSCOL(IND+1, I_c - 1), B(t_bl_index_r, bl_index_c))

            IF (NNB1 .GT. 1) THEN
               CALL cRot(NNB1, NNB2, B(bl_index_r, bl_index_c),
     $              LDB, CSCOL(IND, I_c), LDC)
            ENDIF

            CALL Send(NNB2, 1, LEFT, B(bl_index_r, bl_index_c), LDB,
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)

            CALL Recv(NNB2, 1, LEFT, TMPARR, LDTMP,
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)
                  
            CALL DROT(NNB2, B(bl_index_r, bl_index_c), 1, TMPARR, 1,
     $           CSCOL(IND, I_c - 1), CSCOL(IND + 1, I_c - 1))
         
         ENDIF
         GOTO 999
      ENDIF
      
      IF (NPCOL .EQ. 1) THEN
         bl_owner_r = INDXG2P(I, NB, 0, 0, NPROW)
         bl_owner_c = INDXG2P(I, NB, 0, 0, NPCOL)
         bl_bottom_r = INDXG2P(I + NB1 - 1, NB, 0, 0, NPROW)
*        ..
*        .. The processor who owns B(I,I)
*        ..         
         IF ((MYROW .EQ. bl_owner_r) .AND. (MYCOL .EQ. bl_owner_c)) THEN
            I_r = I + NNB2 - 1 
            I_c = I_r
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            IF ((MOD(I - 1, NB) .NE. 0) .AND. (NNB1 .GT. 0)) THEN

               t_I_c = I + NNB2
               t_bl_index_c = INDXG2L(t_I_c, NB, 0, 0, NPCOL)                  

*              ..
*              .. Exchange elements with down for rowupdate
*              ..
               CALL Send(1, NNB1, DOWN, B(bl_index_r, t_bl_index_c),
     $              LDB, MYROW, MYCOL, NPCOL, NPROW, CTX)
               CALL Recvtmp(1, NNB1, DOWN, TMPARR,
     $              LDTMP, MYROW, MYCOL, NPCOL, NPROW, CTX)
               CALL DROT(NNB1, B(bl_index_r, t_bl_index_c), LDB,
     $              TMPARR, 1,
     $              CSROW(IND, I_r), CSROW(IND + 1, I_r))
*              ..
*              .. Preform rowrotation
*              .. 
               t_I_r = I
               t_bl_index_r = INDXG2L(t_I_r, NB, 0, 0, NPROW)
         
               CALL rRot(NNB2, NNB1, B(t_bl_index_r, t_bl_index_c),
     $              LDB, CSROW(IND, t_I_r), LDC, -1)

*              ..
*              .. Exchange elements with down for rowupdate
*              ..
               CALL Send(1, 1, DOWN, B(bl_index_r, bl_index_c), LDB,
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)

               CALL Recv(1, 1, DOWN, TMPARR, LDTMP,
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)
               CALL DROT(1, B(bl_index_r, bl_index_c), LDB,
     $              TMPARR, 1,
     $              CSROW(IND, I_r), CSROW(IND + 1, I_r))
            ENDIF
*           ..
*           .. Preform rowrotation
*           ..
            I_r = I
            I_c = I
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            CALL droRot(NNB2, B(bl_index_r, bl_index_c),
     $           LDB, CSROW(IND, I_r), LDC)

*           ..
*           .. Preform colrotation
*           ..
            I_c = I + NNB2 - 1
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            IF ((MOD(I - 1, NB) .NE. 0) .AND. (NNB1 .GT. 0)) THEN
*              ..
*              .. Recv cs-colvalues
*              .. 
               CALL BRecv(2, NNB1, COLUMNWISE, CSCOL(IND, t_I_c - 1), 
     $              LDC, bl_bottom_r, bl_owner_c, CTX)
               CALL cRot(NNB1, NNB2, B(bl_index_r, t_bl_index_c),
     $              LDB, CSCOL(IND, t_I_c), LDC)
               CALL DROT(NNB2, B(bl_index_r, t_bl_index_c), 1, 
     $              B(bl_index_r, bl_index_c), 1,
     $              CSCOL(IND, t_I_c - 1), CSCOL(IND + 1, t_I_c - 1))
            ENDIF

            I_c = I
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            CALL dcRot(NNB2, B(bl_index_r, bl_index_c),
     $           LDB, CSCOL(IND, I_c), LDC)

*           ..
*           .. BCast cs-colvalues along column
*           .. 
            IF (NNB2 .GT. 1) THEN
               CALL BSend(2, NNB2 - 1, COLUMNWISE, CSCOL(IND, I),
     $              LDC, CTX)
            ENDIF

         ENDIF

*        ..
*        .. The processor below the one who owns B(I,I)
*        ..
         IF ((MYROW .EQ. bl_bottom_r) .AND. (MYCOL .EQ. bl_owner_c).AND.
     $        (MYROW .NE. bl_owner_r)) THEN
            I_r = I + NNB2
            I_c = I_r - 1
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)            
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

*           ..
*           .. Preform rowrotation
*           ..
            t_I_c = I + NNB2
            t_bl_index_c = INDXG2L(t_I_c, NB, 0, 0, NPCOL)
                  
            CALL droRot(NNB1, B(bl_index_r, t_bl_index_c), LDB,
     $           CSROW(IND, I_r), LDC)
            CALL Send(1, NNB1, UP, B(bl_index_r, t_bl_index_c), LDB,
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)
            CALL Recvtmp(1, NNB1, UP, TMPARR, LDTMP,
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)
            CALL DROT(NNB1, TMPARR, 1,
     $           B(bl_index_r, t_bl_index_c), LDB,
     $           CSROW(IND, I_r - 1), CSROW(IND + 1, I_r - 1))


            CALL Send(1, 1, UP, B(bl_index_r, bl_index_c), LDB,
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)
            CALL Recv(1, 1, UP, TMPARR, LDTMP, MYROW, MYCOL, 
     $           NPCOL, NPROW, CTX)
            CALL DROT(1, TMPARR, 1,
     $           B(bl_index_r, bl_index_c), LDB,
     $           CSROW(IND, I_r - 1), CSROW(IND + 1, I_r - 1))         

*           ..
*           .. Calculate elimination elements and preform colrotation
*           ..
            CALL dcRot(NNB1, B(bl_index_r, t_bl_index_c), LDB,
     $           CSCOL(IND, t_I_c), LDC)

            TEMP = B(bl_index_r, t_bl_index_c)
            CALL DLARTG (TEMP, B(bl_index_r, bl_index_c), 
     $           CSCOL(IND, t_I_c - 1),
     $           CSCOL(IND+1, t_I_c - 1), B(bl_index_r, t_bl_index_c))
            B(bl_index_r, bl_index_c) = ZERO

*           ..
*           .. Send cs-colvalues along column
*           ..
            IF (NNB1 .GT. 0) THEN
               CALL BSend(2, NNB1, COLUMNWISE, CSCOL(IND, t_I_c - 1), 
     $              LDC, CTX)
            ENDIF

            IF (NNB2 .GT. 1) THEN
               CALL BRecv(2, NNB2 - 1, COLUMNWISE, CSCOL(IND, I), 
     $              LDC, bl_owner_r, bl_owner_c, CTX)
            ENDIF

         ELSE IF ((MYCOL .EQ. bl_owner_c) .AND. (MYROW .NE. bl_owner_r) 
     $           .AND. (MYROW .NE. bl_bottom_r)) THEN

            t_I_c = I + NNB2
            IF ((MOD(I - 1, NB) .NE. 0) .AND. (NNB1 .GT. 0)) THEN
               CALL BRecv(2, NNB1, COLUMNWISE, CSCOL(IND, t_I_c - 1), 
     $              LDC, bl_bottom_r, bl_owner_c, CTX)
            ENDIF

            IF (NNB2 .GT. 1) THEN
               CALL BRecv(2, NNB2 - 1, COLUMNWISE, CSCOL(IND, I), 
     $              LDC, bl_owner_r, bl_owner_c, CTX)
            ENDIF
         ENDIF
         GOTO 999
      ENDIF

*     ..
*     .. Else we have at least a 2*2 processor grid
*     ..
      bl_owner_r = INDXG2P(I, NB, 0, 0, NPROW)
      bl_owner_c = INDXG2P(I, NB, 0, 0, NPCOL)
      bl_right_c = INDXG2P(I + NB1 - 1, NB, 0, 0, NPCOL)
      bl_bottom_r = INDXG2P(I + NB1 - 1, NB, 0, 0, NPROW)

      
*     ..
*     .. The processor who owns B(I,I)
*     ..
      IF ((MYROW .EQ. bl_owner_r) .AND. (MYCOL .EQ. bl_owner_c)) THEN
         I_r = I + NNB2 - 1 
         I_c = I_r
         bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
         bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

         IF ((MOD(I - 1, NB) .NE. 0) .AND. (NNB1 .GT. 0)) THEN
*           ..
*           .. Exchange elements for rowupdate
*           ..
            CALL Send(1, 1, DOWN, B(bl_index_r, bl_index_c), LDB,
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)
            TMPARR(1) = 0

            CALL DROT (1, B(bl_index_r,bl_index_c), LDB,
     $           TMPARR, 1,
     $           CSROW(IND, I_r), CSROW(IND + 1, I_r))

         ENDIF
*        ..
*        .. Preform rowrotation
*        ..
         I_r = I
         I_c = I
         bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
         bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

         CALL droRot(NNB2, B(bl_index_r, bl_index_c),
     $        LDB, CSROW(IND, I_r), LDC)

*        ..
*        .. Exchange element with processor to the right for colupdate
*        ..
         IF ((MOD(I - 1, NB) .NE. 0) .AND. (NNB1 .GT. 0)) THEN
            I_c = I + NNB2 - 1
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
            CALL Send(NNB2, 1, RIGHT, B(bl_index_r, bl_index_c),
     $           LDB, MYROW, MYCOL, NPCOL, NPROW, CTX)
            CALL Recv(NNB2, 1, RIGHT, TMPARR, 
     $           LDTMP, MYROW, MYCOL, NPCOL, NPROW, CTX)
*           ..
*           .. Preform colrotation
*           ..
            CALL Recv(2, 1, DOWN, CSCOL(IND, I_c), LDC,
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)

            CALL DROT(NNB2, TMPARR, 1,
     $           B(bl_index_r, bl_index_c), 1,
     $           CSCOL(IND, I_c), CSCOL(IND + 1, I_c))
         ENDIF

         I_c = I
         bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
*        ..
*        .. Diagonal column rot with elimination
*        ..
         CALL dcRot(NNB2, B(bl_index_r, bl_index_c),
     $        LDB, CSCOL(IND, I_c), LDC)

*        ..
*        .. Send cs-colvalues along column
*        .. 

         IF ((MOD(I - 1, NB) .NE. 0) .AND. (NNB1 .GT. 0)) THEN
            IF (NNB2 .GT. 0) THEN
               CALL BSend(2, NNB2, COLUMNWISE, CSCOL(IND, I), LDC, CTX)
            ENDIF
         ELSE IF (NNB2 .GT. 1) THEN
            CALL BSend(2, NNB2 - 1, COLUMNWISE, CSCOL(IND, I), LDC, CTX)
         ENDIF

      ELSE IF ((MYCOL.EQ.bl_owner_c).AND.(MYROW.NE.bl_bottom_r)) THEN
*        ..
*        .. Recv columnwise broadcasted cs-colvalues
*        ..   
         IF ((MOD(I - 1, NB) .NE. 0) .AND. (NNB1 .GT. 0)) THEN
            IF (NNB2 .GT. 0) THEN
               CALL BRecv(2, NNB2, COLUMNWISE, CSCOL(IND, I), LDC,
     $              bl_owner_r, bl_owner_c, CTX)
            ENDIF
         ELSE IF (NNB2 .GT. 1) THEN
            CALL BRecv(2, NNB2 - 1, COLUMNWISE, CSCOL(IND, I), LDC,
     $           bl_owner_r, bl_owner_c, CTX)
         ENDIF         
      ENDIF

*     ..
*     .. The processor to the right of the one who owns B(I,I)
*     ..                 
      IF ((MYROW .EQ. bl_owner_r) .AND. (MYCOL .EQ. bl_right_c).AND.
     $     (MYCOL .NE. bl_owner_c)) THEN

         I_c = I + NNB2
         I_r = I_c - 1
         bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)                  
         bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)                  

         NNB3 = MIN(N - I_c + 1, NB)

*        ..
*        .. Preform rowrotation
*        .. 
         CALL Send(1, NNB3, DOWN, B(bl_index_r, bl_index_c),
     $        LDB, MYROW, MYCOL, NPCOL, NPROW, CTX)
         CALL Recvtmp(1, NNB3, DOWN, TMPARR, LDTMP,
     $        MYROW, MYCOL, NPCOL, NPROW, CTX)


         CALL DROT(NNB3, B(bl_index_r, bl_index_c), LDB, 
     $        TMPARR, 1, CSROW(IND, I_r), CSROW(IND + 1, I_r))

         I_r = I
         bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
         
         CALL rRot(NNB2, NNB3, B(bl_index_r, bl_index_c),
     $        LDB, CSROW(IND, I_r), LDC, -1)
         
*        ..
*        .. Preform colrotation
*        ..

*        ..
*        .. Receive cs-colvalues, broadcasted earlier
*        ..
         CALL BRecv(2, NNB1, COLUMNWISE, CSCOL(IND, I_c - 1),
     $        LDC, bl_bottom_r, bl_right_c, CTX)

         IF (NNB1 .GT. 1) THEN
            CALL cRot(NNB1, NNB2, B(bl_index_r, bl_index_c),
     $           LDB, CSCOL(IND, I_c), LDC)
         ENDIF

         CALL Send(NNB2, 1, LEFT, B(bl_index_r, bl_index_c), LDB,
     $        MYROW, MYCOL, NPCOL, NPROW, CTX)

         CALL Recv(NNB2, 1, LEFT, TMPARR, LDTMP,
     $        MYROW, MYCOL, NPCOL, NPROW, CTX)
                  
         CALL DROT(NNB2, B( bl_index_r, bl_index_c), 1, TMPARR, 1,
     $        CSCOL(IND, I_c - 1), CSCOL(IND + 1, I_c - 1))

      ENDIF

*     ..
*     .. The processor below and to the right of the one who owns B(I,I)
*     ..
      IF ((MYROW .EQ. bl_bottom_r) .AND. (MYCOL .EQ. bl_right_c) .AND.
     $     (MYCOL .NE. bl_owner_c) .AND. (MYROW .NE. bl_owner_r)) THEN

         I_c = I + NNB2
         bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
         I_r = I_c
         bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
         NNB3 = MIN(N - I_c + 1, NB)

*        ..
*        .. Preform rowrotation
*        ..
         IF (NNB1 .GT. 1) THEN
            I_c = I_c + NNB1
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            CALL rRot(NNB1, NNB3-NNB1, B(bl_index_r, bl_index_c),
     $           LDB, CSROW(IND, I_r), LDC, -1)         

            I_c = I + NNB2
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
                  
            CALL droRot(NNB1, B(bl_index_r, bl_index_c), LDB,
     $           CSROW(IND, I_r), LDC)
         ENDIF

         CALL Send(1, NNB3, UP, B(bl_index_r, bl_index_c), LDB,
     $        MYROW, MYCOL, NPCOL, NPROW, CTX)

         CALL Recvtmp(1, NNB3, UP, TMPARR, LDTMP,
     $        MYROW, MYCOL, NPCOL, NPROW, CTX)

         CALL DROT(NNB3, TMPARR, 1,
     $        B(bl_index_r, bl_index_c), LDB,
     $        CSROW(IND, I_r - 1), CSROW(IND + 1, I_r - 1))

*        ..
*        .. Diagonal column rot with elimination
*        ..
         IF (NNB1 .GT. 1) THEN
            CALL dcRot(NNB1, B(bl_index_r, bl_index_c), LDB,
     $           CSCOL(IND, I_c), LDC)
         ENDIF

         CALL Send(1, 1, LEFT, B(bl_index_r, bl_index_c), LDB,
     $        MYROW, MYCOL, NPCOL, NPROW, CTX)

         CALL Recv(1, 1, LEFT, TMPARR, LDTMP, MYROW, MYCOL, 
     $        NPCOL, NPROW, CTX)

         TEMP = B(bl_index_r, bl_index_c)
         CALL DLARTG (TEMP, TMPARR, CSCOL(IND, I_c - 1),
     $        CSCOL(IND+1, I_c - 1), B(bl_index_r, bl_index_c))
*        ..
*        .. Broadcast cs-colvalues along column
*        ..
         CALL BSend(2, NNB1, COLUMNWISE, CSCOL(IND, I_c - 1),
     $        LDC, CTX)

      ELSE IF ((MYCOL.EQ.bl_right_c).AND.(MYROW.NE.bl_owner_r).AND.
     $        (MYCOL.NE.bl_owner_c)) THEN
*        ..
*        .. Receive broadcasted cs-colvalues
*        ..
         I_c = I + NNB2
         CALL BRecv(2, NNB1, COLUMNWISE, CSCOL(IND, I_c - 1),
     $        LDC, bl_bottom_r, bl_right_c, CTX)

      ENDIF


*     ..
*     .. The processor below the one who owns B(I,I)
*     ..
      IF ((MYROW .EQ. bl_bottom_r) .AND. (MYCOL.EQ.bl_owner_c) .AND.
     $     (MYROW .NE. bl_owner_r)) THEN

         I_r = I + NNB2
         I_c = I_r - 1
         bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)            
         bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

         CALL Recv(1, 1, UP, TMPARR, LDTMP, MYROW, MYCOL, 
     $        NPCOL, NPROW, CTX)

*        ..
*        .. Preform rowrotation
*        ..
         CALL DROT(1, TMPARR, 1,
     $        B(bl_index_r, bl_index_c), LDB,
     $        CSROW(IND, I_r - 1), CSROW(IND + 1, I_r - 1))
*        ..
*        .. Calculate elimination elements and preform colrotation
*        .. 
         CALL Send(1, 1, RIGHT, B(bl_index_r, bl_index_c), LDB,
     $        MYROW, MYCOL, NPCOL, NPROW, CTX)
         CALL Recv(1, 1, RIGHT, TMPARR, LDTMP,
     $        MYROW, MYCOL, NPCOL, NPROW, CTX)

         TEMP = TMPARR(1)
         CALL DLARTG (TEMP, B(bl_index_r, bl_index_c), CSCOL(IND, I_c),
     $        CSCOL(IND + 1, I_c), TMPARR)
         B(bl_index_r, bl_index_c) = ZERO

         CALL Send(2, 1, UP, CSCOL(IND, I_c), LDC, MYROW, MYCOL, 
     $           NPCOL, NPROW, CTX)
        
         IF (NNB2 .GT. 0) THEN
            CALL BRecv(2, NNB2, COLUMNWISE, CSCOL(IND, I), LDC,
     $           bl_owner_r, bl_owner_c, CTX)
         ENDIF
      ENDIF

*     ..
*     .. Update the remaining cols I:N and rows I:I+NB of B
*     ..
      I_c = I + NB
      IF (MOD(I_c - 1, NB) .NE. 0) THEN
         I_c = I_c + NB - MOD(I_c - 1, NB)
      ENDIF

      IF (I_c .GT. N) THEN
         GOTO 999
      ENDIF

      DO 10 t_I_c = I_c, N, NB
         NNB3 = MIN(N - t_I_c + 1, NB)
         bl_owner_c = INDXG2P(t_I_c, NB, 0, 0, NPCOL)
         IF ((MYROW .EQ. bl_owner_r) .AND. (MYCOL .EQ. bl_owner_c)) THEN                        
            bl_index_c = INDXG2L(t_I_c, NB, 0, 0, NPCOL)
            IF ((MOD(I - 1, NB) .NE. 0) .AND. (NNB1 .GT. 0)) THEN
               I_r = I + NNB2 - 1
               bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
               CALL Send(1, NNB3, DOWN, B(bl_index_r, bl_index_c),
     $              LDB, MYROW, MYCOL, NPCOL, NPROW, CTX)         
               CALL Recvtmp(1, NNB3, DOWN, TMPARR, LDTMP,
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)
               CALL DROT(NNB3,B(bl_index_r, bl_index_c), LDB, 
     $              TMPARR, 1,
     $              CSROW(IND, I_r), CSROW(IND + 1, I_r))
            ENDIF
            I_r = I
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            CALL rRot(NNB2, NNB3, B(bl_index_r, bl_index_c),
     $           LDB, CSROW(IND, I_r), LDC, -1)
         ENDIF

         IF ((MYROW .EQ. bl_bottom_r) .AND. (MYCOL .EQ. bl_owner_c) 
     $        .AND. (MYROW .NE. bl_owner_r)) THEN
            I_r = I + NNB2
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            bl_index_c = INDXG2L(t_I_c, NB, 0, 0, NPCOL)
            CALL rRot(NNB1, NNB3, B(bl_index_r, bl_index_c),
     $           LDB, CSROW(IND, I_r), LDC, -1)

            CALL Send(1, NNB3, UP, B(bl_index_r, bl_index_c), LDB,
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)

            CALL Recvtmp(1, NNB3, UP, TMPARR, LDTMP,
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)

            CALL DROT(NNB3, TMPARR, 1,
     $           B(bl_index_r, bl_index_c), LDB,
     $           CSROW(IND, I_r - 1), CSROW(IND + 1, I_r - 1))

         ENDIF

 10   CONTINUE

     
 999  RETURN
      END

      SUBROUTINE pUPDB(B,LDB,CSROW,CSCOL,LDC,I,JC,MAXI,NB,NB1,
     $     NPROW,NPCOL,IND,TMPARR,LDTMP,MYROW,MYCOL,CTX)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER           I, JC, MAXI, NB, NB1, NPROW, NPCOL, LDB, LDC 
      INTEGER           IND, LDTMP,MYROW, MYCOL, CTX

*     .. Array Arguments ..
      DOUBLE PRECISION   B(LDB,*)
      DOUBLE PRECISION   CSCOL(LDC,*), CSROW(LDC,*)
      DOUBLE PRECISION   TMPARR(LDTMP)

*     .. Parameters ..
      INTEGER            DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )

*     .. Local scalars ..
      INTEGER            bl_right_c, bl_r, bl_c, NNB1,NNB2
      INTEGER            I_c, I_r, bl_index_r, bl_index_c
      INTEGER            t_I_c, t_I_r, t_bl_index_r, t_bl_index_c

*     .. External Functions ..
      EXTERNAL           INDXG2L, INDXG2P
      INTEGER            INDXG2L, INDXG2P

*     .. Start execution ..

      NNB1 = MOD(I - 1, NB) - (NB - NB1)
      NNB2 = NB - MOD(I - 1, NB)

      IF (NNB1 .LT. 0) THEN
         NNB2 = NNB2 + NNB1
      ENDIF

      bl_c = INDXG2P(I, NB, 0, 0, NPCOL)
      bl_right_c = INDXG2P(I + NB1 - 1, NB, 0, 0, NPCOL)
      
      IF (MOD(I - 1, NB) .EQ. 0 ) THEN
         bl_r = INDXG2P(I, NB, 0, 0, NPROW)
         IF ((MYROW .EQ. bl_r) .AND. (MYCOL .EQ. bl_c)) THEN
*           ..
*           .. Rowrotate toprow with row from block above
*           ..
            IF (((I - 1) .GT. JC) .AND. 
     $           (MOD(I - 1, NB) .NE. MOD(JC,NB)) 
     $           .AND. (NPCOL .EQ. 1)) THEN
               bl_index_r = INDXG2L(I, NB, 0, 0, NPROW)
               bl_index_c = INDXG2L(I, NB, 0, 0, NPCOL)
               IF (NPROW .GT. 1) THEN
                  CALL Send(1, NNB2, UP, B(bl_index_r, bl_index_c), LDB,
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)

                  CALL Recvtmp(1, NNB2, UP, TMPARR, LDTMP,
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)
               
                  CALL DROT(NNB2, TMPARR, 1, B(bl_index_r, bl_index_c),
     $                 LDB, CSROW(IND, I - 1),
     $                 CSROW(IND + 1, I - 1))    
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      bl_r = INDXG2P(MAXI, NB, 0, 0, NPROW)
*     ..
*     .. If NPROW = 1 then NPCOL > 1 and vice versa
*     ..

      IF (NPROW .EQ. 1) THEN
*        ..
*        .. The processor who owns B(MAXI,I)
*        ..
         IF (MYCOL .EQ. bl_c) THEN
*           ..
*           .. Preform colrotation
*           ..
            I_r = MAXI - MOD(MAXI - 1, NB)
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            I_c = I + NNB2 - 1
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            IF (NNB1 .GT. 0) THEN
               CALL Send(MOD(MAXI - 1, NB) + 1, 1, RIGHT, 
     $              B(bl_index_r, bl_index_c), LDB, MYROW, MYCOL, 
     $              NPCOL, NPROW, CTX)
        
               CALL Recv(MOD(MAXI - 1, NB) + 1, 1, RIGHT, TMPARR, LDTMP, 
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)

               CALL DROT(MOD(MAXI - 1, NB) + 1, 
     $              TMPARR, 1,
     $              B(bl_index_r, bl_index_c), 1, 
     $              CSCOL(IND, I_c), CSCOL(IND + 1, I_c))
            ENDIF

            I_c = I
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            IF (NNB2 .GT. 1) THEN
               CALL cRot(NNB2, MOD(MAXI - 1, NB) + 1,
     $              B(bl_index_r, bl_index_c), LDB, CSCOL(IND, I_c),
     $              LDC) 
            ENDIF

*           ..
*           .. Preform rowrotation
*           ..
            IF (MOD(I - 1, NB) .EQ. 0) THEN
               I_r = MAXI
               bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
               IF ((I_r .GT. JC) .AND. (MOD(I_r, NB) .NE. MOD(JC,NB))) 
     $              THEN
                  t_I_r = I
                  t_bl_index_r = INDXG2L(t_I_r, NB, 0, 0, NPROW)
                  CALL DROT(NNB2, B(bl_index_r, bl_index_c), LDB, 
     $                 B(t_bl_index_r, bl_index_c), LDB, 
     $                 CSROW(IND, I - 1), CSROW(IND + 1, I - 1))
               ENDIF
               I_r = MAXI - MOD(MAXI - 1, NB)
               bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            ENDIF
            
            IF (MOD(MAXI - 1, NB) .GT. 0) THEN
               IF (I_r .LE. JC) THEN
                  
                  bl_index_r = INDXG2L(JC + 1, NB, 0, 0, NPROW)
                  IF ((MAXI - JC) .GT. 1) THEN
                     CALL rRot(MAXI - JC, NNB2,
     $                    B(bl_index_r, bl_index_c), LDB,
     $                    CSROW(IND, JC + 1), LDC, -1)
                  ENDIF

               ELSE
                  CALL rRot(MOD(MAXI - 1, NB) + 1, NNB2,
     $                 B(bl_index_r, bl_index_c), LDB,
     $                 CSROW(IND, I_r), LDC,
     $                 MOD(JC, NB))

               ENDIF
            ENDIF
         ENDIF
*        ..
*        .. The processor to the right of the one who owns B(MAXI,I)
*        ..
         IF ((MYCOL .EQ. bl_right_c) .AND. (MYCOL .NE. bl_c)) THEN
*           ..
*           .. Preform colrotation
*           ..
            I_r = MAXI-MOD(MAXI-1,NB)
            I_c = I + NNB2
            bl_index_c = INDXG2L(I_c,NB,0,0,NPCOL)
            bl_index_r = INDXG2L(I_r,NB,0,0,NPROW)

            CALL cRot(NNB1, MOD(MAXI - 1, NB) + 1,
     $           B(bl_index_r, bl_index_c), LDB, CSCOL(IND, I_c),
     $           LDC)

            CALL Send(MOD(MAXI - 1, NB) + 1, 1, LEFT, 
     $           B(bl_index_r, bl_index_c), LDB, 
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)

            CALL Recv(MOD(MAXI - 1, NB) + 1, 1, LEFT, TMPARR, LDTMP, 
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)

            CALL DROT(MOD(MAXI - 1, NB) + 1,
     $           B(bl_index_r, bl_index_c), 1, TMPARR, 1,
     $           CSCOL(IND, I_c - 1), CSCOL(IND + 1, I_c - 1) )

*           ..
*           .. Preform rowrotation
*           ..
            IF (MOD(MAXI - 1, NB) .GT. 0) THEN
               IF (I_r .LE. JC) THEN
                  
                  bl_index_r = INDXG2L(JC + 1, NB, 0, 0, NPROW)
                  IF ((MAXI - JC) .GT. 1) THEN
                     CALL rRot(MAXI - JC, NNB1,
     $                    B(bl_index_r, bl_index_c), LDB,
     $                    CSROW(IND, JC + 1), LDC, -1)
                  ENDIF

               ELSE

                  CALL rRot(MOD(MAXI - 1, NB) + 1, NNB1,
     $                 B(bl_index_r, bl_index_c), LDB,
     $                 CSROW(IND, I_r), LDC,
     $                 MOD(JC, NB))

               ENDIF
            ENDIF
         ENDIF

         DO 2 I_r = MAXI - MOD(MAXI - 1, NB) - 1, NB, -NB          
            IF (MYCOL .EQ. bl_c) THEN
*              ..
*              .. Preform colrotation
*              ..
               bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)

               I_c = I + NNB2 - 1
               bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
 
               IF (NNB1 .GT. 0) THEN
                  CALL Send(NB, 1, RIGHT, B(bl_index_r, bl_index_c),
     $                 LDB, MYROW, MYCOL, NPCOL, NPROW, CTX)

                  CALL Recv(NB, 1, RIGHT, TMPARR, LDTMP,
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)

                  CALL DROT(NB, TMPARR, 1,
     $                 B(bl_index_r, bl_index_c), 1, 
     $                 CSCOL(IND, I_c), CSCOL(IND + 1, I_c))
               ENDIF
            
               I_c = I
               bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)           

               IF (NNB2 .GT. 1) THEN
                  CALL cRot(NNB2, NB,
     $                 B(bl_index_r, bl_index_c), LDB, CSCOL(IND, I_c),
     $                 LDC)
               ENDIF

               bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)  

*              ..
*              .. Preform rowrotation
*              ..
               IF ((I_r .GT. JC) .AND. (MOD(I_r,NB) .NE. MOD(JC,NB))) 
     $              THEN
                  t_I_r = I_r + 1
                  t_bl_index_r = INDXG2L(t_I_r, NB, 0, 0, NPROW)

                  CALL DROT(NNB2, B(bl_index_r, bl_index_c), LDB, 
     $                 B(t_bl_index_r, bl_index_c), LDB,
     $                 CSROW(IND, I_r), CSROW(IND + 1, I_r))
               ENDIF

               bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)

               IF ((I_r - NB) .LE. JC) THEN
                  
                  bl_index_r = INDXG2L(JC + 1, NB, 0, 0, NPROW)
                  IF ((I_r - JC) .GT. 1) THEN
                     CALL rRot(I_r - JC, NNB2,
     $                    B(bl_index_r, bl_index_c), LDB,
     $                    CSROW(IND, JC + 1), LDC, -1)
                  ENDIF

               ELSE 

                  CALL rRot(NB, NNB2,
     $                 B(bl_index_r, bl_index_c), LDB,
     $                 CSROW(IND, I_r - NB + 1), LDC,
     $                 MOD(JC, NB))

               ENDIF
            ENDIF
            
            IF ((MYCOL .EQ. bl_right_c) .AND. (MYCOL .NE. bl_c)) THEN
*              ..
*              .. Preform colrotation
*              ..
               I_c = I + NNB2
               bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
               bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)  
               CALL cRot(NNB1, NB,
     $              B(bl_index_r, bl_index_c), LDB, CSCOL(IND, I_c),
     $              LDC)

               CALL Send(NB, 1, LEFT, B(bl_index_r, bl_index_c), 
     $              LDB, MYROW, MYCOL,
     $              NPCOL, NPROW, CTX)

               CALL Recv(NB, 1, LEFT, TMPARR, LDTMP, MYROW, MYCOL,
     $              NPCOL, NPROW, CTX)

               CALL DROT(NB, B( bl_index_r, bl_index_c),
     $              1, TMPARR, 1, 
     $              CSCOL(IND, I_c - 1), CSCOL(IND + 1, I_c - 1) )

*              ..
*              .. Preform rowrotation
*              ..                     
               bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)  

               IF ((I_r .GT. JC) .AND. (MOD(I_r,NB) .NE. MOD(JC,NB))) 
     $              THEN
                  t_I_r = I_r + 1
                  t_bl_index_r = INDXG2L(t_I_r, NB, 0, 0, NPROW)

                  CALL DROT(NNB1, B(bl_index_r, bl_index_c), LDB, 
     $                 B(t_bl_index_r, bl_index_c), LDB,
     $                 CSROW(IND, I_r), CSROW(IND + 1, I_r))
               ENDIF 
           
               bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)   

               IF ((I_r - NB) .LE. JC) THEN
                  
                  bl_index_r = INDXG2L(JC+1, NB, 0, 0, NPROW)
                  IF ((I_r - JC) .GT. 1) THEN
                     CALL rRot(I_r - JC, NNB1,
     $                    B(bl_index_r, bl_index_c), LDB,
     $                    CSROW(IND, JC + 1), LDC, -1)
                  ENDIF
                  
               ELSE
                  
                  CALL rRot(NB, NNB1,
     $                 B(bl_index_r, bl_index_c), LDB,
     $                 CSROW(IND, I_r - NB + 1), LDC,
     $                 MOD(JC, NB))
                  
               ENDIF            
            ENDIF
 2       CONTINUE
         GOTO 999
      ENDIF

      IF (NPCOL .EQ. 1) THEN
*        ..
*        .. The processor who owns B(MAXI,I)
*        ..
         IF (MYROW .EQ. bl_r) THEN
            I_r = MAXI - MOD(MAXI - 1, NB)
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            I_c = I + NNB2 - 1
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
            
*           ..
*           .. Preform colrotation
*           ..
            IF (NNB1 .GT. 0) THEN
               t_I_c = I + NNB2
               t_bl_index_c = INDXG2L(t_I_c, NB, 0, 0, NPCOL)

               CALL cRot(NNB1, MOD(MAXI - 1, NB) + 1,
     $              B(bl_index_r, t_bl_index_c), LDB, CSCOL(IND, t_I_c),
     $              LDC)

               CALL DROT(MOD(MAXI - 1, NB) + 1,
     $              B(bl_index_r, t_bl_index_c), 1, 
     $              B(bl_index_r, bl_index_c), 1,
     $              CSCOL(IND, t_I_c - 1), CSCOL(IND + 1, t_I_c - 1) )
            ENDIF

            I_c = I
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            IF (NNB2 .GT. 1) THEN
               CALL cRot(NNB2, MOD(MAXI - 1, NB) + 1,
     $              B(bl_index_r, bl_index_c), LDB, CSCOL(IND, I_c),
     $              LDC) 
            ENDIF


*           ..
*           .. Preform rowrotation
*           ..
            IF (MOD(I - 1, NB) .EQ. 0) THEN
               I_r = MAXI
               bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
               IF ((I_r .GT. JC) .AND. (MOD(I_r, NB) .NE. MOD(JC,NB))) 
     $              THEN
                  CALL Send(1, NNB2, DOWN, B(bl_index_r, bl_index_c),
     $                 LDB, MYROW, MYCOL, NPCOL, NPROW, CTX) 

                  CALL Recvtmp(1, NNB2, DOWN, TMPARR, LDTMP,
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)

                  CALL DROT(NNB2, B(bl_index_r, bl_index_c), LDB, 
     $                 TMPARR, 1,
     $                 CSROW(IND, I_r), CSROW(IND + 1, I_r))     
               ENDIF
               I_r = MAXI - MOD(MAXI - 1, NB)
               bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            ENDIF
            
            IF (MOD(MAXI - 1, NB) .GT. 0) THEN
               IF (I_r .LE. JC) THEN

                  bl_index_r = INDXG2L(JC + 1, NB, 0, 0, NPROW)
                  IF ((MAXI - JC) .GT. 1) THEN
                     CALL rRot(MAXI - JC, NNB2,
     $                    B(bl_index_r, bl_index_c), LDB,
     $                    CSROW(IND, JC + 1), LDC, -1)
                     IF (NNB1 .GT. 0) THEN
                        CALL rRot(MAXI - JC, NNB1,
     $                       B(bl_index_r, t_bl_index_c), LDB,
     $                       CSROW(IND, JC + 1), LDC, -1)
                     ENDIF
                  ENDIF
                  bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
                  
               ELSE
                  CALL rRot(MOD(MAXI - 1, NB) + 1, NNB2,
     $                 B(bl_index_r, bl_index_c), LDB,
     $                 CSROW(IND, I_r), LDC,
     $                 MOD(JC, NB))
                  IF (NNB1 .GT. 0) THEN
                     CALL rRot(MOD(MAXI - 1, NB) + 1, NNB1,
     $                    B(bl_index_r, t_bl_index_c), LDB,
     $                    CSROW(IND, I_r), LDC,
     $                    MOD(JC, NB))
                  ENDIF

               ENDIF
            ENDIF

*           ..
*           .. Update rowvalues and send to next block above
*           ..
            IF (((I_r - 1) .GT. JC) .AND.
     $           (MOD(I_r - 1, NB) .NE. MOD(JC, NB))) THEN
               IF (NNB1 .GT. 0) THEN
                  CALL Send(1, NNB1, UP, B(bl_index_r, t_bl_index_c),
     $                 LDB, MYROW, MYCOL, NPCOL, NPROW, CTX)

                  CALL Recvtmp(1, NNB1, UP, TMPARR, LDTMP,
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)
                     
                  CALL DROT(NNB1, TMPARR, 1, 
     $                 B(bl_index_r, t_bl_index_c), LDB, 
     $                 CSROW(IND, I_r - 1), CSROW(IND + 1, I_r - 1))

               ENDIF
               CALL Send(1, NNB2, UP, B(bl_index_r, bl_index_c),
     $              LDB, MYROW, MYCOL, NPCOL, NPROW, CTX)

               CALL Recvtmp(1, NNB2, UP, TMPARR, LDTMP,
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)
               
               CALL DROT(NNB2, TMPARR, 1, B(bl_index_r, bl_index_c),
     $              LDB, CSROW(IND, I_r - 1), CSROW(IND + 1,I_r - 1))

            ENDIF

         ENDIF
         DO 3 I_r = MAXI - MOD(MAXI - 1, NB) - 1, NB, -NB
            bl_r = INDXG2P(I_r, NB, 0, 0, NPROW)         
            
            IF (MYROW .EQ. bl_r) THEN
               
               bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)

               t_I_c = I + NNB2
               t_bl_index_c = INDXG2L(t_I_c, NB, 0, 0, NPCOL)

               I_c = I + NNB2 - 1
               bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
*              ..
*              .. Preform colrotation
*              ..
               IF (NNB1 .GT. 0) THEN
                  CALL cRot(NNB1, NB,
     $                 B(bl_index_r, t_bl_index_c), LDB, 
     $                 CSCOL(IND, t_I_c),
     $                 LDC)
                  
                  CALL DROT(NB, B(bl_index_r, t_bl_index_c), 1, 
     $                 B(bl_index_r, bl_index_c), 1, 
     $                 CSCOL(IND, t_I_c - 1), 
     $                 CSCOL(IND + 1, t_I_c - 1) )
                  
               ENDIF
             
               I_c = I
               bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)           
               
               IF (NNB2 .GT. 1) THEN
                  CALL cRot(NNB2, NB,
     $                 B(bl_index_r, bl_index_c), LDB, 
     $                 CSCOL(IND, I_c), LDC)
               ENDIF
               
               bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)  
               
*                 ..
*                 .. Preform rowrotation
*                 ..
               IF ((I_r .GT. JC) .AND. (MOD(I_r, NB) .NE. MOD(JC, NB))) 
     $              THEN
                  IF (NNB1 .GT. 0) THEN
                     CALL Send(1, NNB1, DOWN, 
     $                    B(bl_index_r, t_bl_index_c), LDB,
     $                    MYROW, MYCOL, NPCOL, NPROW, CTX)
                     
                     CALL Recvtmp(1, NNB1, DOWN, 
     $                    TMPARR, LDTMP,
     $                    MYROW, MYCOL, NPCOL, NPROW, CTX)

                     CALL DROT(NNB1, B(bl_index_r, t_bl_index_c), LDB,
     $                    TMPARR, 1,
     $                    CSROW(IND, I_r), 
     $                    CSROW(IND + 1, I_r))
                  ENDIF
                  
                  CALL Send(1, NNB2, DOWN, B(bl_index_r, bl_index_c),
     $                 LDB, MYROW, MYCOL, NPCOL, NPROW, CTX)
                  
                  CALL Recvtmp(1, NNB2, DOWN, TMPARR, LDTMP,
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)

                  CALL DROT(NNB2, B(bl_index_r, bl_index_c), LDB,
     $                 TMPARR, 1,
     $                 CSROW(IND, I_r), 
     $                 CSROW(IND + 1, I_r))
               ENDIF

               bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)
               
               IF ((I_r - NB) .LE. JC) THEN
                  
                  bl_index_r = INDXG2L(JC + 1, NB, 0, 0, NPROW)
                  IF ((I_r - JC) .GT. 1) THEN
                     CALL rRot(I_r - JC, NNB2,
     $                    B(bl_index_r, bl_index_c), LDB,
     $                    CSROW(IND, JC + 1), LDC, -1)
                     IF (NNB1 .GT. 0) THEN
                        CALL rRot(I_r - JC, NNB1,
     $                       B(bl_index_r, t_bl_index_c), LDB,
     $                       CSROW(IND, JC + 1), LDC, -1)
                     ENDIF
                  ENDIF
                  bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)
                  
               ELSE 

                  CALL rRot(NB, NNB2,
     $                 B(bl_index_r, bl_index_c), LDB,
     $                 CSROW(IND, I_r - NB + 1), LDC,
     $                 MOD(JC, NB))
                  
                  IF (NNB1 .GT. 0) THEN
                     CALL rRot(NB, NNB1,
     $                    B(bl_index_r, t_bl_index_c), LDB,
     $                    CSROW(IND, I_r - NB + 1), LDC,
     $                    MOD(JC, NB))
                  ENDIF
                  
               ENDIF

               
*              ..
*              .. Send updated rowvalues to next block above
*              ..
               IF (((I_r - NB) .GT. JC) .AND.
     $              (MOD(I_r, NB) .NE. MOD(JC, NB))) THEN
                  IF (NNB1 .GT. 0) THEN
                     CALL Send(1, NNB1, UP, B(bl_index_r, t_bl_index_c),
     $                    LDB, MYROW, MYCOL, NPCOL, NPROW, CTX)

                     CALL Recvtmp(1, NNB1, UP, TMPARR, LDTMP,
     $                    MYROW, MYCOL, NPCOL, NPROW, CTX)
                     
                     CALL DROT(NNB1, TMPARR, 1, 
     $                    B(bl_index_r, t_bl_index_c),
     $                    LDB, CSROW(IND, I_r - NB), 
     $                    CSROW(IND + 1, I_r - NB))
                  ENDIF

                  CALL Send(1, NNB2, UP, B(bl_index_r, bl_index_c), LDB,
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)

                  CALL Recvtmp(1, NNB2, UP, TMPARR, LDTMP,
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)
                  
                  CALL DROT(NNB2, TMPARR, 1, 
     $                 B(bl_index_r, bl_index_c), LDB,
     $                 CSROW(IND, I_r - NB), 
     $                 CSROW(IND + 1, I_r - NB))
               ENDIF
            ENDIF
 3       CONTINUE

         GOTO 999
      ENDIF



      bl_r = INDXG2P(MAXI, NB, 0, 0, NPROW)

*     ..
*     .. Else we have at least a 2*2 processor grid
*     ..
      IF ((MYROW .EQ. bl_r) .AND. (MYCOL .EQ. bl_c)) THEN
*        ..
*        .. Preform colrotation
*        ..
         I_r = MAXI - MOD(MAXI - 1, NB)
         bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
         I_c = I + NNB2 - 1
         bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

         IF (NNB1 .GT. 0) THEN
            CALL Send(MOD(MAXI - 1, NB) + 1, 1, RIGHT, 
     $           B(bl_index_r, bl_index_c), LDB, MYROW, MYCOL, 
     $           NPCOL, NPROW, CTX)
        
            CALL Recv(MOD(MAXI - 1, NB) + 1, 1, RIGHT, 
     $           TMPARR, LDTMP, MYROW, MYCOL, 
     $           NPCOL, NPROW, CTX)

            CALL DROT(MOD(MAXI - 1, NB) + 1, 
     $           TMPARR, 1,
     $           B(bl_index_r, bl_index_c), 1, 
     $           CSCOL(IND, I_c), CSCOL(IND + 1, I_c))

         ENDIF

         I_c = I
         bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

         IF (NNB2 .GT. 1) THEN
            CALL cRot(NNB2, MOD(MAXI - 1, NB) + 1,
     $           B(bl_index_r, bl_index_c), LDB, CSCOL(IND, I_c),
     $           LDC) 
         ENDIF

      ENDIF
      
      IF ((MYROW.EQ.bl_r).AND.(MYCOL.EQ.bl_right_c).AND.(MYCOL.NE.bl_c)) 
     $     THEN
*        ..
*        .. Preform colrotation
*        ..
         I_r = MAXI-MOD(MAXI-1,NB)
         I_c = I + NNB2
         bl_index_c = INDXG2L(I_c,NB,0,0,NPCOL)
         bl_index_r = INDXG2L(I_r,NB,0,0,NPROW)
         CALL cRot(NNB1, MOD(MAXI - 1, NB) + 1,
     $        B(bl_index_r, bl_index_c), LDB, CSCOL(IND, I_c),
     $        LDC)

         CALL Send(MOD(MAXI - 1, NB) + 1, 1, LEFT, 
     $        B(bl_index_r, bl_index_c), LDB, 
     $        MYROW, MYCOL, NPCOL, NPROW, CTX)

         CALL Recv(MOD(MAXI - 1, NB) + 1, 1, LEFT, TMPARR, LDTMP, 
     $        MYROW, MYCOL, NPCOL, NPROW, CTX)

         CALL DROT(MOD(MAXI - 1, NB) + 1,
     $        B(bl_index_r, bl_index_c), 1, TMPARR, 1,
     $        CSCOL(IND, I_c - 1), CSCOL(IND + 1, I_c - 1))

      ENDIF


      DO 10 I_r = MAXI - MOD(MAXI - 1, NB) - 1, NB, -NB
         bl_r = INDXG2P(I_r, NB, 0, 0, NPROW)
         
         IF ((MYROW .EQ. bl_r) .AND. (MYCOL .EQ. bl_c)) THEN
*           ..
*           .. Preform colrotation
*           ..
            bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)

            I_c = I + NNB2 - 1
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            IF (NNB1 .GT. 0) THEN
               CALL Send(NB, 1, RIGHT, B(bl_index_r, bl_index_c), LDB,
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)

               CALL Recv(NB, 1, RIGHT, TMPARR, LDTMP,
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)

               CALL DROT(NB, 
     $              TMPARR, 1,
     $              B(bl_index_r, bl_index_c), 1, 
     $              CSCOL(IND, I_c), CSCOL(IND + 1, I_c))               
            ENDIF
            
            I_c = I
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)           

            IF (NNB2 .GT. 1) THEN
               CALL cRot(NNB2, NB,
     $              B(bl_index_r, bl_index_c), LDB, CSCOL(IND, I_c),
     $              LDC)
            ENDIF

         ENDIF
                  
         IF ((MYROW.EQ.bl_r).AND.(MYCOL.EQ.bl_right_c).AND.
     $        (MYCOL.NE.bl_c)) THEN
*           ..
*           .. Preform colrotation
*           .. 
            I_c = I + NNB2
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
            bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)  

            CALL cRot(NNB1, NB,
     $           B(bl_index_r, bl_index_c), LDB, CSCOL(IND, I_c),
     $           LDC)

            CALL Send(NB, 1, LEFT, B(bl_index_r, bl_index_c), 
     $           LDB, MYROW, MYCOL,
     $           NPCOL, NPROW, CTX)

            CALL Recv(NB, 1, LEFT, TMPARR, LDTMP, MYROW, MYCOL,
     $           NPCOL, NPROW, CTX)

            CALL DROT(NB, B( bl_index_r, bl_index_c),
     $           1, TMPARR, 1, 
     $           CSCOL(IND, I_c - 1), CSCOL(IND + 1, I_c - 1) )

         ENDIF
10    CONTINUE

 999  RETURN
      END


      SUBROUTINE pUPDA(A,LDA,CSROW,CSCOL,LDC,I,JC,MAXI,NB,NB1,NPROW,
     $     NPCOL,IND,TMPARR,LDTMP,MYROW,MYCOL,CTX)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            I, MAXI, NB, NB1, NPROW, NPCOL, LDA, LDC, IND 
      INTEGER            LDTMP, MYROW, MYCOL, CTX, JC

*     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*)
      DOUBLE PRECISION   CSCOL(LDC,*), CSROW(LDC,*)
      DOUBLE PRECISION   TMPARR(LDTMP)

*     .. Parameters ..
      INTEGER            DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )

*     .. Local scalars ..
      INTEGER            bl_right_c, bl_r, bl_c, NNB1,NNB2
      INTEGER            I_c, I_r, bl_index_r, bl_index_c
      INTEGER            t_I_c, t_I_r, t_bl_index_r, t_bl_index_c

*     .. External Functions ..
      EXTERNAL           INDXG2L, INDXG2P
      INTEGER            INDXG2L,INDXG2P


*     ..
*     .. Start execution
*     ..
      NNB1 = MOD(I - 1, NB) - (NB - NB1)
      NNB2 = NB - MOD(I - 1, NB)

      IF (NNB1 .LT. 0) THEN
         NNB2 = NNB2 + NNB1
      ENDIF

*     ..
*     .. Apply rotations to block above diagonal for A
*     ..
      bl_r = INDXG2P(MAXI, NB, 0, 0, NPROW)
      bl_c = INDXG2P(I, NB, 0, 0, NPCOL)
      bl_right_c = INDXG2P(I + NB1 - 1, NB, 0, 0, NPCOL)

*     ..
*     .. If NPROW = 1 then NPCOL > 1 and vice versa
*     ..
      IF (NPROW .EQ. 1) THEN
*        ..
*        .. The processor who owns A(MAXI,I)
*        ..
         IF (MYCOL .EQ. bl_c) THEN

*           ..
*           .. Preform colrotation
*           ..
            I_r = MAXI - MOD(MAXI - 1, NB)
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)

            I_c = I + NNB2 - 1 
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            IF (NNB1 .GT. 0) THEN
               CALL Send(MOD(MAXI - 1,NB) + 1, 1, RIGHT, 
     $              A(bl_index_r, bl_index_c),
     $              LDA, MYROW, MYCOL, NPCOL, NPROW, CTX)

               CALL Recv(MOD(MAXI - 1, NB) + 1, 1, RIGHT, 
     $              TMPARR, LDTMP,
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)

               CALL DROT(MOD(MAXI - 1, NB) + 1, 
     $              TMPARR, 1,
     $              A(bl_index_r, bl_index_c), 1, 
     $              CSCOL(IND, I_c), CSCOL(IND + 1, I_c))
            ENDIF

            I_c = I
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
            CALL cRot(NNB2, MOD(MAXI - 1, NB) + 1,
     $           A(bl_index_r, bl_index_c), LDA, CSCOL(IND, I_c),
     $           LDC) 

*           ..
*           .. Preform rowrotation
*           ..
            IF (MOD(MAXI - 1, NB) .GT. 0) THEN
               IF (I_r .LE. JC) THEN
                  
                  bl_index_r = INDXG2L(JC + 1, NB, 0, 0, NPROW)
                  IF ((MAXI - JC) .GT. 1) THEN
                     CALL rRot(MAXI - JC, NNB2,
     $                    A(bl_index_r, bl_index_c), LDA,
     $                    CSROW(IND, JC + 1), LDC, -1)
                  ENDIF
                  bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
                  
               ELSE

                  CALL rRot(MOD(MAXI - 1, NB) + 1, NNB2,
     $                 A(bl_index_r, bl_index_c), LDA,
     $                 CSROW(IND, I_r), LDC,
     $                 MOD(JC, NB))

               ENDIF
            ENDIF
         ENDIF

*        ..
*        .. The processor to the right of the one who owns A(MAXI,I)
*        ..
         IF ((MYCOL .EQ. bl_right_c) .AND. (MYCOL .NE. bl_c)) THEN
*           ..
*           .. Preform colrotation
*           ..
            I_r = MAXI - MOD(MAXI - 1, NB)
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            
            I_c = I + NNB2
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            CALL cRot(NNB1, MOD(MAXI - 1, NB) + 1,
     $           A(bl_index_r, bl_index_c), LDA, CSCOL(IND, I_c),
     $           LDC)

            CALL Send(MOD(MAXI - 1, NB) + 1, 1, LEFT, 
     $           A(bl_index_r, bl_index_c), LDA,
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)
            
            CALL Recv(MOD(MAXI - 1, NB) + 1, 1, LEFT, TMPARR, LDTMP,
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)
            
            CALL DROT(MOD(MAXI - 1, NB) + 1, A( bl_index_r, bl_index_c),
     $           1, TMPARR, 1,
     $           CSCOL(IND, I_c - 1), CSCOL(IND + 1, I_c - 1))
            
*           ..
*           .. Preform rowrotation
*           ..
            IF (MOD(MAXI - 1, NB) .GT. 0) THEN
               IF (I_r .LE. JC) THEN
                  
                  bl_index_r = INDXG2L(JC + 1, NB, 0, 0, NPROW)
                  IF ((MAXI - JC) .GT. 1) THEN
                     CALL rRot(MAXI - JC, NNB1,
     $                    A(bl_index_r, bl_index_c), LDA,
     $                    CSROW(IND, JC + 1), LDC, -1)
                  ENDIF
                  bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
                  
               ELSE
                  CALL rRot(MOD(MAXI - 1, NB) + 1, NNB1,
     $                 A(bl_index_r, bl_index_c), LDA,
     $                 CSROW(IND, I_r), LDC,
     $                 MOD(JC, NB))
               ENDIF
            ENDIF            
         ENDIF

         DO 2 I_r = MAXI - MOD(MAXI - 1, NB) - 1, NB, -NB
            IF (MYCOL .EQ. bl_c) THEN
*              ..
*              .. Preform colrotation
*              ..
               bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)  
               I_c = I + NNB2 - 1                     
               bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

               IF (NNB1 .GT. 0) THEN
                  CALL Send(NB, 1, RIGHT, A(bl_index_r, bl_index_c),
     $                 LDA, MYROW, MYCOL, NPCOL, NPROW, CTX)

                  CALL Recv(NB, 1, RIGHT, TMPARR, LDTMP, 
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)
                  
                  CALL DROT(NB, TMPARR, 1,
     $                 A(bl_index_r, bl_index_c), 1, 
     $                 CSCOL(IND, I_c), CSCOL(IND + 1, I_c))
               ENDIF


               I_c = I
               bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
               CALL cRot(NNB2, NB, A(bl_index_r, bl_index_c), LDA, 
     $              CSCOL(IND,I_c), LDC)

*              ..
*              .. Preform rowrotation
*              ..
               bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)  
               IF ((I_r .GT. JC) .AND. (MOD(I_r,NB) .NE. MOD(JC,NB))) 
     $              THEN
                  t_I_r = I_r + 1
                  t_bl_index_r = INDXG2L(t_I_r, NB, 0, 0, NPROW)

                  CALL DROT(NNB2, A(bl_index_r, bl_index_c), LDA, 
     $                 A(t_bl_index_r, bl_index_c), LDA,
     $                 CSROW(IND, I_r), CSROW(IND + 1, I_r))

               ENDIF
               bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)
               
               IF ((I_r - NB) .LE. JC) THEN               

                  bl_index_r = INDXG2L(JC + 1, NB, 0, 0, NPROW)
                  IF ((I_r - JC) .GT. 1) THEN
                     CALL rRot(I_r - JC, NNB2,
     $                    A(bl_index_r, bl_index_c), LDA,
     $                    CSROW(IND, JC + 1), LDC, -1)
                  ENDIF
                  bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)

               ELSE

                  CALL rRot(NB, NNB2,
     $                 A(bl_index_r, bl_index_c), LDA,
     $                 CSROW(IND, I_r - NB + 1), LDC,
     $                 MOD(JC, NB))
                  
               ENDIF
            ENDIF

            IF ((MYCOL .EQ. bl_right_c) .AND. (MYCOL .NE. bl_c)) THEN
*              ..
*              .. Preform colrotation
*              ..
               I_c = I + NNB2
               bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
               bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)

               CALL cRot(NNB1, NB, A(bl_index_r, bl_index_c), LDA, 
     $              CSCOL(IND, I_c), LDC)

               CALL Send(NB, 1, LEFT, A(bl_index_r, bl_index_c), LDA, 
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)

               CALL Recv(NB, 1, LEFT, TMPARR, LDTMP, MYROW, MYCOL,
     $              NPCOL, NPROW, CTX)

               CALL DROT(NB, A( bl_index_r, bl_index_c),
     $              1, TMPARR, 1,
     $              CSCOL(IND, I_c - 1), CSCOL(IND + 1, I_c - 1))

*              ..
*              .. Preform rowrotation
*              ..
               bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
               IF ((I_r .GT. JC) .AND. (MOD(I_r,NB) .NE. MOD(JC,NB))) 
     $              THEN
                  t_I_r = I_r + 1
                  t_bl_index_r = INDXG2L(t_I_r, NB, 0, 0, NPROW)

                  CALL DROT(NNB1, A(bl_index_r, bl_index_c), LDA, 
     $                 A(t_bl_index_r, bl_index_c), LDA,
     $                 CSROW(IND, I_r), CSROW(IND + 1, I_r))
               ENDIF
               

               bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)   
               
               IF ((I_r - NB) .LE. JC) THEN
                  
                  bl_index_r = INDXG2L(JC + 1, NB, 0, 0, NPROW)
                  IF ((I_r - JC) .GT. 1) THEN
                     CALL rRot(I_r - JC, NNB1,
     $                    A(bl_index_r, bl_index_c), LDA,
     $                    CSROW(IND, JC + 1), LDC, -1)
                  ENDIF
                  bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)
                  
               ELSE
               
                  CALL rRot(NB, NNB1,
     $                 A(bl_index_r, bl_index_c), LDA,
     $                 CSROW(IND, I_r - NB + 1), LDC,
     $                 MOD(JC, NB))
                  
               ENDIF
            ENDIF            
 2       CONTINUE  
         GOTO 999
      ENDIF

      IF (NPCOL .EQ. 1) THEN
*        ..
*        .. The processor who owns A(MAXI,I)
*        ..
         IF (MYROW .EQ. bl_r) THEN
*           ..
*           .. Preform colrotation
*           ..
            I_r = MAXI - MOD(MAXI - 1, NB)
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)

            I_c = I + NNB2 - 1 
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            IF (NNB1 .GT. 0) THEN
               t_I_c = I + NNB2
               t_bl_index_c = INDXG2L(t_I_c, NB, 0, 0, NPCOL)

               CALL cRot(NNB1, MOD(MAXI - 1, NB) + 1,
     $              A(bl_index_r, t_bl_index_c), LDA, CSCOL(IND, t_I_c),
     $              LDC)

               CALL DROT(MOD(MAXI - 1, NB) + 1,
     $              A(bl_index_r, t_bl_index_c), 1, 
     $              A(bl_index_r, bl_index_c), 1,
     $              CSCOL(IND, t_I_c - 1), CSCOL(IND + 1, t_I_c - 1))

            ENDIF


            I_c = I
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
            CALL cRot(NNB2, MOD(MAXI - 1, NB) + 1,
     $           A(bl_index_r, bl_index_c), LDA, CSCOL(IND, I_c),
     $           LDC)


*           ..
*           .. Preform rowrotation
*           ..
            IF (MOD(MAXI - 1, NB) .GT. 0) THEN
               IF (I_r .LE. JC) THEN

                  bl_index_r = INDXG2L(JC + 1, NB, 0, 0, NPROW)
                  IF ((MAXI - JC) .GT. 1) THEN
                     CALL rRot(MAXI - JC, NNB2,
     $                    A(bl_index_r, bl_index_c), LDA,
     $                    CSROW(IND, JC + 1), LDC, -1)
                  
                     IF (NNB1 .GT. 0) THEN
                        CALL rRot(MAXI - JC, NNB1,
     $                       A(bl_index_r, t_bl_index_c), LDA,
     $                       CSROW(IND, JC + 1), LDC, -1)
                     ENDIF
                  ENDIF
                  bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
                  
               ELSE

                  CALL rRot(MOD(MAXI - 1, NB) + 1, NNB2,
     $                 A(bl_index_r, bl_index_c), LDA,
     $                 CSROW(IND, I_r), LDC,
     $                 MOD(JC, NB))
                  IF (NNB1 .GT. 0) THEN
                     CALL rRot(MOD(MAXI - 1, NB) + 1, NNB1,
     $                    A(bl_index_r, t_bl_index_c), LDA,
     $                    CSROW(IND, I_r), LDC,
     $                    MOD(JC, NB))
                  ENDIF
                  
               ENDIF
            ENDIF


*           ..
*           .. Send updated rowvalues to next block above
*           ..
            IF (((I_r - 1) .GT. JC) .AND.
     $           (MOD(I_r - 1, NB) .NE. MOD(JC, NB))) THEN
               IF (NNB1 .GT. 0) THEN
                  CALL Send(1, NNB1, UP, A(bl_index_r, t_bl_index_c),
     $                 LDA, MYROW, MYCOL, NPCOL, NPROW, CTX)

                  CALL Recvtmp(1, NNB1, UP, TMPARR, LDTMP,
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)
                     
                  CALL DROT(NNB1, TMPARR, 1, 
     $                 A(bl_index_r, t_bl_index_c), LDA, 
     $                 CSROW(IND, I_r - 1), CSROW(IND + 1, I_r - 1))
               ENDIF
               CALL Send(1, NNB2, UP, A(bl_index_r, bl_index_c),
     $              LDA, MYROW, MYCOL, NPCOL, NPROW, CTX)

               CALL Recvtmp(1, NNB2, UP, TMPARR, LDTMP,
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)
            
               CALL DROT(NNB2, TMPARR, 1, A(bl_index_r, bl_index_c),
     $              LDA, CSROW(IND, I_r - 1), CSROW(IND + 1, I_r - 1))
            ENDIF
         ENDIF

         DO 3 I_r = MAXI - MOD(MAXI - 1, NB) - 1, NB, -NB
            bl_r = INDXG2P(I_r, NB, 0, 0, NPROW)
            IF (MYROW .EQ. bl_r) THEN
*              ..
*              .. Preform colrotation
*              ..
               bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)
               I_c = I + NNB2 - 1
               bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

               IF (NNB1 .GT. 0) THEN
                  t_I_c = I + NNB2
                  t_bl_index_c = INDXG2L(t_I_c, NB, 0, 0, NPCOL)
                  CALL cRot(NNB1, NB,
     $                 A(bl_index_r, t_bl_index_c), LDA,
     $                 CSCOL(IND, t_I_c),
     $                 LDC)

                  CALL DROT(NB, A(bl_index_r, t_bl_index_c), 1,
     $                 A(bl_index_r, bl_index_c), 1,
     $                 CSCOL(IND, t_I_c - 1),
     $                 CSCOL(IND + 1, t_I_c - 1))
               ENDIF
               

               I_c = I
               bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
               CALL cRot(NNB2, NB, A(bl_index_r, bl_index_c), LDA,
     $              CSCOL(IND, I_c), LDC)

*              ..
*              .. Preform rowrotation
*              ..
               bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)  
               IF ((I_r .GT. JC) .AND. (MOD(I_r, NB) .NE. MOD(JC, NB))) 
     $              THEN
                  IF (NNB1 .GT. 0) THEN
                     CALL Send(1, NNB1, DOWN,
     $                    A(bl_index_r, t_bl_index_c), LDA,
     $                    MYROW, MYCOL, NPCOL, NPROW, CTX) 
         
                     CALL Recvtmp(1, NNB1, DOWN, 
     $                    TMPARR, LDTMP,
     $                    MYROW, MYCOL, NPCOL, NPROW, CTX)

                     CALL DROT(NNB1, A(bl_index_r, t_bl_index_c), LDA,
     $                    TMPARR, 1,
     $                    CSROW(IND, I_r), CSROW(IND + 1, I_r))
                  ENDIF
                  CALL Send(1, NNB2, DOWN, A(bl_index_r,bl_index_c),
     $                 LDA, MYROW, MYCOL, NPCOL, NPROW, CTX)

                  CALL Recvtmp(1, NNB2, DOWN, TMPARR, LDTMP,
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)

                  CALL DROT(NNB2, A(bl_index_r, bl_index_c), LDA,
     $                 TMPARR, 1, CSROW(IND, I_r), CSROW(IND + 1, I_r))
               ENDIF
               

               bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)   
               
               IF ((I_r - NB) .LE. JC) THEN                  
                  bl_index_r = INDXG2L(JC + 1, NB, 0, 0, NPROW)
                  IF ((I_r - JC) .GT. 1) THEN
                     CALL rRot(I_r - JC, NNB2,
     $                    A(bl_index_r, bl_index_c), LDA,
     $                    CSROW(IND, JC + 1), LDC, -1)
                  
                     IF (NNB1 .GT. 0) THEN
                        CALL rRot(I_r - JC, NNB1,
     $                       A(bl_index_r, t_bl_index_c), LDA,
     $                       CSROW(IND, JC + 1), LDC, -1)
                     ENDIF
                  ENDIF
                  bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)                  
               ELSE
                  CALL rRot(NB, NNB2,
     $                 A(bl_index_r, bl_index_c), LDA,
     $                 CSROW(IND, I_r - NB + 1), LDC,
     $                 MOD(JC, NB))
                  IF (NNB1 .GT. 0) THEN
                     CALL rRot(NB, NNB1,
     $                    A(bl_index_r, t_bl_index_c), LDA,
     $                    CSROW(IND, I_r - NB + 1), LDC,
     $                    MOD(JC, NB))
                  ENDIF 
               ENDIF
*              ..
*              .. Send updated rowvalues to next block above
*              .. 
               IF (((I_r - NB) .GT. JC) .AND.
     $              (MOD(I_r, NB) .NE. MOD(JC, NB))) THEN
                  IF (NNB1 .GT. 0) THEN
                     CALL Send(1, NNB1, UP, A(bl_index_r, t_bl_index_c),
     $                    LDA, MYROW, MYCOL, NPCOL, NPROW, CTX)

                     CALL Recvtmp(1, NNB1, UP, TMPARR, LDTMP,
     $                    MYROW, MYCOL, NPCOL, NPROW, CTX)                     

                     CALL DROT(NNB1, TMPARR, 1, 
     $                    A(bl_index_r, t_bl_index_c),
     $                    LDA, CSROW(IND, I_r - NB), 
     $                    CSROW(IND + 1, I_r - NB))                     
                  ENDIF
                  CALL Send(1, NNB2, UP, A(bl_index_r, bl_index_c),
     $                 LDA, MYROW, MYCOL, NPCOL, NPROW, CTX)

                  CALL Recvtmp(1, NNB2, UP, TMPARR, LDTMP,
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)
             
                  CALL DROT(NNB2, TMPARR, 1, A(bl_index_r, bl_index_c),
     $                 LDA, CSROW(IND, I_r - NB), 
     $                 CSROW(IND + 1, I_r - NB))
               ENDIF               
            ENDIF
 3       CONTINUE
         GOTO 999
      ENDIF

*     ..
*     .. Else we have at least a 2*2 processor grid
*     ..
      bl_r = INDXG2P(MAXI, NB, 0, 0, NPROW)
      bl_c = INDXG2P(I, NB, 0, 0, NPCOL)
      bl_right_c = INDXG2P(I + NB1 - 1, NB, 0, 0, NPCOL)

               
      IF ((MYROW.EQ.bl_r).AND.(MYCOL.EQ.bl_c)) THEN
*        ..
*        .. Preform colrotation
*        ..
         I_r = MAXI - MOD(MAXI - 1, NB)
         bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)

         I_c = I + NNB2 - 1 
         bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

         IF (NNB1 .GT. 0) THEN
            CALL Send(MOD(MAXI - 1,NB) + 1, 1, RIGHT, 
     $           A(bl_index_r, bl_index_c),
     $           LDA, MYROW, MYCOL, NPCOL, NPROW, CTX)

            CALL Recv(MOD(MAXI - 1, NB) + 1, 1, RIGHT, 
     $           TMPARR, LDTMP,
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)

            CALL DROT(MOD(MAXI - 1, NB) + 1, 
     $           TMPARR, 1,
     $           A(bl_index_r, bl_index_c), 1, 
     $           CSCOL(IND, I_c), CSCOL(IND + 1, I_c))
         ENDIF

         I_c = I
         bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
         CALL cRot(NNB2, MOD(MAXI - 1, NB) + 1,
     $        A(bl_index_r, bl_index_c), LDA, CSCOL(IND, I_c),
     $        LDC) 

      ENDIF
      
      IF ((MYROW.EQ.bl_r).AND.(MYCOL.EQ.bl_right_c).AND.(MYCOL.NE.bl_c)) 
     $     THEN
*        ..
*        .. Preform colrotation
*        ..
         I_r = MAXI - MOD(MAXI - 1, NB)
         bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)

         I_c = I + NNB2
         bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

         CALL cRot(NNB1, MOD(MAXI - 1, NB) + 1,
     $        A(bl_index_r, bl_index_c), LDA, CSCOL(IND, I_c),
     $        LDC)

         CALL Send(MOD(MAXI - 1, NB) + 1, 1, LEFT, 
     $        A(bl_index_r, bl_index_c), LDA,
     $        MYROW, MYCOL, NPCOL, NPROW, CTX)

         CALL Recv(MOD(MAXI - 1, NB) + 1, 1, LEFT, TMPARR, LDTMP,
     $        MYROW, MYCOL, NPCOL, NPROW, CTX)

         CALL DROT(MOD(MAXI - 1, NB) + 1, A( bl_index_r, bl_index_c),
     $        1, TMPARR, 1,
     $        CSCOL(IND, I_c - 1), CSCOL(IND + 1, I_c - 1))

      ENDIF


      DO 10 I_r = MAXI - MOD(MAXI - 1, NB) - 1, NB, -NB                  
         bl_r = INDXG2P(I_r, NB, 0, 0, NPROW)

         IF ((MYROW.EQ.bl_r).AND.(MYCOL.EQ.bl_c)) THEN
*           ..
*           .. Preform colrotation
*           ..
            bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)  
            I_c = I + NNB2 - 1                     
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            IF (NNB1 .GT. 0) THEN
               CALL Send(NB, 1, RIGHT, A(bl_index_r, bl_index_c), LDA,
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)

               CALL Recv(NB, 1, RIGHT, TMPARR, LDTMP, 
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)

               CALL DROT(NB, 
     $              TMPARR, 1,
     $              A(bl_index_r, bl_index_c), 1, 
     $              CSCOL(IND, I_c), CSCOL(IND + 1, I_c))
            ENDIF

            I_c = I
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
            CALL cRot(NNB2, NB, A(bl_index_r, bl_index_c), LDA, 
     $           CSCOL(IND,I_c), LDC)

         ENDIF
                  
         IF ((MYROW.EQ.bl_r).AND.(MYCOL.EQ.bl_right_c).AND.
     $        (MYCOL.NE.bl_c)) THEN
*           ..
*           .. Preform colrotation
*           ..
            I_c = I + NNB2
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
            bl_index_r = INDXG2L(I_r - NB + 1, NB, 0, 0, NPROW)

            CALL cRot(NNB1, NB, A(bl_index_r, bl_index_c), LDA, 
     $           CSCOL(IND, I_c), LDC)

            CALL Send(NB, 1, LEFT, A(bl_index_r, bl_index_c), LDA, 
     $           MYROW, MYCOL,
     $           NPCOL, NPROW, CTX)

            CALL Recv(NB, 1, LEFT, TMPARR, LDTMP, MYROW, MYCOL,
     $           NPCOL, NPROW, CTX)

            CALL DROT(NB, A(bl_index_r, bl_index_c),
     $           1, TMPARR, 1,
     $           CSCOL(IND, I_c - 1), CSCOL(IND + 1, I_c - 1))
         ENDIF
 10   CONTINUE

 999  RETURN
      END

      SUBROUTINE pMVBA(A,LDA,CSROW,CSCOL,LDC,I,NB,NB1,N,NPROW,NPCOL,IND,
     $     TMPARR,LDTMP,COMPQ,MYROW,MYCOL,CTX)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            I, NB, NB1, NPROW, NPCOL, LDA, LDC, IND, LDTMP
      INTEGER            N, MYROW, MYCOL, CTX
	LOGICAL			   COMPQ

*     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*)
      DOUBLE PRECISION   CSCOL(LDC,*), CSROW(LDC,*)
      DOUBLE PRECISION   TMPARR(LDTMP)

*     .. Parameters ..
      DOUBLE             PRECISION ZERO, TEMP
      PARAMETER          ( ZERO = 0.0D+0 )
      INTEGER            DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
      INTEGER            ROWWISE, COLUMNWISE,ALL
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2,ALL=3)

*     .. Local scalars ..
      INTEGER            bl_owner_r, bl_owner_c, bl_right_c, bl_bottom_r
      INTEGER            I_c, I_r, bl_index_r, bl_index_c
      INTEGER            t_I_c, t_I_r, t_bl_index_r, t_bl_index_c
	INTEGER			   bl_right_c2
      INTEGER            NNB_RT, NNB_DW, NNB_LT, NNB_UP, NNB3
*     .. External Functions ..
      EXTERNAL INDXG2L, INDXG2P
      INTEGER INDXG2L,INDXG2P

*     ..
*     .. Start execution
*     ..
      NNB_DW = MOD(I - 1, NB) - (NB - NB1)
      NNB_RT = MOD(I - 1, NB) 
      NNB_LT = NB - MOD(I - 1, NB)
      NNB_UP = NB - MOD(I - 1, NB)

      IF (NNB_DW .LT. 0) THEN
         NNB_UP = NNB_UP + NNB_DW
	   NNB_DW = 0
      ENDIF


      bl_owner_c = INDXG2P(I, NB, 0, 0, NPCOL)
      bl_right_c = INDXG2P(I+ NB - 1, NB, 0, 0, NPCOL)
      bl_bottom_r = INDXG2P(I + NB + NB1 - 1, NB, 0, 0, NPROW)      
      bl_owner_r = INDXG2P(I + NB, NB, 0, 0, NPROW)
      bl_right_c2 = INDXG2P(I+ NB, NB, 0, 0, NPCOL)

*     ..
*     .. If NPROW = 1 then NPCOL > 1 and vice versa
*     ..
      IF (NPROW .EQ. 1) THEN
*        ..
*        .. The processor who owns A(I+NB,I)
*        ..
         IF (MYCOL .EQ. bl_owner_c) THEN
            I_r = I + NB
            I_c = I + NNB_LT - 1

            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            IF (MOD(I - 1, NB) .NE. 0) THEN
*              ..
*              .. Send elements to right for colupdate
*              ..
               t_I_r = I + NB
               t_I_r = t_I_r + NNB_UP

               t_bl_index_r = INDXG2L(t_I_r, NB, 0, 0, NPROW)

               CALL Send(NNB_UP, 1, RIGHT, A(bl_index_r, bl_index_c),
     $              LDA, MYROW, MYCOL, NPCOL, NPROW, CTX)

               CALL Recv(NNB_UP, 1, RIGHT, TMPARR, LDTMP, MYROW, MYCOL,
     $              NPCOL, NPROW, CTX)

               CALL DROT(NNB_UP, TMPARR, 1, 
     $              A( bl_index_r,bl_index_c), 1, 
     $              CSCOL(IND, I_c), CSCOL(IND + 1, I_c))

               IF (NNB_DW .GT. 0) THEN
                  CALL Recv(1, 1, RIGHT, TMPARR, 1,
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)
                  CALL DROT(1, TMPARR, 1,
     $                 A(t_bl_index_r, bl_index_c), 1, 
     $                 CSCOL(IND, I_c),CSCOL(IND + 1, I_c))
               ENDIF

            ENDIF

            I_r = I + NB
            I_r = I_r + NNB_UP - 1
            I_c = I + NNB_LT - 1

            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            IF ((MOD(I - 1, NB) .NE. 0) .AND. (NNB_DW .GT. 0)) THEN
*              .. 
*              .. Preform elimination and calculate new rowvalues
*              ..
               TEMP = A(bl_index_r, bl_index_c)
               CALL DLARTG (TEMP, A(t_bl_index_r, bl_index_c),
     $              CSROW(IND, t_I_r - 1), CSROW(IND + 1, t_I_r - 1), 
     $              A(bl_index_r, bl_index_c))
               A(t_bl_index_r, bl_index_c) = ZERO
            ENDIF

            IF (NNB_DW .GT. 1) THEN
               CALL BRecv(2, NNB_DW - 1, ROWWISE, CSROW(IND, t_I_r),
     $              LDC, bl_owner_r, bl_right_c, CTX)
            ENDIF
*           ..
*           .. Preform altering colrotation and rowelimination
*           .. 
            I_c = I
            I_r = I + NB
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            IF (NNB_UP .NE. NNB_LT) THEN
               IF (NNB_LT .GT. 1) THEN
                  I_c = I + NNB_UP - 1
                  bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
                  CALL cRot(NNB_LT - NNB_UP + 1, NNB_UP,
     $                 A(bl_index_r, bl_index_c), LDA, CSCOL(IND, I_c),
     $                 LDC)
                  I_c = I
                  bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)               
               ENDIF
               IF (NNB_UP .GT. 1) THEN
                  CALL dcrRot(NNB_UP, NNB_LT, A(bl_index_r, bl_index_c),
     $                 LDA, CSROW(IND, I_r), CSCOL(IND, I_c), LDC)
               ENDIF
            ELSE
               IF (NNB_UP .GT. 1) THEN
                  CALL dcrRot(NNB_UP, NNB_LT, A(bl_index_r, bl_index_c),
     $                 LDA, CSROW(IND, I_r), CSCOL(IND, I_c), LDC)
               ENDIF
            ENDIF
*           ..
*           .. Broadcast cos and sin values along current row
*           ..
            IF (NNB_DW .GT. 0) THEN
               CALL BSend(2, NNB_UP + NNB_DW - 1, ROWWISE, 
     $              CSROW(IND, I_r), LDC, CTX)
            ELSE
               CALL BSend(2, NNB_UP - 1, ROWWISE, CSROW(IND, I_r),
     $              LDC, CTX)     
            ENDIF

         ELSE IF (MYCOL .NE. bl_right_c) THEN
*           ..
*           .. Recv cos and sin values in a broadcast along the current row
*           ..
            I_r = I + NB 
            I_r = I_r + NNB_UP
            IF (NNB_DW .GT. 1) THEN
               CALL BRecv(2, NNB_DW - 1, ROWWISE, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_right_c, CTX)
            ENDIF

            I_r = I + NB
            IF (NNB_DW .GT. 0) THEN
               CALL BRecv(2, NNB_UP + NNB_DW - 1, ROWWISE, 
     $              CSROW(IND, I_r), LDC, bl_owner_r, bl_owner_c, CTX)
            ELSE
               CALL BRecv(2, NNB_UP - 1, ROWWISE, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
            ENDIF
         ENDIF
*        ..
*        .. The processor to the right of the one who owns A(I+NB,I)
*        ..
         IF ((MYCOL .EQ. bl_right_c) .AND. (MYCOL .NE. bl_owner_c)) THEN
            I_r = I + NB
            I_c = I + NNB_LT
            t_I_r = I + NB
            t_I_r = t_I_r + NNB_UP
            t_bl_index_r = INDXG2L(t_I_r, NB, 0, 0, NPROW)
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
*           ..
*           .. Preform colrotation
*           ..
            IF (NNB_DW .NE. NNB_RT) THEN
               IF (NNB_DW .GT. 0) THEN
                  I_c = I + NNB_LT + NNB_DW - 1
                  bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)               
                  CALL cRot(NNB_RT - NNB_DW + 1, NNB_DW,
     $                 A(t_bl_index_r, bl_index_c), LDA, 
     $                 CSCOL(IND, I_c), LDC)
                  I_c = I + NNB_LT
                  bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
               ENDIF             
               IF (NNB_DW .GT. 1) THEN
                  CALL dcrRot(NNB_DW, NNB_RT, 
     $                 A(t_bl_index_r, bl_index_c),
     $                 LDA, CSROW(IND, t_I_r), CSCOL(IND, I_c), LDC)
               ENDIF
            ELSE
               IF (NNB_DW .GT. 1) THEN
                  CALL dcrRot(NNB_DW, NNB_RT, 
     $                 A(t_bl_index_r, bl_index_c),
     $                 LDA, CSROW(IND, t_I_r), CSCOL(IND, I_c), LDC)
               ENDIF
            ENDIF

            IF (NNB_RT .GT. 1) THEN
               CALL cRot(NNB_RT, NNB_UP,
     $              A(bl_index_r, bl_index_c), LDA, CSCOL(IND, I_c),
     $              LDC)
            ENDIF
            CALL Send(NNB_UP, 1, LEFT, A(bl_index_r, bl_index_c), LDA, 
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)

            CALL Recv(NNB_UP, 1, LEFT, TMPARR, LDTMP, MYROW, MYCOL,
     $           NPCOL, NPROW, CTX)    

            CALL DROT(NNB_UP, A(bl_index_r, bl_index_c), 1, 
     $           TMPARR, 1, CSCOL(IND, I_c - 1), 
     $           CSCOL(IND + 1, I_c - 1))


            I_r = I + NB
            I_r = I_r + NNB_UP - 1 

            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)                  

            IF (NNB_DW .GT. 0) THEN
               TMPARR(1) = 0
               CALL Send(1, 1, LEFT, A(t_bl_index_r, bl_index_c), LDA, 
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)

               CALL DROT(1, A(t_bl_index_r, bl_index_c), 1, TMPARR, 1,
     $              CSCOL(IND, I_c - 1), CSCOL(IND + 1, I_c - 1))
            ENDIF
*           ..
*           .. Broadcast cos and sin values along current row
*           ..
            IF (NNB_DW .GT. 1) THEN
               CALL BSend(2, NNB_DW - 1, ROWWISE, CSROW(IND, t_I_r),
     $              LDC, CTX)
            ENDIF

*           ..
*           .. Recv cos and sin values sent in an earlier broadcast
*           ..
            I_r = I + NB

            IF (NNB_DW .GT. 0) THEN
               CALL BRecv(2, NNB_UP + NNB_DW - 1, ROWWISE, 
     $              CSROW(IND, I_r), LDC, bl_owner_r, bl_owner_c, CTX)
            ELSE
               CALL BRecv(2, NNB_UP - 1, ROWWISE, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
            ENDIF
*           ..
*           .. Preform rowrotation
*           .. 
            IF (NNB_DW .GT. 0) THEN
               CALL DROT(NNB_RT, A(bl_index_r, bl_index_c), LDA,
     $              A(t_bl_index_r, bl_index_c), LDA,
     $              CSROW(IND, t_I_r - 1), CSROW(IND + 1, t_I_r - 1))
            ENDIF

            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)

            CALL rRot(NNB_UP, NNB_RT,
     $           A(bl_index_r, bl_index_c), LDA, CSROW(IND, I_r),
     $           LDC, -1)
            
         ENDIF
         GOTO 999
      ENDIF

      IF (NPCOL .EQ. 1) THEN
*        ..
*        .. The processor who owns A(I+NB,I)
*        ..
         IF (MYROW .EQ. bl_owner_r) THEN
            I_r = I + NB
            I_c = I + NNB_LT - 1
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            t_I_c = I + NNB_LT
            t_bl_index_c = INDXG2L(t_I_c, NB, 0, 0, NPCOL)

            IF (MOD(I - 1, NB) .NE. 0) THEN
*              ..
*              .. Preform colrotation
*              ..
               IF (NNB_RT .GT. 1) THEN
                  CALL cRot(NNB_RT, NNB_UP,
     $                 A(bl_index_r, t_bl_index_c), LDA, 
     $                 CSCOL(IND, t_I_c), LDC)
               ENDIF
               IF (NNB_RT .GT. 0) THEN
                  CALL DROT(NNB_UP, 
     $                 A(bl_index_r, t_bl_index_c), 1, 
     $                 A(bl_index_r, bl_index_c), 1, 
     $                 CSCOL(IND, t_I_c - 1), CSCOL(IND+1, t_I_c - 1))
               ENDIF
            ENDIF

*           ..
*           .. Exchange values with down to allow that block to be rowrotated
*           ..
            I_r = I + NB + NNB_UP - 1
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)

            IF ((MOD(I - 1, NB) .NE. 0) .AND. (NNB_DW .GT. 0)) THEN
               CALL Send(1, 1, DOWN, A(bl_index_r, bl_index_c), LDA,
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)

               CALL Recv(1, 1, DOWN, TMPARR, LDTMP,
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)                 
*              .. 
*              .. Preform elimination and calculate new cs-rowvalues
*              ..
               TEMP = A(bl_index_r, bl_index_c)
               CALL DLARTG (TEMP, TMPARR,
     $              CSROW(IND, I_r), CSROW(IND + 1, I_r), 
     $              A(bl_index_r, bl_index_c))

               CALL Send(1, NNB_RT, DOWN, A(bl_index_r, t_bl_index_c),
     $              LDA, MYROW, MYCOL, NPCOL, NPROW, CTX)

               CALL Recvtmp(1, NNB_RT, DOWN, TMPARR, LDTMP,
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)  

               CALL DROT(NNB_RT, A(bl_index_r, t_bl_index_c), LDA,
     $              TMPARR, 1, CSROW(IND, I_r), CSROW(IND + 1, I_r))              
            ENDIF

*           ..
*           .. Preform altering colrotation and rowelimination
*           .. 
            I_c = I
            I_r = I + NB
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            IF (NNB_UP .NE. NNB_LT) THEN
               IF (NNB_LT .GT. 1) THEN
                  I_c = I + NNB_UP - 1
                  bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
                  CALL cRot(NNB_LT - NNB_UP + 1, NNB_UP,
     $                 A(bl_index_r, bl_index_c), LDA, CSCOL(IND, I_c),
     $                 LDC)
                  I_c = I
                  bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)               
               ENDIF
               IF (NNB_UP .GT. 1) THEN
                  CALL dcrRot(NNB_UP, NNB_LT, A(bl_index_r, bl_index_c),
     $                 LDA, CSROW(IND, I_r), CSCOL(IND, I_c), LDC)
               ENDIF
            ELSE
               IF (NNB_UP .GT. 1) THEN
                  CALL dcrRot(NNB_UP, NNB_LT, A(bl_index_r, bl_index_c),
     $                 LDA, CSROW(IND, I_r), CSCOL(IND, I_c), LDC)
               ENDIF
            ENDIF
            IF (NNB_UP .GT. 1) THEN
               CALL rRot(NNB_UP, NNB_RT,
     $              A(bl_index_r, t_bl_index_c), LDA, CSROW(IND, I_r),
     $              LDC, -1)
            ENDIF
            IF ((NNB_DW .GT. 0).AND.(COMPQ)) THEN
               CALL BSend(2, NNB_UP, COLUMNWISE, CSROW(IND, I_r),
     $              LDC, CTX)
            ELSE IF ((NNB_UP .GT. 1).AND.(COMPQ)) THEN
			 CALL BSend(2, NNB_UP - 1, COLUMNWISE, CSROW(IND, I_r),
     $              LDC, CTX)  
            ENDIF
	      IF (COMPQ) THEN	 
		     I_r = I + NB 
               I_r = I_r + NNB_UP

               IF (NNB_DW .GT. 0) THEN
               CALL BRecv(2, NNB_DW, COLUMNWISE, CSROW(IND, I_r - 1),
     $              LDC, bl_bottom_r, bl_owner_c, CTX)
               ENDIF
  	      END IF


	   ELSE IF ((MYROW .NE. bl_bottom_r).AND. 
     $		(COMPQ)) THEN
            I_r = I + NB
            IF (NNB_UP .GT. 0) THEN
              IF (NNB_DW .GT. 0) THEN
               CALL BRecv(2, NNB_UP, COLUMNWISE, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
              ELSE IF (NNB_UP .GT. 1) THEN
               CALL BRecv(2, NNB_UP - 1, COLUMNWISE, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
              ENDIF
            ENDIF
         ENDIF
*        ..
*        .. The processor below the one who owns A(I+NB,I)
*        ..
         IF ((MYROW.EQ.bl_bottom_r).AND.(MYROW.NE.bl_owner_r)) 
     $        THEN
            I_r = I + NB + NNB_UP
            I_c = I + NNB_LT - 1

            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)                  
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)                  

            t_I_c = I + NNB_LT
            t_bl_index_c = INDXG2L(t_I_c, NB, 0, 0, NPCOL)

*           ..
*           .. Preform alterning colrotation and rowelimination
*           ..
            IF (NNB_DW .NE. NNB_RT) THEN
               t_I_c = I + NNB_LT + NNB_DW - 1
               t_bl_index_c = INDXG2L(t_I_c, NB, 0, 0, NPCOL)               
               CALL cRot(NNB_RT - NNB_DW + 1, NNB_DW,
     $              A(bl_index_r, t_bl_index_c), LDA, CSCOL(IND, t_I_c),
     $              LDC)
               t_I_c = I + NNB_LT
               t_bl_index_c = INDXG2L(t_I_c, NB, 0, 0, NPCOL)
               IF (NNB_DW .GT. 1) THEN
                  CALL dcrRot(NNB_DW, NNB_RT, 
     $                 A(bl_index_r, t_bl_index_c), LDA, 
     $                 CSROW(IND, I_r), CSCOL(IND, t_I_c), LDC)
               ENDIF
            ELSE
               IF (NNB_DW .GT. 1) THEN
                  CALL dcrRot(NNB_DW, NNB_RT, 
     $                 A(bl_index_r, t_bl_index_c), LDA, 
     $                 CSROW(IND, I_r), CSCOL(IND, t_I_c), LDC)
               ENDIF
            ENDIF
*           ..
*           .. Preform colrotation on top-left element
*           ..            
            CALL DROT(1, A(bl_index_r, t_bl_index_c), 1, 
     $           A(bl_index_r, bl_index_c), 1,
     $           CSCOL(IND, t_I_c - 1),CSCOL(IND + 1, t_I_c - 1))

*           ..
*           .. Exchange element with above to calculate new cs-rowvalues
*           ..     
            CALL Send(1, 1, UP, A(bl_index_r, bl_index_c), LDA, 
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)
            
            CALL Recv(1, 1, UP, TMPARR, LDTMP, MYROW, MYCOL, 
     $           NPCOL, NPROW, CTX)
*           .. 
*           .. Preform elimination and calculate new cs-rowvalues
*           ..
            TEMP = TMPARR(1)
            CALL DLARTG (TEMP, A(bl_index_r, bl_index_c),
     $           CSROW(IND, I_r - 1), CSROW(IND + 1, I_r - 1), TMPARR)
            A(bl_index_r, bl_index_c) = ZERO            

*           ..
*           .. Exchange elements with above to be able to preform rowrotation on toprow
*           ..
            CALL Send(1, NNB_RT, UP, A(bl_index_r, t_bl_index_c), LDA, 
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)

            CALL Recvtmp(1, NNB_RT, UP, TMPARR, LDTMP, MYROW, MYCOL,
     $           NPCOL, NPROW, CTX)

            CALL DROT(NNB_RT, TMPARR, 1,
     $           A(bl_index_r, t_bl_index_c), LDA,
     $           CSROW(IND, I_r - 1), CSROW(IND + 1, I_r - 1))

            I_r = I + NB
            IF ((NNB_UP .GT. 0).AND.(COMPQ)) THEN
              IF (NNB_DW .GT. 0) THEN
               CALL BRecv(2, NNB_UP, COLUMNWISE, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
              ELSE IF (NNB_UP .GT. 1) THEN
               CALL BRecv(2, NNB_UP - 1, COLUMNWISE, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
              ENDIF
            ENDIF   

            IF ((COMPQ).AND.(NNB_DW .GT. 0))  THEN
			 I_r = I + NB + NNB_UP
		     CALL BSend(2, NNB_DW, COLUMNWISE, CSROW(IND, I_r - 1),
     $           LDC, CTX)
		  END IF		       
          ELSE IF ((MYROW.NE.bl_owner_r).AND.(COMPQ)) THEN
			I_r = I + NB + NNB_UP

			IF (NNB_DW .GT. 0) THEN
                CALL BRecv(2, NNB_DW, COLUMNWISE, CSROW(IND, I_r - 1),
     $             LDC, bl_bottom_r, bl_owner_c, CTX)          
  			END IF
		END IF
		GOTO 999
	ENDIF
	

*     ..
*     .. Else we have at least a 2*2 processor grid
*     ..


*     ..
*     .. The processor who owns A(I+NB,I)
*     ..
      IF ((MYROW.EQ.bl_owner_r).AND.(MYCOL.EQ.bl_owner_c)) 
     $     THEN
         I_r = I + NB
         I_c = I + NNB_LT - 1

         bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
         bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

         IF (NNB_RT .GT. 0) THEN
*           ..
*           .. Exchange elements with right for colupdate
*           ..
            CALL Send(NNB_UP, 1, RIGHT, A(bl_index_r, bl_index_c),
     $           LDA, MYROW, MYCOL, NPCOL, NPROW, CTX)

            CALL Recv(NNB_UP, 1, RIGHT, TMPARR, LDTMP, MYROW, MYCOL,
     $           NPCOL, NPROW, CTX)

            CALL DROT(NNB_UP, TMPARR, 1, 
     $           A( bl_index_r,bl_index_c), 1, 
     $           CSCOL(IND, I_c), CSCOL(IND + 1, I_c))
         ENDIF

*        ..
*        .. Send values downwards to allow that part of the block to be rowrotated
*        ..

         IF (NNB_DW .GT. 0) THEN
		   I_r = I + NB + NNB_UP - 1
		   I_c = I + NNB_LT - 1

		   bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
		   bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

            CALL Send(1, 1, DOWN, A(bl_index_r, bl_index_c), LDA,
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)

            CALL Recv(1, 1, DOWN, TMPARR, LDTMP,
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)
*           .. 
*           .. Preform elimination and calculate new rowvalues
*           ..
            TEMP = A(bl_index_r, bl_index_c)
            CALL DLARTG (TEMP, TMPARR(1),
     $           CSROW(IND, I_r), CSROW(IND + 1, I_r), 
     $           A(bl_index_r, bl_index_c))

         ENDIF

*        ..
*        .. Preform altering colrotation and rowelimination
*        .. 
         I_c = I
         I_r = I + NB
         bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
         bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

         IF (NNB_UP .NE. NNB_LT) THEN
            IF (NNB_LT .GT. 1) THEN
               I_c = I + NNB_UP - 1
               bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
               CALL cRot(NNB_LT - NNB_UP + 1, NNB_UP,
     $              A(bl_index_r, bl_index_c), LDA, CSCOL(IND, I_c),
     $              LDC)
               I_c = I
               bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)               
            ENDIF
            IF (NNB_UP .GT. 1) THEN
               CALL dcrRot(NNB_UP, NNB_LT, A(bl_index_r, bl_index_c),
     $              LDA, CSROW(IND, I_r), CSCOL(IND, I_c), LDC)
            ENDIF
         ELSE
            IF (NNB_UP .GT. 1) THEN
               CALL dcrRot(NNB_UP, NNB_LT, A(bl_index_r, bl_index_c),
     $              LDA, CSROW(IND, I_r), CSCOL(IND, I_c), LDC)
            ENDIF
         ENDIF



*        ..
*        .. Broadcast cos and sin values along current row
*        ..
         IF (NNB_UP .GT. 0) THEN
		  I_r = I + NB
		  IF (COMPQ) THEN
            IF (NNB_DW .GT. 0) THEN	        
               CALL BSend(2, NNB_UP, ALL, CSROW(IND, I_r),
     $              LDC, CTX)
			 I_r = I + NB +NNB_UP
               CALL BRecv(2, NNB_DW, ALL, CSROW(IND, I_r - 1),
     $              LDC, bl_bottom_r, bl_right_c, CTX)

            ELSE IF (NNB_UP .GT. 0) THEN
	         I_r = I_r - MOD(I-1, NB)
			 NNB_UP = MIN(NB, N-I_r+1)
               CALL BSend(2, NNB_UP, ALL, CSROW(IND, I_r),
     $              LDC, CTX)  
            ENDIF
		  
		  
	      ELSE
		  

            IF (NNB_DW .GT. 0) THEN
			 CALL BSend(2, NNB_UP, ROWWISE, CSROW(IND, I_r),
     $              LDC, CTX)

            ELSE IF (NNB_UP .GT. 0) THEN
	         I_r = I_r - MOD(I-1, NB)
			 NNB_UP = MIN(NB, N-I_r+1)
               CALL BSend(2, NNB_UP, ROWWISE, 
     $              CSROW(IND, I_r), LDC, CTX)  
            ENDIF
	
	   	  END IF

         ENDIF

      ELSE IF ((MYROW.EQ.bl_owner_r).AND.(MYCOL.NE.bl_right_c)) THEN
*        ..
*        .. Recv cos and sin values in a broadcast along the current row
*        ..
	   I_r = I + NB
	   IF (COMPQ) THEN

         IF (NNB_UP .GT. 0) THEN
            IF (NNB_DW .GT. 0) THEN
               CALL BRecv(2, NNB_UP, ALL, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
		     I_r = I + NB +NNB_UP
               CALL BRecv(2, NNB_DW, ALL, CSROW(IND, I_r - 1),
     $              LDC, bl_bottom_r, bl_right_c, CTX)

            ELSE IF (NNB_UP .GT. 0) THEN
	         I_r = I_r - MOD(I-1, NB)
			 NNB_UP = MIN(NB, N-I_r+1)
               CALL BRecv(2, NNB_UP, ALL, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
            ENDIF
         ENDIF


	   ELSE 

         IF (NNB_UP .GT. 0) THEN
            IF (NNB_DW .GT. 0) THEN
               CALL BRecv(2, NNB_UP, ROWWISE, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)

            ELSE IF (NNB_UP .GT. 0) THEN
			 I_r = I + NB
	         I_r = I_r - MOD(I-1, NB)
			 NNB_UP = MIN(NB, N-I_r+1)

               CALL BRecv(2, NNB_UP, ROWWISE, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
            ENDIF    	
         ENDIF

	   END IF
    
         
	ELSE IF ((MYCOL.EQ.bl_owner_c).AND.(MYROW .NE. bl_bottom_r)) THEN 
         I_r = I + NB
	   IF (COMPQ) THEN
         IF (NNB_UP .GT. 0) THEN
            IF (NNB_DW .GT. 0) THEN
               CALL BRecv(2, NNB_UP, ALL, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
			 I_r = I + NB +NNB_UP
               CALL BRecv(2, NNB_DW, ALL, CSROW(IND, I_r - 1),
     $              LDC, bl_bottom_r, bl_right_c, CTX)

            ELSE IF (NNB_UP .GT. 0) THEN
	         I_r = I_r - MOD(I-1, NB)
			 NNB_UP = MIN(NB, N-I_r+1)
               CALL BRecv(2, NNB_UP, ALL, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
            ENDIF
         ENDIF
	   ENDIF
	ELSE IF ((MYCOL.NE.bl_owner_c).AND.(MYROW .NE. bl_bottom_r).AND.
     $	(MYCOL.NE.bl_right_c).AND.(MYROW .NE. bl_owner_r)) THEN
         I_r = I + NB
	   IF (COMPQ) THEN
         IF (NNB_UP .GT. 0) THEN
            IF (NNB_DW .GT. 0) THEN
               CALL BRecv(2, NNB_UP, ALL, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
		     I_r = I + NB +NNB_UP
               CALL BRecv(2, NNB_DW, ALL, CSROW(IND, I_r - 1),
     $              LDC, bl_bottom_r, bl_right_c, CTX)

            ELSE IF (NNB_UP .GT. 0) THEN
	         I_r = I_r - MOD(I-1, NB)
			 NNB_UP = MIN(NB, N-I_r+1)
               CALL BRecv(2, NNB_UP, ALL, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
            ENDIF
         ENDIF
	   ENDIF	

      ENDIF

*     ..
*     .. The processor to the right of the one who owns A(I+NB,I)
*     ..
      IF ((MYROW .EQ. bl_owner_r) .AND. (MYCOL .EQ. bl_right_c) .AND.
     $     (MYCOL .NE. bl_owner_c)) THEN

         I_r = I + NB
         I_c = I + NNB_LT
         NNB3 = MIN(N - I_c + 1, NB)

*        ..
*        .. Preform colrotation
*        ..

         bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
         bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

         IF (NNB_RT .GT. 1) THEN
            CALL cRot(NNB_RT, NNB_UP,
     $           A(bl_index_r, bl_index_c), LDA, CSCOL(IND, I_c),
     $           LDC)
         ENDIF

         CALL Send(NNB_UP, 1, LEFT, A(bl_index_r, bl_index_c), LDA, 
     $        MYROW, MYCOL, NPCOL, NPROW, CTX)

         CALL Recv(NNB_UP, 1, LEFT, TMPARR, LDTMP, MYROW, MYCOL,
     $        NPCOL, NPROW, CTX)
                     
         CALL DROT(NNB_UP, A( bl_index_r,bl_index_c),
     $        1, TMPARR, 1, 
     $        CSCOL(IND,I_c - 1), CSCOL(IND+1,I_c - 1))



*        ..
*        .. Send bottomelemnts downwards so that that block can be rowrotated
*        ..
	   IF (NNB_DW .GT. 0) THEN
		  I_r = I + NB + NNB_UP - 1
		  bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)                  
         
            CALL Send(1, NNB3, DOWN, A(bl_index_r, bl_index_c),
     $           LDA, MYROW, MYCOL, NPCOL, NPROW, CTX)
            CALL Recvtmp(1, NNB3, DOWN, TMPARR, LDTMP, 
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)
         ENDIF
*        ..
*        .. Recv cos and sin values sent in an earlier broadcast
*        ..  000401 Missade att s„tta I_r = I + NB!


	   I_r = I + NB
	   IF (COMPQ) THEN

         IF (NNB_UP .GT. 0) THEN
            IF (NNB_DW .GT. 0) THEN
               CALL BRecv(2, NNB_UP, ALL, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
		     I_r = I + NB +NNB_UP
               CALL BRecv(2, NNB_DW, ALL, CSROW(IND, I_r - 1),
     $              LDC, bl_bottom_r, bl_right_c, CTX)

            ELSE IF (NNB_UP .GT. 0) THEN
	         I_r = I_r - MOD(I-1, NB)
			 NNB_UP = MIN(NB, N-I_r+1)
               CALL BRecv(2, NNB_UP, ALL, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
            ENDIF
         ENDIF


	   ELSE 

         IF (NNB_UP .GT. 0) THEN
            IF (NNB_DW .GT. 0) THEN
               CALL BRecv(2, NNB_UP, ROWWISE, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
            ELSE IF (NNB_UP .GT. 0) THEN		
	         I_r = I_r - MOD(I-1, NB)
			 NNB_UP = MIN(NB, N-I_r+1)

               CALL BRecv(2, NNB_UP, ROWWISE, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
            ENDIF    	
         ENDIF

	   END IF

	   NNB_DW = MOD(I - 1, NB) - (NB - NB1)
	   NNB_UP = NB - MOD(I-1,NB) 
	   IF (NNB_DW.LT.0) THEN
		NNB_UP = NNB_UP + NNB_DW
		NNB_DW = 0
	   END IF
*        ..
*        .. Preform rowrotation
*        .. 
         IF (NNB_DW .GT. 0) THEN
            I_r = I + NB + NNB_UP - 1
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)

            CALL DROT(NNB3, A(bl_index_r, bl_index_c), LDA,
     $           TMPARR, 1, CSROW(IND, I_r), CSROW(IND + 1, I_r))
         ENDIF       
         IF (NNB_UP .GT. 1) THEN
            I_r  = I + NB
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            
            CALL rRot(NNB_UP, NNB3,
     $           A(bl_index_r, bl_index_c), LDA, CSROW(IND, I_r),
     $           LDC, -1)
         ENDIF

      ENDIF

*     ..
*     .. The processor below the one who owns A(I+NB,I)
*     ..
      IF ((MYROW .EQ. bl_bottom_r) .AND. (MYCOL .EQ. bl_owner_c) .AND.
     $     (MYROW .NE. bl_owner_r)) THEN
         IF (NNB_DW .GT. 0) THEN
	   I_r = I + NB + NNB_UP
         I_c = I + NNB_LT - 1

         bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
         bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)


	   	
         CALL Recv(1, 1, RIGHT, TMPARR, LDTMP,
     $         MYROW, MYCOL, NPCOL, NPROW, CTX)

         CALL DROT(1, TMPARR, 1,
     $        A(bl_index_r, bl_index_c), 1, 
     $        CSCOL(IND, I_c),CSCOL(IND + 1, I_c))

*        ..
*        .. Exchange elements with above to be able to preform rowrotation
*        ..   
         CALL Send(1, 1, UP, A(bl_index_r, bl_index_c), LDA,
     $         MYROW, MYCOL, NPCOL, NPROW, CTX)

         CALL Recv(1, 1, UP, TMPARR, LDTMP, MYROW, MYCOL, 
     $        NPCOL, NPROW, CTX)

*        .. 
*        .. Preform elimination and calculate new rowvalues
*        ..
         TEMP = TMPARR(1)
         CALL DLARTG (TEMP, A(bl_index_r, bl_index_c),
     $        CSROW(IND, I_r - 1), CSROW(IND + 1, I_r - 1), TMPARR)
         A(bl_index_r, bl_index_c) = ZERO
	   
*        ..
*        .. Send calculated cos and sin values to right
*        ..
         CALL Send(2, 1, RIGHT, CSROW(IND, I_r - 1), LDC, 
     $        MYROW, MYCOL, NPCOL, NPROW, CTX)

	   
         
         ENDIF


	   I_r = I + NB
	   IF (COMPQ) THEN
            CALL BRecv(2, NNB_UP, ALL, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
	      I_r = I + NB +NNB_UP
            CALL BRecv(2, NNB_DW, ALL, CSROW(IND, I_r - 1),
     $              LDC, bl_bottom_r, bl_right_c, CTX)

	   ELSE 

	      I_r = I + NB + NNB_UP
            CALL BRecv(2, NNB_DW, ROWWISE, CSROW(IND, I_r - 1),
     $           LDC, bl_bottom_r, bl_right_c, CTX)
	   END IF

      ENDIF

*     ..
*     .. The processor below and to the right of the one who owns A(I+NB,I)
*     ..
      IF ((MYROW .EQ. bl_bottom_r) .AND. (MYCOL .EQ. bl_right_c) .AND.
     $     (MYROW .NE. bl_owner_r) .AND. (MYCOL .NE. bl_owner_c))  THEN

	   I_c = I + NNB_LT
         I_r = I + NB + NNB_UP

         NNB3 = MIN(N - I_c + 1, NB)

         bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
         bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

*        ..
*        .. Preform alterning colrotation and rowelimination
*        ..
         IF (NNB_DW .NE. NNB_RT) THEN
               I_c = I + NNB_LT + NNB_DW - 1
               bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
               CALL cRot(NNB_RT - NNB_DW + 1, NNB_DW,
     $              A(bl_index_r, bl_index_c), LDA, CSCOL(IND, I_c),
     $              LDC)
               I_c = I + NNB_LT
               bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
               IF (NNB_DW .GT. 1) THEN
                  CALL dcrRot(NNB_DW, NNB3, A(bl_index_r, bl_index_c),
     $                 LDA, CSROW(IND, I_r), CSCOL(IND, I_c), LDC)
               ENDIF
         ELSE
            IF (NNB_DW .GT. 1) THEN
               CALL dcrRot(NNB_DW, NNB3, A(bl_index_r, bl_index_c),
     $              LDA, CSROW(IND, I_r), CSCOL(IND, I_c), LDC)
            ENDIF
         ENDIF
*        ..
*        .. Preform colrotation on top-left element
*        ..
         TMPARR(1) = 0
         CALL Send(1, 1, LEFT, A(bl_index_r, bl_index_c), LDA, 
     $        MYROW, MYCOL, NPCOL, NPROW, CTX)
         
         CALL DROT(1, A(bl_index_r, bl_index_c), 1, TMPARR, 1,
     $        CSCOL(IND, I_c - 1),CSCOL(IND + 1, I_c - 1))

*        ..
*        .. Receive elements from above to be able to preform rowrotation on toprow
*        ..
         CALL Send(1, NNB3, UP, A(bl_index_r, bl_index_c), LDA, 
     $        MYROW, MYCOL, NPCOL, NPROW, CTX)

         CALL Recvtmp(1, NNB3, UP, TMPARR, LDTMP, MYROW, MYCOL,
     $        NPCOL, NPROW, CTX)

*        ..
*        .. Receive cos and sin values from left
*        ..
         CALL Recv(2, 1, LEFT, CSROW(IND, I_r - 1), LDC,
     $        MYROW, MYCOL, NPCOL, NPROW, CTX)

                     
         CALL DROT(NNB3, TMPARR, 1,
     $        A(bl_index_r, bl_index_c), LDA,
     $        CSROW(IND, I_r - 1), CSROW(IND + 1, I_r - 1))
*        ..
*        .. Broadcast cos and sin values along current row
*        ..
         IF (COMPQ) THEN
		  I_r = I + NB
            CALL BRecv(2, NNB_UP, ALL, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
	      I_r = I + NB +NNB_UP
		  CALL BSend(2, NNB_DW, ALL, CSROW(IND, I_r - 1),
     $           LDC, CTX)
	   
	   
	   ELSE 
            CALL BSend(2, NNB_DW, ROWWISE, CSROW(IND, I_r - 1),
     $           LDC, CTX)            
	   ENDIF



      ELSE IF ((MYROW .EQ. bl_bottom_r) .AND. (MYROW .NE. bl_owner_r) 
     $        .AND. (MYCOL.NE.bl_owner_c)) THEN
*        ..
*        .. Receive broadcasted cos and sin values along current row
*        ..
         I_r = I + NB
         IF (COMPQ) THEN
           CALL BRecv(2, NNB_UP, ALL, CSROW(IND, I_r),
     $         LDC, bl_owner_r, bl_owner_c, CTX)

	   	 I_r = I + NB + NNB_UP	 
           CALL BRecv(2, NNB_DW, ALL, CSROW(IND, I_r - 1),
     $           LDC, bl_bottom_r, bl_right_c, CTX)
         ELSE
	     I_r = I + NB + NNB_UP
           CALL BRecv(2, NNB_DW, ROWWISE, CSROW(IND, I_r - 1),
     $           LDC, bl_bottom_r, bl_right_c, CTX)


	   END IF
	ELSE IF ((MYCOL.EQ.bl_right_c).AND.(MYROW.NE.bl_owner_r).AND.
     $	(MYROW.NE.bl_bottom_r).AND.(MYCOL.NE.bl_owner_c)) THEN

         I_r = I + NB
         IF (COMPQ) THEN	
         IF (NNB_UP .GT. 0) THEN
            IF (NNB_DW .GT. 0) THEN
               CALL BRecv(2, NNB_UP, ALL, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
		     I_r = I + NB +NNB_UP
               CALL BRecv(2, NNB_DW, ALL, CSROW(IND, I_r - 1),
     $              LDC, bl_bottom_r, bl_right_c, CTX)

            ELSE IF (NNB_UP .GT. 0) THEN
	         I_r = I_r - MOD(I-1, NB)
			 NNB_UP = MIN(NB, N-I_r+1)
               CALL BRecv(2, NNB_UP, ALL, CSROW(IND, I_r),
     $              LDC, bl_owner_r, bl_owner_c, CTX)
            ENDIF

         ENDIF
	   ENDIF

      ENDIF



*     ..
*     .. Update the remaining cols I:N and rows I:I+NB of A (now that they have the csrow values)
*     ..


	NNB_DW = MOD(I - 1, NB) - (NB - NB1)
	NNB_UP = NB - MOD(I-1,NB) 
	IF (NNB_DW.LT.0) THEN
		NNB_UP = NNB_UP + NNB_DW
		NNB_DW = 0
	END IF

      I_c = I + NB
      IF (MOD(I_c - 1, NB) .NE. 0) THEN
         I_c = I_c + NB - MOD(I_c - 1, NB)
      ENDIF

      IF (I_c .GT. N) THEN
         GOTO 999
      ENDIF

      DO 10 t_I_c = I_c, N, NB
         NNB3 = MIN(N - t_I_c + 1, NB)
         bl_owner_c = INDXG2P(t_I_c, NB, 0, 0, NPCOL)
         IF ((MYROW .EQ. bl_owner_r) .AND. (MYCOL .EQ. bl_owner_c)) THEN

            bl_index_c = INDXG2L(t_I_c, NB, 0, 0, NPCOL)
            IF (NNB_DW .GT. 0) THEN
               I_r = I + NB + NNB_UP - 1
               bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
               CALL Send(1, NNB3, DOWN, A(bl_index_r, bl_index_c),
     $              LDA, MYROW, MYCOL, NPCOL, NPROW, CTX)
         
               CALL Recvtmp(1, NNB3, DOWN, TMPARR, LDTMP,
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)

               CALL DROT(NNB3, A(bl_index_r, bl_index_c), LDA,
     $              TMPARR, 1,
     $              CSROW(IND, I_r), CSROW(IND + 1, I_r))
            ENDIF

            IF (NNB_UP .GT. 1) THEN
               I_r = I + NB
               bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
               CALL rRot(NNB_UP, NNB3, A(bl_index_r, bl_index_c),
     $              LDA, CSROW(IND, I_r), LDC, -1)
            ENDIF
         ENDIF

         IF ((MYROW .EQ. bl_bottom_r) .AND. (MYCOL .EQ. bl_owner_c) 
     $        .AND. (MYROW .NE. bl_owner_r)) THEN
            I_r = I + NB + NNB_UP
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            bl_index_c = INDXG2L(t_I_c, NB, 0, 0, NPCOL)

            IF (NNB_DW .GT. 1) THEN
               CALL rRot(NNB_DW, NNB3, A(bl_index_r, bl_index_c),
     $              LDA, CSROW(IND, I_r), LDC, -1)
            ENDIF

            CALL Send(1, NNB3, UP, A(bl_index_r, bl_index_c), LDA, 
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)
            CALL Recvtmp(1, NNB3, UP, TMPARR, LDTMP,
     $           MYROW, MYCOL, NPCOL, NPROW, CTX)

            CALL DROT(NNB3, TMPARR, 1,
     $           A(bl_index_r, bl_index_c), LDA,
     $           CSROW(IND, I_r - 1), CSROW(IND + 1, I_r - 1))
         ENDIF

 10   CONTINUE

 999  RETURN
      END


      SUBROUTINE pKILLCOL(A,LDA,CSROW,CSCOL,LDC,J,L,NB,N,NPROW,NPCOL,
     $     IND,TMPARR,LDTMP,COMPQ,MYROW,MYCOL,CTX)
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER J, N, NB, NPROW, NPCOL, LDA, LDC, LDTMP, L, IND,
     $     MYROW, MYCOL, CTX
	LOGICAL COMPQ
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*)
      DOUBLE PRECISION CSCOL(LDC,*), CSROW(LDC,*)
      DOUBLE PRECISION TMPARR(LDTMP)

*     .. Parameters ..
      DOUBLE PRECISION   ZERO, TEMP
      PARAMETER          ( ZERO = 0.0D+0 )
      INTEGER            DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
      INTEGER            ROWWISE, COLUMNWISE,ALL
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2,ALL=3)

*     .. Local scalars ..
      INTEGER            bl_owner_r, bl_owner_c, bl_bottom_r,bl_right_c
      INTEGER            I_c, I_r, bl_index_r, bl_index_c, LL
      INTEGER            t_I_c, t_I_r, t_bl_index_r
      INTEGER            NNB1, NNB2, NNB3
*     .. External Functions ..
      EXTERNAL           INDXG2L, INDXG2P
      INTEGER            INDXG2L, INDXG2P

*     ..
*     .. Start execution
*     ..
      I_c = J
      LL = L-J+1

      NNB1 = NB - MOD(J, NB)
      NNB2 = LL - NNB1

      IF (NNB2 .LT. 0) THEN
         NNB1 = NNB1 + NNB2
      ENDIF


      bl_owner_c = INDXG2P(J, NB, 0, 0, NPCOL)
	bl_right_c = INDXG2P(J+1, NB, 0, 0, NPCOL)
      bl_owner_r = INDXG2P(J + 1, NB, 0, 0, NPROW)
      bl_bottom_r = INDXG2P(J + LL, NB, 0, 0, NPROW)

      NNB3 = MIN(NB - MOD(I_c - 1, NB), N - I_c + 1)
      IF ((NPROW .EQ. 1) .OR. (NPCOL .EQ. 1)) THEN
         NNB3 = 1
      ENDIF

*     ..
*     .. the J+LL part
*     ..
      IF ((MYROW .EQ. bl_bottom_r) .AND. (MYCOL .EQ. bl_owner_c)
     $     .AND. (MYROW .NE. bl_owner_r)) THEN

         I_r = J + 1 + NNB1

         bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
         bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)

         IF (NNB2 .GT. 1) THEN
            CALL KILLCOL(NNB2, A(bl_index_r, bl_index_c), LDA,
     $           CSROW(IND, I_r), LDC)
         ENDIF

         IF ((NNB3 .GT. 1) .AND. (NNB2 .GT. 1)) THEN
            I_c = J + 1
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)            
            CALL rRot(NNB2, NNB3 - 1,
     $           A(bl_index_r, bl_index_c), LDA, CSROW(IND, I_r),
     $           LDC, -1)
            I_c = J
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
         ENDIF

         CALL Send(1, NNB3, UP, A(bl_index_r, bl_index_c), LDA,
     $        MYROW, MYCOL, NPCOL, NPROW, CTX)
         CALL Recvtmp(1, NNB3, UP, TMPARR, LDTMP, MYROW, MYCOL, 
     $        NPCOL, NPROW, CTX)

         TEMP = TMPARR(1)
         CALL DLARTG (TEMP, A(bl_index_r, bl_index_c),
     $        CSROW(IND, I_r - 1), CSROW(IND + 1, I_r - 1), TMPARR)
         A(bl_index_r, bl_index_c) = ZERO


*        ..
*        .. Send csrowvalues in a rowwise broadcast
*        ..
         IF ((NNB2.GT.0).AND.(NPCOL.GT.1).AND.(.NOT.COMPQ)) THEN
            CALL BSend(2, NNB2, ROWWISE, CSROW(IND, I_r - 1),
     $           LDC, CTX)
	   ELSE IF ((NNB2.GT.0).AND.(COMPQ)) THEN
            CALL BSend(2, NNB2, ALL, CSROW(IND, I_r - 1),
     $           LDC, CTX)
         ENDIF





         IF (NNB3 .GT. 1) THEN
            I_c = J + 1
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)            
            CALL DROT(NNB3 - 1, TMPARR(2), 1,
     $           A(bl_index_r, bl_index_c), LDA,
     $           CSROW(IND, I_r - 1), CSROW(IND + 1, I_r - 1))
         ENDIF

	   I_r = J + 1
         IF ((NNB1.GT.0).AND.COMPQ) THEN
               CALL BRecv(2, NNB1, ALL,
     $              CSROW(IND, I_r), LDC, bl_owner_r, bl_owner_c, CTX)
         ENDIF


      ELSE IF ((MYROW .EQ. bl_bottom_r).AND.(MYROW.NE.bl_owner_r))
     $        THEN
*        ..
*        .. Recv csrowvalues from a rowwise broadcast
*        ..         
         I_r = J + 1 + NNB1
         IF ((COMPQ).AND.(NNB2.GT.0)) THEN
            CALL BRecv(2, NNB2, ALL, CSROW(IND, I_r - 1),
     $           LDC, bl_bottom_r, bl_owner_c, CTX)
	   ELSE IF ((NNB2 .GT. 0) .AND. (NPCOL .GT. 1)) THEN
            CALL BRecv(2, NNB2, ROWWISE, CSROW(IND, I_r - 1),
     $           LDC, bl_bottom_r, bl_owner_c, CTX)
         ENDIF

	   I_r = J + 1
         IF ((NNB1.GT.0).AND.COMPQ) THEN
               CALL BRecv(2, NNB1, ALL,
     $              CSROW(IND, I_r), LDC, bl_owner_r, bl_owner_c, CTX)
         ENDIF

	ELSE IF ((MYCOL.EQ.bl_owner_c).AND.(MYROW.NE.bl_owner_r)) THEN
         
	   IF (COMPQ) THEN
		IF (NNB2 .GT. 0) THEN
            I_r = J + 1 + NNB1
	      CALL BRecv(2, NNB2, ALL, CSROW(IND, I_r - 1),
     $           LDC, bl_bottom_r, bl_owner_c, CTX)
		  I_r = J + 1
		  CALL BRecv(2, NNB1, ALL,
     $              CSROW(IND, I_r), LDC, bl_owner_r, bl_owner_c, CTX)
	    ELSE
			I_r = J + 1
	        CALL BRecv(2, NNB1 - 1, ALL,
     $              CSROW(IND, I_r), LDC, bl_owner_r, bl_owner_c, CTX)
		END IF

         ENDIF
	ELSE IF ((MYCOL.NE.bl_owner_c).AND.(MYROW.NE.bl_owner_r).AND.
     $	(MYROW.NE.bl_bottom_r)) THEN
         IF (COMPQ) THEN
     	    IF (NNB2 .GT. 0) THEN
	      I_r = J + 1 + NNB1
            CALL BRecv(2, NNB2, ALL, CSROW(IND, I_r - 1),
     $           LDC, bl_bottom_r, bl_owner_c, CTX)
		  I_r = J + 1
		  CALL BRecv(2, NNB1, ALL,
     $              CSROW(IND, I_r), LDC, bl_owner_r, bl_owner_c, CTX)
	    ELSE
			I_r = J + 1
	        CALL BRecv(2, NNB1 - 1, ALL,
     $             CSROW(IND, I_r), LDC, bl_owner_r, bl_owner_c, CTX)
		END IF
	   END IF	         
      ENDIF

*     ..
*     .. the J+1 part
*     ..
      IF ((MYROW.EQ.bl_owner_r) .AND. (MYCOL.EQ.bl_owner_c)) THEN        

         bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)
         
         IF (NNB2 .GT. 0) THEN

            I_r = J + 1 + NNB1 - 1
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)

            IF (NPROW .EQ. 1) THEN
               t_I_r = J + 1 + NNB1
               t_bl_index_r = INDXG2L(t_I_r, NB, 0, 0, NPROW)
               IF (NNB2 .GT. 1) THEN
                  CALL KILLCOL(NNB2, A(t_bl_index_r, bl_index_c), LDA,
     $                 CSROW(IND, t_I_r), LDC)
               ENDIF

               TEMP = A(bl_index_r, bl_index_c)
               CALL DLARTG (TEMP, A(t_bl_index_r, bl_index_c),
     $              CSROW(IND, t_I_r - 1), CSROW(IND + 1, t_I_r - 1), 
     $              A(bl_index_r, bl_index_c))
               A(t_bl_index_r, bl_index_c) = ZERO
            ELSE
               CALL Send(1, NNB3, DOWN, A(bl_index_r, bl_index_c), LDA,
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)
               CALL Recvtmp(1, NNB3, DOWN, TMPARR, LDTMP,
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)
               TEMP = A(bl_index_r, bl_index_c)
               CALL DLARTG (TEMP, TMPARR(1),
     $              CSROW(IND, I_r), CSROW(IND + 1, I_r), 
     $              A(bl_index_r, bl_index_c))
            ENDIF
         ENDIF

         I_r = J + 1
         bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
*        ..
*        .. Delete values from bl_index_r to NB
*        ..               
         CALL KILLCOL(NNB1, A(bl_index_r, bl_index_c),
     $        LDA, CSROW(IND, I_r), LDC)

         I_r = J + 1 + NNB1
         IF ((NNB2 .GT. 0) .AND. (NPROW .GT. 1).AND.(COMPQ)) THEN
            CALL BRecv(2, NNB2, ALL, CSROW(IND, I_r - 1),
     $           LDC, bl_bottom_r, bl_owner_c, CTX)
         ENDIF

	   I_r = J + 1
*        ..
*        .. Send csrowvalues in a rowwise broadcast
*        ..
         IF ((.NOT.COMPQ).AND.(NNB1 .GT. 0) .AND. (NPCOL .GT. 1)) THEN
            IF (NNB2 .GT. 0) THEN
               CALL BSend(2, NNB1, ROWWISE, 
     $              CSROW(IND, I_r), LDC, CTX)
            ELSE
               CALL BSend(2, NNB1-1, ROWWISE, 
     $              CSROW(IND, I_r), LDC, CTX)
            ENDIF
	   ELSE IF ((COMPQ).AND.(NNB1.GT.0)) THEN
            IF (NNB2 .GT. 0) THEN
               CALL BSend(2, NNB1, ALL, 
     $              CSROW(IND, I_r), LDC, CTX)
            ELSE
               CALL BSend(2, NNB1 - 1, ALL, 
     $              CSROW(IND, I_r), LDC, CTX)
            ENDIF

         ENDIF

         IF ((NPROW .EQ. 1) .AND. (NNB2 .GT. 0)) THEN
            CALL BSend(2, NNB2, ROWWISE, CSROW(IND, t_I_r - 1),
     $           LDC, CTX)
         ENDIF

         IF (NNB3 .GT. 1) THEN
            I_c = J + 1
            bl_index_c = INDXG2L(I_c, NB, 0, 0, NPCOL)            

            IF (NNB2 .GT. 0) THEN
               I_r = J + 1 + NNB1 - 1
               bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
               CALL DROT(NNB3 - 1, A(bl_index_r, bl_index_c), LDA,
     $              TMPARR(2), 1,
     $              CSROW(IND, I_r), CSROW(IND + 1, I_r))
            ENDIF

            I_r = J + 1
            bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
            CALL rRot(NNB1, NNB3 - 1,
     $           A(bl_index_r, bl_index_c), LDA, CSROW(IND, I_r),
     $           LDC, -1)
            
         ENDIF

      ELSE IF (MYROW.EQ.bl_owner_r) THEN
*        ..
*        .. Recv csrowvalues from a rowwise broadcast
*        ..
         I_r = J + 1 + NNB1
         IF ((COMPQ).AND.(NNB2.GT.0).AND.(NPROW.GT.1)) THEN
            CALL BRecv(2, NNB2, ALL, CSROW(IND, I_r - 1),
     $           LDC, bl_bottom_r, bl_owner_c, CTX)
         ENDIF

         I_r = J + 1
         IF ((.NOT.COMPQ).AND.(NNB1 .GT. 0) .AND. (NPCOL .GT. 1)) THEN
            IF (NNB2 .GT. 0) THEN
               CALL BRecv(2, NNB1, ROWWISE,
     $              CSROW(IND, I_r), LDC, bl_owner_r, bl_owner_c, CTX)
            ELSE
               CALL BRecv(2, NNB1-1, ROWWISE,
     $              CSROW(IND, I_r), LDC, bl_owner_r, bl_owner_c, CTX)
            ENDIF
         ELSE IF (COMPQ.AND.(NNB1 .GT. 0)) THEN
            IF (NNB2 .GT. 0) THEN
               CALL BRecv(2, NNB1, ALL,
     $              CSROW(IND, I_r), LDC, bl_owner_r, bl_owner_c, CTX)
            ELSE
               CALL BRecv(2, NNB1-1, ALL,
     $              CSROW(IND, I_r), LDC, bl_owner_r, bl_owner_c, CTX)
            ENDIF

   	   END IF
	   		

         IF ((NPROW .EQ. 1) .AND. (NNB2 .GT. 0)) THEN
            I_r = J + 1 + NNB1
            CALL BRecv(2, NNB2, ROWWISE, CSROW(IND, I_r - 1),
     $           LDC, bl_owner_r, bl_owner_c, CTX)        
         ENDIF
      ENDIF




      IF ((NPROW .GT. 1) .AND. (NPCOL .GT. 1)) THEN
         I_c = J + 1
         IF (MOD(J, NB) .NE. 0) THEN
            I_c = I_c + NB - MOD(J, NB)
         ENDIF

         IF (I_c .GT. N) THEN
            GOTO 999
         ENDIF

         DO 10 t_I_c = I_c, N, NB
            NNB3 = MIN(N - t_I_c + 1, NB)
            bl_owner_c = INDXG2P(t_I_c, NB, 0, 0, NPCOL)
            IF ((MYROW .EQ. bl_owner_r) .AND. (MYCOL .EQ. bl_owner_c)) 
     $           THEN
               bl_index_c = INDXG2L(t_I_c, NB, 0, 0, NPCOL)
               IF ((MOD(J, NB) .NE. 0) .AND. (NNB2 .GT. 0)) THEN
                  I_r = J + 1 + NNB1 - 1
                  bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
                  CALL Send(1, NNB3, DOWN, A(bl_index_r, bl_index_c),
     $                 LDA, MYROW, MYCOL, NPCOL, NPROW, CTX)
                  
                  CALL Recvtmp(1, NNB3, DOWN, TMPARR, LDTMP,
     $                 MYROW, MYCOL, NPCOL, NPROW, CTX)
               
                  CALL DROT(NNB3, A(bl_index_r, bl_index_c), LDA,
     $                 TMPARR, 1,
     $                 CSROW(IND, I_r), CSROW(IND + 1, I_r))
               ENDIF
               
               I_r = J + 1
               bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
               CALL rRot(NNB1, NNB3, A(bl_index_r, bl_index_c),
     $              LDA, CSROW(IND, I_r), LDC, -1)
            ENDIF
            
            IF ((MYROW .EQ. bl_bottom_r) .AND. (MYCOL .EQ. bl_owner_c) 
     $           .AND. (MYROW .NE. bl_owner_r)) THEN
               I_r = J + 1 + NNB1
               bl_index_r = INDXG2L(I_r, NB, 0, 0, NPROW)
               bl_index_c = INDXG2L(t_I_c, NB, 0, 0, NPCOL)

               IF (NNB2 .GT. 1) THEN
                  CALL rRot(NNB2, NNB3, A(bl_index_r, bl_index_c),
     $                 LDA, CSROW(IND, I_r), LDC, -1)
               ENDIF

               CALL Send(1, NNB3, UP, A(bl_index_r, bl_index_c),
     $              LDA, MYROW, MYCOL, NPCOL, NPROW, CTX)
               
               CALL Recvtmp(1, NNB3, UP, TMPARR, LDTMP,
     $              MYROW, MYCOL, NPCOL, NPROW, CTX)
               
               CALL DROT(NNB3, TMPARR, 1,
     $              A(bl_index_r, bl_index_c), LDA,
     $              CSROW(IND, I_r - 1), CSROW(IND + 1, I_r - 1))
            ENDIF

 10      CONTINUE
      ENDIF
 999  RETURN
      END

	SUBROUTINE UPDABQZ(A, B, Q, Z,
     $     LDA, LDB, LDQ, LDZ,
     $     TMPW, TMPW2,
     $	 LDTMP,
     $     TMPQ, TMPZ,
     $	 LDTMPQ, LDTMPZ,
     $     DESCA, DESCB, DESCQ, DESCZ,
     $     IFIRST, ILASTM,
     $     ILO, NB, N, NN,
     $	 ILQ, ILZ,WSIZE)



      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
*     ..
      INTEGER            LDA, LDB, LDQ, LDZ,LDTMP
      INTEGER            IFIRST, ILASTM
      INTEGER            ILO, NB, N, NN,WSIZE
	INTEGER			   LDTMPQ, LDTMPZ
      INTEGER            MAXIT, RET
      LOGICAL            ILQ, ILZ
*     ..
*     .. Array Arguments ..
*     ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
      DOUBLE PRECISION   Q( LDQ, * ), Z( LDZ, * )

 
	DOUBLE PRECISION   TMPQ(LDTMPQ,*), TMPZ(LDTMPZ,*)
      DOUBLE PRECISION   TMPW(LDTMP, *), TMPW2(LDTMP,*)
      INTEGER            DESCA(*), DESCB(*)
      INTEGER            DESCQ(*), DESCZ(*)


*


*     Purpose
*     =======
*
*
*
*     ..
*     .. Parameters ..
*     ..
      DOUBLE PRECISION   HALF, ZERO, ONE, SAFETY
      PARAMETER          ( HALF = 0.5D+0, ZERO = 0.0D+0, ONE = 1.0D+0,
     $                   SAFETY = 1.0D+2 )
      INTEGER            DOWN, UP, RIGHT, LEFT
      PARAMETER          ( DOWN = 1, UP = 2, RIGHT = 3, LEFT = 4 )
      INTEGER            ROWWISE, COLUMNWISE, ALL
      PARAMETER          ( ROWWISE = 1, COLUMNWISE = 2, ALL=3)
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_, 
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
*     ..
      INTEGER            IACOL, IAROW
      INTEGER            LROW1, LROW2, LROW3
	INTEGER            LCOL1, LCOL2, LCOL3
      INTEGER            IACOL2, IAROW2
      INTEGER            JJ, JC, J, JJC
      INTEGER            JR, JJR, INCJ, INCI
      INTEGER            SCNT, DCNT, RCNT      
	INTEGER			   ICTXT, MYROW, MYCOL, NPROW, NPCOL

     

*     ..
*     .. Local Arrays
*     ..
	

*     ..
*     .. External Subroutines ..
*     ..
      EXTERNAL           BLACS_BARRIER, PDGATTER, PDLACP3,
     $                   DLAG2, DLARTG
*     ..
*     .. Intrinsic Functions ..
*     ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT, MOD
*     ..
*     .. External Functions ..
*     ..
      INTEGER            INDXG2L, INDXG2P
      EXTERNAL           INDXG2L, INDXG2P
      DOUBLE PRECISION   DLAPY3, DLAPY2
      EXTERNAL           DLAPY3, DLAPY2

	ICTXT = DESCA(CTXT_)
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )  



      J = IFIRST

      IAROW = INDXG2P(J, NB, 0, 0, NPROW)
      IACOL = INDXG2P(J, NB, 0, 0, NPCOL)



	IF (NPROW.EQ.1) THEN
		DCNT = 0
	ELSE IF (INDXG2P(J,NB,0,0,NPROW).EQ.
     $	(INDXG2P(J+NN-1,NB,0,0,NPROW))) THEN
		DCNT = 0
	ELSE
		DCNT = MOD(J+NN-1, NB)
		
	END IF

	IF (NPCOL.EQ.1) THEN
		RCNT = 0
	ELSE IF (INDXG2P(J,NB,0,0,NPCOL).EQ.
     $	(INDXG2P(J+NN-1,NB,0,0,NPCOL))) THEN
		RCNT = 0
	ELSE
		RCNT = MOD(J+NN-1, NB)
	    
			
	END IF




	JC = J+WSIZE
*	Update remaing columns of A, B
	IF ((JC.LE.ILASTM).AND.(J.GT.0)) THEN
		SCNT = MIN(NB - MOD(JC-1,NB), ILASTM-JC+1)
		INCJ = NB

		JR = J
		IACOL = INDXG2P(JC, NB, 0, 0, NPCOL)
		LCOL1 = INDXG2L(JC, NB, 0, 0, NPCOL)
		LROW1 = INDXG2L(JR, NB, 0, 0, NPROW)
		IAROW = INDXG2P(JR, NB, 0, 0, NPROW)


*	..  Multiply TMPQ on the left of A, B
		IF (DCNT.EQ.0) THEN 				
			IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
				DO 500 JJC = JC, JC+SCNT-1, INCJ
					LCOL1 = INDXG2L(JJC, NB, 0, 0, NPCOL)
					CALL DGEMM('T','N',NN,SCNT,NN,ONE,
     $					TMPQ,LDTMPQ, A(LROW1,LCOL1),LDA,
     $					ZERO,TMPW,LDTMP)
					CALL DLACPY('A',NN,SCNT,TMPW, LDTMP, 
     $					A(LROW1, LCOL1), LDA)
					CALL DGEMM('T','N',NN,SCNT,NN,ONE,
     $					TMPQ,LDTMPQ,B(LROW1,LCOL1),LDB,
     $					ZERO,TMPW,LDTMP)
					CALL DLACPY('A',NN,SCNT,TMPW, LDTMP, 
     $					B(LROW1, LCOL1), LDB)
 500				CONTINUE
			END IF
			DO 501 JJC=JC+SCNT, ILASTM, INCJ
				SCNT = MIN(NB, ILASTM-JJC+1)

				IACOL = INDXG2P(JJC, NB, 0, 0, NPCOL)
				IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL))THEN
					LCOL1 = INDXG2L(JJC, NB, 0, 0, NPCOL)
					CALL DGEMM('T','N',NN,SCNT,NN,ONE,
     $					TMPQ,LDTMPQ,A(LROW1,LCOL1),LDA,
     $					ZERO,TMPW,LDTMP)
					CALL DLACPY('A',NN,SCNT,TMPW, LDTMP, 
     $					A(LROW1, LCOL1), LDA)

					CALL DGEMM('T','N',NN,SCNT,NN,ONE,
     $					TMPQ,LDTMPQ,B(LROW1, LCOL1),LDB,
     $					ZERO,TMPW,LDTMP)
					CALL DLACPY('A',NN,SCNT,TMPW, LDTMP, 
     $					B(LROW1, LCOL1), LDB) 
				END IF
 501			CONTINUE
		ELSE
			IAROW2 = INDXG2P(JR+NN-1, NB, 0, 0, NPROW)
			LCOL1 = INDXG2L(JC, NB, 0, 0, NPCOL)
			IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN				
				LROW1 = INDXG2L(JR, NB, 0, 0, NPROW)
*				Send upper part of A
				CALL Send(NN-DCNT,SCNT,DOWN,A(LROW1, LCOL1),
     $				 LDA, MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				
				
				CALL DLACPY('A',NN-DCNT,SCNT,A(LROW1, LCOL1),
     $				 LDA, TMPW2, LDTMP)
*				Recv lower part of A				
				CALL Recv(DCNT,SCNT,DOWN,TMPW2(NN-DCNT+1,1), 
     $				 LDTMP, MYROW, MYCOL, NPCOL, NPROW, ICTXT)

											
				CALL DGEMM('T','N',NN,SCNT, NN,ONE,
     $				TMPQ,LDTMPQ, TMPW2,LDTMP,
     $				ZERO,TMPW,LDTMP)


*				Send upper part of B
				CALL Send(NN-DCNT,SCNT,DOWN,B(LROW1, LCOL1),
     $				 LDB, MYROW, MYCOL, NPCOL, NPROW, ICTXT)

				CALL DLACPY('A',NN-DCNT,SCNT,B(LROW1, LCOL1),
     $				 LDB, TMPW2, LDTMP)

*				Restore upper port of A				
				CALL DLACPY('A',NN-DCNT,SCNT,TMPW,
     $				 LDTMP, A(LROW1, LCOL1),LDA )



*				Recv lower part of B				
				CALL Recv(DCNT,SCNT,DOWN,TMPW2(NN-DCNT+1,1), 
     $				 LDTMP, MYROW, MYCOL, NPCOL, NPROW, ICTXT)


				CALL DGEMM('T','N',NN,SCNT, NN,ONE,
     $				TMPQ,LDTMPQ, TMPW2,LDTMP,
     $				ZERO,TMPW,LDTMP)


*				Restore upper port of B
				CALL DLACPY('A',NN-DCNT,SCNT,TMPW, LDTMP,
     $					B(LROW1, LCOL1), LDB)


			
			ELSE IF ((MYROW.EQ.IAROW2).AND.
     $				(MYCOL.EQ.IACOL)) THEN
				LROW1 = INDXG2L(JR+(NN-DCNT), NB, 0,0,NPROW)

	
*				Send lower part of A	
				Call Send(DCNT, SCNT, UP, A(LROW1, LCOL1),
     $				LDA, MYROW, MYCOL, NPCOL, NPROW, ICTXT)

				CALL DLACPY('A',DCNT,SCNT,A(LROW1, LCOL1),
     $				 LDA, TMPW2(NN-DCNT+1,1), LDTMP)

*				Recv upper part of A
				Call Recv(NN-DCNT, SCNT, UP,TMPW2,
     $				LDTMP, MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				
							
				CALL DGEMM('T','N',NN,SCNT, NN,ONE,
     $				TMPQ,LDTMPQ, TMPW2,LDTMP,
     $				ZERO,TMPW,LDTMP)

 
 
*				Send lower part of B
				Call Send(DCNT, SCNT, UP, B(LROW1, LCOL1),
     $				LDB, MYROW, MYCOL, NPCOL, NPROW, ICTXT)

*				Restore A lower part
				CALL DLACPY('A',DCNT,SCNT,TMPW(NN-DCNT+1,1),
     $				 LDTMP, A(LROW1, LCOL1),LDA)



				CALL DLACPY('A',DCNT,SCNT,B(LROW1, LCOL1),
     $				 LDB, TMPW2(NN-DCNT+1,1), LDTMP)

*				Recv upper part of B
				Call Recv(NN-DCNT, SCNT, UP,TMPW2,
     $				LDTMP, MYROW, MYCOL, NPCOL, NPROW, ICTXT)



				CALL DGEMM('T','N',NN,SCNT, NN,ONE,
     $				TMPQ,LDTMPQ, TMPW2,LDTMP,
     $				ZERO,TMPW,LDTMP)

*				Restore B lower part
				CALL DLACPY('A',DCNT,SCNT,TMPW(NN-DCNT+1,1),
     $				 LDTMP, B(LROW1, LCOL1),LDB)


			END IF

			DO 520 JJC=JC+SCNT, ILASTM, NB
				SCNT = MIN(NB, ILASTM-JJC+1)
				IACOL = INDXG2P(JJC, NB, 0, 0, NPCOL)
				LCOL1 = INDXG2L(JJC, NB, 0, 0, NPCOL)

				IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
	
					LROW1 = INDXG2L(JR, NB, 0, 0, NPROW)

*					Send upper part of A
					CALL Send(NN-DCNT,SCNT,DOWN,A(LROW1, LCOL1),
     $					LDA, MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				
				
					CALL DLACPY('A',NN-DCNT,SCNT,A(LROW1, LCOL1),
     $					LDA, TMPW2, LDTMP)
*					Recv lower part of A				
					CALL Recv(DCNT,SCNT,DOWN,TMPW2(NN-DCNT+1,1), 
     $					LDTMP, MYROW, MYCOL, NPCOL, NPROW, ICTXT)

											
					CALL DGEMM('T','N',NN,SCNT, NN,ONE,
     $					TMPQ,LDTMPQ, TMPW2,LDTMP,
     $					ZERO,TMPW,LDTMP)


*					Send upper part of B
					CALL Send(NN-DCNT,SCNT,DOWN,B(LROW1, LCOL1),
     $					 LDB, MYROW, MYCOL, NPCOL, NPROW, ICTXT)

					CALL DLACPY('A',NN-DCNT,SCNT,B(LROW1, LCOL1),
     $					LDB, TMPW2, LDTMP)

*					Restore upper port of A				
					CALL DLACPY('A',NN-DCNT,SCNT,TMPW,
     $					LDTMP, A(LROW1, LCOL1),LDA )

*					Recv lower part of A				
					CALL Recv(DCNT,SCNT,DOWN,TMPW2(NN-DCNT+1,1), 
     $					LDTMP, MYROW, MYCOL, NPCOL, NPROW, ICTXT)


					CALL DGEMM('T','N',NN,SCNT, NN,ONE,
     $					TMPQ,LDTMPQ, TMPW2,LDTMP,
     $					ZERO,TMPW,LDTMP)


*					Restore upper port of B
					CALL DLACPY('A',NN-DCNT,SCNT,TMPW, LDTMP,
     $					B(LROW1, LCOL1), LDB)


				ELSE IF ((MYROW.EQ.IAROW2).AND.
     $					(MYCOL.EQ.IACOL)) THEN	
			
					LROW1 =INDXG2L(JR+(NN-DCNT),NB,0,0,NPROW)
	
*					Send lower part of A	
					Call Send(DCNT, SCNT, UP, A(LROW1, LCOL1),
     $					LDA, MYROW, MYCOL, NPCOL, NPROW, ICTXT)

					CALL DLACPY('A',DCNT,SCNT,A(LROW1, LCOL1),
     $					LDA, TMPW2(NN-DCNT+1,1), LDTMP)

*					Recv upper part of A
					Call Recv(NN-DCNT, SCNT, UP,TMPW2,
     $					LDTMP, MYROW, MYCOL, NPCOL, NPROW, ICTXT)
				
							
					CALL DGEMM('T','N',NN,SCNT, NN,ONE,
     $					TMPQ,LDTMPQ, TMPW2,LDTMP,
     $					ZERO,TMPW,LDTMP)

 
 
*					Send lower part of B
					Call Send(DCNT, SCNT, UP, B(LROW1, LCOL1),
     $					LDB, MYROW, MYCOL, NPCOL, NPROW, ICTXT)

*					Restore A lower part
					CALL DLACPY('A',DCNT,SCNT,TMPW(NN-DCNT+1,1),
     $					 LDTMP, A(LROW1, LCOL1),LDA)



					CALL DLACPY('A',DCNT,SCNT,B(LROW1, LCOL1),
     $					 LDB, TMPW2(NN-DCNT+1,1), LDTMP)

*					Recv upper part of B
					Call Recv(NN-DCNT, SCNT, UP,TMPW2,
     $					LDTMP, MYROW, MYCOL, NPCOL, NPROW, ICTXT)



					CALL DGEMM('T','N',NN,SCNT, NN,ONE,
     $					TMPQ,LDTMPQ, TMPW2,LDTMP,
     $					ZERO,TMPW,LDTMP)

*					Restore B lower part
					CALL DLACPY('A',DCNT,SCNT,TMPW(NN-DCNT+1,1),
     $					LDTMP, B(LROW1, LCOL1),LDB)

				END IF
 520			CONTINUE
		END IF
	END IF
		
*	Update remaining rows of A, B
	IF (J.GT.0) THEN
		JC = J
		IACOL = INDXG2P(JC, NB, 0, 0, NPCOL)
		LCOL1 = INDXG2L(JC, NB, 0, 0, NPCOL)
		INCJ = NB

*		Multiply TMPZ on the right of A, B
		IF (RCNT.EQ.0) THEN
			DO 540 JR = 1, J-1, INCJ

				SCNT = MIN(NB, J-1-JR+1)

				IAROW = INDXG2P(JR, NB, 0, 0, NPROW)
				LROW1 = INDXG2L(JR, NB, 0, 0, NPROW)
					
				IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL))
     $				THEN
					CALL DGEMM('N','N',SCNT,NN,NN,ONE,
     $					A(LROW1, LCOL1),LDA,TMPZ,
     $					LDTMPZ,ZERO,
     $					TMPW,LDTMP)
					CALL DLACPY('A',SCNT,NN,TMPW, LDTMP, 
     $					A(LROW1, LCOL1), LDA)
					CALL DGEMM('N','N',SCNT,NN,NN,ONE,
     $					B(LROW1, LCOL1),LDB,TMPZ,
     $					LDTMPZ,ZERO,
     $					TMPW,LDTMP)
					CALL DLACPY('A',SCNT,NN,TMPW, LDTMP, 
     $					B(LROW1, LCOL1), LDB)
				END IF
 540			CONTINUE
		ELSE
			DO 550 JR = 1, J-1, NB
				SCNT = MIN(NB, J-1-JR+1)
				IAROW = INDXG2P(JR, NB, 0, 0, NPROW)
				IACOL2 = INDXG2P(JC+NN-1, NB, 0, 0, NPCOL)
				IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) THEN
					LCOL1 = INDXG2L(JC, NB, 0, 0, NPCOL)
					LROW1 = INDXG2L(JR, NB, 0, 0, NPROW)
					CALL DLACPY('A',SCNT, NN-RCNT,
     $					A(LROW1, LCOL1),LDA, TMPW2, LDTMP)
		
					CALL Recv(SCNT, RCNT, RIGHT, 
     $					TMPW2(1,NN-RCNT+1), LDTMP, MYROW,
     $					MYCOL, NPCOL, NPROW, ICTXT)

					CALL DGEMM('N','N',SCNT,NN,NN,ONE,
     $						TMPW2,LDTMP,TMPZ,
     $						LDTMPZ,ZERO,
     $						TMPW,LDTMP)
					CALL DLACPY('A',SCNT,NN-RCNT,TMPW, LDTMP, 
     $					A(LROW1, LCOL1),LDA)
					CALL Send(SCNT, RCNT, RIGHT, 
     $					TMPW(1,NN-RCNT+1), LDTMP, MYROW, 
     $					MYCOL, NPCOL, NPROW, ICTXT)

					CALL DLACPY('A',SCNT, NN-RCNT,
     $					B(LROW1, LCOL1),LDB, TMPW2, LDTMP)
			
					CALL Recv(SCNT, RCNT, RIGHT, 
     $					TMPW2(1,NN-RCNT+1), LDTMP, MYROW,
     $					MYCOL, NPCOL, NPROW, ICTXT)

					CALL DGEMM('N','N',SCNT,NN,NN,ONE,
     $						TMPW2,LDTMP,TMPZ,
     $						LDTMPZ,ZERO,
     $						TMPW,LDTMP)


					CALL DLACPY('A',SCNT,NN-RCNT,TMPW, LDTMP, 
     $					B(LROW1, LCOL1),LDB)
		
					CALL Send(SCNT, RCNT, RIGHT, 
     $					TMPW(1,NN-RCNT+1), LDTMP, MYROW, 
     $					MYCOL, NPCOL, NPROW, ICTXT)


				ELSE IF((MYROW.EQ.IAROW).AND.
     $						(MYCOL.EQ.IACOL2)) THEN
					LCOL1 = INDXG2L(JC+(NN-RCNT),NB,0,0,NPCOL)
					LROW1 = INDXG2L(JR, NB, 0, 0, NPROW)
	
					Call Send(SCNT, RCNT, LEFT,A(LROW1,LCOL1),
     $					LDA, MYROW, MYCOL, NPCOL, NPROW,ICTXT)

					Call Recv(SCNT, RCNT, LEFT,A(LROW1,LCOL1),
     $					LDA, MYROW, MYCOL, NPCOL, NPROW,ICTXT)

					Call Send(SCNT, RCNT, LEFT,B(LROW1,LCOL1),
     $					LDB, MYROW, MYCOL, NPCOL, NPROW,ICTXT)

					Call Recv(SCNT, RCNT, LEFT,B(LROW1, LCOL1),
     $					LDB, MYROW, MYCOL, NPCOL,NPROW,ICTXT)

				END IF
 550			CONTINUE

		END IF
	END IF
			
		
*	. Update Q on the right with TMPQ
	IF (ILQ.AND.(J.LE.ILASTM).AND.(J.GT.0)) THEN
		JC = J
		IACOL = INDXG2P(JC, NB, 0, 0, NPCOL)
		LCOL1 = INDXG2L(JC, NB, 0, 0, NPCOL)
		INCJ = NB
		
		IF (RCNT.EQ.0) THEN
			DO 530 JR = 1, N, INCJ
				SCNT = MIN(NB, N-JR+1)

				IAROW = INDXG2P(JR, NB, 0, 0, NPROW)
				LROW1 = INDXG2L(JR, NB, 0, 0, NPROW)
				IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL))THEN
					CALL DGEMM('N','N',SCNT,NN,NN,ONE,
     $					Q(LROW1,LCOL1),LDQ,TMPQ,LDTMPQ,
     $					ZERO,TMPW,LDTMP)
					CALL DLACPY('A',SCNT,NN,TMPW, LDTMP, 
     $					Q(LROW1, LCOL1), LDQ)
				END IF
 530			CONTINUE
		ELSE
			IACOL2 = INDXG2P(JC+NN-1, NB, 0, 0, NPCOL)
			DO 535 JR = 1, N, NB
				SCNT = MIN(NB, N-JR+1)
				IAROW = INDXG2P(JR, NB, 0, 0, NPROW)
				LROW1 = INDXG2L(JR, NB, 0, 0, NPROW)
				IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) 
     $				THEN
					LCOL1 = INDXG2L(JC, NB, 0, 0, NPCOL)
					CALL DLACPY('A',SCNT, NN-RCNT,
     $					Q(LROW1, LCOL1),LDQ, TMPW2, LDTMP)
					CALL Recv(SCNT, RCNT, RIGHT, 
     $					TMPW2(1,NN-RCNT+1), LDTMP, MYROW,
     $					MYCOL, NPCOL, NPROW, ICTXT)
					CALL DGEMM('N','N',SCNT,NN,NN,ONE,
     $					TMPW2,LDTMP,TMPQ,
     $					LDTMPQ,ZERO,
     $					TMPW,LDTMP)
					CALL DLACPY('A',SCNT,NN-RCNT,TMPW, LDTMP, 
     $					Q(LROW1, LCOL1),LDQ)
					CALL Send(SCNT, RCNT, RIGHT, 
     $					TMPW(1,NN-RCNT+1), LDTMP, MYROW, 
     $					MYCOL, NPCOL, NPROW, ICTXT)

				ELSE IF ((MYROW.EQ.IAROW).AND.
     $				(MYCOL.EQ.IACOL2)) THEN
						LCOL1 = INDXG2L(JC+NN-RCNT, NB, 0, 0, NPCOL)
					Call Send(SCNT, RCNT, LEFT,
     $					Q(LROW1,LCOL1),LDQ, MYROW, MYCOL,
     $					NPCOL, NPROW,ICTXT) 
					Call Recv(SCNT, RCNT, LEFT,
     $					Q(LROW1,LCOL1),LDQ, MYROW, MYCOL,
     $					NPCOL, NPROW,ICTXT)
				END IF
  535 		CONTINUE
		END IF
	END IF

*	.. Update Z on the right with TMPZ
	IF (ILZ.AND.(J.LE.ILASTM).AND.(J.GT.0)) THEN
		JC = J
		IACOL = INDXG2P(JC, NB, 0, 0, NPCOL)
		LCOL1 = INDXG2L(JC, NB, 0, 0, NPCOL)
		INCJ = NB


		IF (RCNT.EQ.0) THEN
			DO 560 JR = 1, N, INCJ
				SCNT = MIN(NB, N-JR+1)

				IAROW = INDXG2P(JR, NB, 0, 0, NPROW)
				LROW1 = INDXG2L(JR, NB, 0, 0, NPROW)
				IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL))
     $				THEN
					CALL DGEMM('N','N',SCNT,NN,NN,ONE,
     $					Z(LROW1,LCOL1),LDZ,TMPZ,LDTMPZ,
     $					ZERO,TMPW,LDTMP)
					CALL DLACPY('A',SCNT,NN,TMPW, LDTMP, 
     $					Z(LROW1, LCOL1), LDZ)
				END IF
 560			CONTINUE
		ELSE
			IACOL2 = INDXG2P(JC+NN-1, NB, 0, 0, NPCOL)
			DO 570 JR = 1, N, NB
				SCNT = MIN(NB, N-JR+1)
				LROW1 = INDXG2L(JR, NB, 0, 0, NPROW)
				IAROW = INDXG2P(JR, NB, 0, 0, NPROW)
				IF ((MYROW.EQ.IAROW).AND.(MYCOL.EQ.IACOL)) 
     $				THEN
					LCOL1 = INDXG2L(JC, NB, 0, 0, NPCOL)
					CALL DLACPY('A',SCNT, NN-RCNT,
     $					Z(LROW1, LCOL1),LDZ, TMPW2, LDTMP)
		
					CALL Recv(SCNT, RCNT, RIGHT, 
     $					TMPW2(1,NN-RCNT+1), LDTMP, MYROW,
     $					MYCOL, NPCOL, NPROW, ICTXT)

					CALL DGEMM('N','N',SCNT,NN,NN,ONE,
     $						TMPW2,LDTMP,TMPZ,
     $						LDTMPZ,ZERO,
     $						TMPW,LDTMP)

					CALL DLACPY('A',SCNT,NN-RCNT,TMPW, LDTMP, 
     $					Z(LROW1, LCOL1),LDZ)
					CALL Send(SCNT, RCNT, RIGHT, 
     $					TMPW(1,NN-RCNT+1), LDTMP, MYROW, 
     $					MYCOL, NPCOL, NPROW, ICTXT)
				ELSE IF ((MYROW.EQ.IAROW).AND.
     $				(MYCOL.EQ.IACOL2)) THEN
					LCOL1 = INDXG2L(JC+NN-RCNT, NB, 0, 0, NPCOL)
					Call Send(SCNT, RCNT, LEFT,
     $					Z(LROW1,LCOL1),LDZ, MYROW, MYCOL,
     $					NPCOL, NPROW,ICTXT) 
					Call Recv(SCNT, RCNT, LEFT,
     $					Z(LROW1,LCOL1),LDZ, MYROW, MYCOL,
     $					NPCOL, NPROW,ICTXT)
				END IF
 570			CONTINUE
		END IF
	END IF
	RETURN


*
 800  FORMAT(20F7.2)
      END

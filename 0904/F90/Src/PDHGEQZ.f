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
      INFO = 99
      RETURN
      END

	






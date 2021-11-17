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
      INFO = 99
      RETURN
      END




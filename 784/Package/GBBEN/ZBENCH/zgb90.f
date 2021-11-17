      LOGICAL FUNCTION ZGB90 ( IP, DIM1, DIM2 )
*     .. Scalar Arguments ..
      INTEGER                IP, DIM1, DIM2
*     ..
*
*  Purpose
*  =======
*
*  ZGB90 determines which of two alternative code sections in a GEMM-
*  Based Level 3 BLAS routine that will be the fastest for a particular
*  problem. If the problem is considered large enough ZGB90 returns
*  .TRUE., otherwise .FALSE. is returned. The input parameter IP
*  specifies the calling routine and a break point for alternative code
*  sections. The input parameters DIM1 and DIM2 are matrix dimensions.
*  The returned value is a function of the input parameters and the
*  performance characteristics of the two alternative code sections.
*
*  In this simple implementation, the returned values are determined by
*  looking at only one of the two dimensions DIM1 and DIM2. It may be
*  rewarding to rewrite the logical expressions in ZGB90 so that both
*  dimensions are involved. The returned values should effectively
*  reflect the performance characteristics of the underlying BLAS
*  routines.
*
*
*  Input
*  =====
*
*  IP     - INTEGER
*           On entry, IP specifies which routine and which alternative
*           code sections that the decision is intended for.
*           Unchanged on exit.
*
*  DIM1   - INTEGER.
*           On entry, DIM1 specifies the first dimension in the calling
*           sequence of the Level 3 routine specified by IP.
*           Unchanged on exit.
*
*  DIM2   - INTEGER.
*           On entry, DIM2 specifies the second dimension in the
*           calling sequence of the Level 3 routine specified by IP.
*           Unchanged on exit.
*
*
*  -- Written in May-1994.
*     GEMM-Based Level 3 BLAS.
*     Per Ling, Institute of Information Processing,
*     University of Umea, Sweden.
*
*
*     .. User specified parameters for ZGB90 ..
      INTEGER            ZIP41, ZIP42,
     $                   ZIP51, ZIP52,
     $                   ZIP81, ZIP82, ZIP83,
     $                   ZIP91, ZIP92, ZIP93
      PARAMETER        ( ZIP41 = 4, ZIP42 = 3,
     $                   ZIP51 = 4, ZIP52 = 3,
     $                   ZIP81 = 4, ZIP82 = 3, ZIP83 = 4,
     $                   ZIP91 = 4, ZIP92 = 3, ZIP93 = 4 )
*     ..
*     .. Executable Statements ..
      IF( IP.EQ.41 )THEN
         ZGB90 = DIM1.GE.ZIP41
      ELSE IF( IP.EQ.42 )THEN
         ZGB90 = DIM2.GE.ZIP42
      ELSE IF( IP.EQ.51 )THEN
         ZGB90 = DIM1.GE.ZIP51
      ELSE IF( IP.EQ.52 )THEN
         ZGB90 = DIM2.GE.ZIP52
      ELSE IF( IP.EQ.81 )THEN
         ZGB90 = DIM2.GE.ZIP81
      ELSE IF( IP.EQ.82 )THEN
         ZGB90 = DIM2.GE.ZIP82
      ELSE IF( IP.EQ.83 )THEN
         ZGB90 = DIM1.GE.ZIP83
      ELSE IF( IP.EQ.91 )THEN
         ZGB90 = DIM2.GE.ZIP91
      ELSE IF( IP.EQ.92 )THEN
         ZGB90 = DIM2.GE.ZIP92
      ELSE IF( IP.EQ.93 )THEN
         ZGB90 = DIM1.GE.ZIP93
      ELSE
         ZGB90 = .FALSE.
      END IF
*
      RETURN
*
*     End of ZGB90.
*
      END

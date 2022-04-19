      LOGICAL FUNCTION GETWRD( LN, LLN, IB, IE )
*     .. Scalar Arguments ..
      INTEGER            LLN, IB, IE
*     .. Array Arguments ..
      CHARACTER          LN( LLN )
*
*  Read the first non-blank word from the character string LN. Set
*  the indices IB and IE to the beginning and end of the word,
*  respectively. Return .TRUE. if a word was found and .FALSE. if no
*  word was found.
*
*
*  -- Written in December-1993.
*     GEMM-Based Level 3 BLAS.
*     Per Ling, Institute of Information Processing,
*     University of Umea, Sweden.
*
*  -- Modified in February-1996.
*     Per Ling, Department of Computing Science,
*     Umea University, Sweden.
*
*     ..
*     .. Executable Statements ..
*
*     Find the beginning of the word.
*
      IB = 1
   10 IF( ( LN( IB ).EQ.' ' ).AND.( IB.LT.LLN ) )THEN
         IB = IB + 1
         GO TO 10
      END IF
*
*     Find the end of the word.
*
      IE = IB
   20 IF( IE.LT.LLN )THEN
         IF( LN( IE+1 ).NE.' ' )THEN
            IE = IE + 1
            GO TO 20
         END IF
      END IF
*
*     Check if any word was found.
*
      IF( LN( IB ).NE.' ' )THEN
         GETWRD = .TRUE.
      ELSE
         GETWRD = .FALSE.
      END IF
*
      RETURN
*
*     End of GETWRD.
*
      END

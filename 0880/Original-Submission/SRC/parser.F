FUNCTION PARSER( RECORD ) RESULT( OUTPUT_FROM_PARSER )
!
USE GSTEDEFINITIONS
!
!.. Scalar Argument ..
CHARACTER( * ) :: RECORD
!
!.. Function Result ..
TYPE( DATA_FROM_RECORD ), POINTER :: OUTPUT_FROM_PARSER
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  PARSER parsers one record (line) of the input file.                         !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  RECORD  (input) CHARACTER( * )                                              !
!          One record (line) of the input file.                                !
!                                                                              !
!==============================================================================!
!
!.. Local Scalars ..
CHARACTER( LEN=FIELD_LENGTH ) :: TEMP
CHARACTER( LEN=1 ) :: ONECHAR
INTEGER :: ICHAR, JCHAR
!
!.. Local Derived Data Type ..
TYPE( DATA_FROM_RECORD ), POINTER :: LIST, LIST_CURRENT, LIST_LAST
!
!.. Executable Statements ......................................................
!
ICHAR = 0; JCHAR = 0; TEMP = ''
ALLOCATE( LIST ); NULLIFY( LIST%NEXT ); LIST_CURRENT => LIST
!
DO WHILE ( ICHAR < LEN( RECORD ) )
   ICHAR = ICHAR + 1
   ONECHAR = RECORD( ICHAR:ICHAR )
   IF      ( ONECHAR == '%' ) THEN
           EXIT
   ELSE IF ( ONECHAR /= ' ' .AND. ONECHAR /= ',' ) THEN
           JCHAR = JCHAR + 1 
           TEMP( JCHAR:JCHAR ) = ONECHAR
   ELSE IF ( JCHAR /= 0 ) THEN
           LIST_LAST => LIST_CURRENT
           LIST_CURRENT%FIELD = TEMP; TEMP = ''
           ALLOCATE( LIST_CURRENT%NEXT ); NULLIFY( LIST_CURRENT%NEXT%NEXT )
           LIST_CURRENT => LIST_CURRENT%NEXT; LIST_CURRENT%FIELD = ''
           JCHAR = 0
   END IF
END DO
!
IF ( LIST_CURRENT%FIELD == '' ) NULLIFY( LIST_LAST%NEXT )
ALLOCATE( OUTPUT_FROM_PARSER )
OUTPUT_FROM_PARSER = LIST
!
END FUNCTION PARSER

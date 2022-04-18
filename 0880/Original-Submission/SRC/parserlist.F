FUNCTION PARSERLIST( STRING, LIST ) RESULT( OUTPUT_FROM_PARSERLIST )
!
USE GSTEDEFINITIONS
USE GSTEINTERFACES1, ONLY : GETINTGR, PARSER
!
!.. Scalar Argument ..
CHARACTER( * ) :: STRING
!
!.. Derived Data Type Argument ..
TYPE( DATA_FROM_RECORD ), TARGET :: LIST
!
!.. Function Result ..
TYPE( DATA_FROM_RECORD ), POINTER :: OUTPUT_FROM_PARSERLIST
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  PARSERLIST parsers the fields of list.                                      !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  STRING  (input) CHARACTER                                                   !
!          String to be printed in case of error.                              !
!                                                                              !
!  LIST    (input) DATA_FROM_RECORD                                            !
!          List of data as interpreted by PARSER.                              !
!                                                                              !
!==============================================================================!
!
!.. Local Scalars ..
CHARACTER( LEN=FIELD_LENGTH ) :: TEMP_FIELD
LOGICAL :: COLON
INTEGER :: ICHAR, N, NINC, NMAX, NMIN, NINDEX
!
!.. Static Array ..
INTEGER :: INDEX( 3 )
!
!.. Local Derived Data Type ..
TYPE( DATA_FROM_RECORD ), POINTER :: LIST_IN
TYPE( DATA_FROM_RECORD ), POINTER :: LIST_OUT, LIST_OUT_CURRENT, LIST_OUT_LAST
TYPE( DATA_FROM_RECORD ), POINTER :: TEMP_LIST
!
!.. External Function ..
INTEGER, EXTERNAL :: LISTLENGHT
!
!.. Executable Statements ......................................................
!
LIST_IN => LIST
ALLOCATE( LIST_OUT ); NULLIFY( LIST_OUT%NEXT ); LIST_OUT_CURRENT => LIST_OUT
!
DO WHILE ( ASSOCIATED( LIST_IN ) )
!
!  Check for loop constructions (defined by colons).
!
   INDEX = 0
   ICHAR = 0  
   COLON = .FALSE.
   TEMP_FIELD = LIST_IN%FIELD
   DO WHILE ( ICHAR < FIELD_LENGTH )
      ICHAR = ICHAR + 1
      IF      ( TEMP_FIELD( ICHAR:ICHAR ) == ' ' ) THEN
              EXIT
      ELSE IF ( TEMP_FIELD( ICHAR:ICHAR ) /= ':' ) THEN
              COLON = .FALSE.
      ELSE 
              IF ( COLON ) CALL HANDLER( 1, STRING )
              TEMP_FIELD( ICHAR:ICHAR ) = ' '
              COLON = .TRUE.
      END IF
   END DO
   IF ( TEMP_FIELD( 1:1 ) == ':' .OR. COLON ) CALL HANDLER( 1, STRING )
!
   TEMP_LIST => PARSER( TEMP_FIELD ); NINDEX = LISTLENGHT( TEMP_LIST )
   INDEX( 1:NINDEX ) = GETINTGR( STRING, TEMP_LIST, NINDEX )
!
!  Set loop bounds and increment.
!
   IF      ( NINDEX == 1 ) THEN
           NMIN = INDEX( 1 )
           NMAX = INDEX( 1 )
           NINC = 1
   ELSE IF ( NINDEX == 2 ) THEN
           NMIN = INDEX( 1 )
           NMAX = INDEX( 2 )
           NINC = 1
   ELSE IF ( NINDEX == 3 ) THEN
           NMIN = INDEX( 1 )
           NMAX = INDEX( 3 )
           NINC = INDEX( 2 )
   ELSE 
           CALL HANDLER( 1, STRING )
   END IF
   IF ( NMIN < 1 )    CALL HANDLER( 1, STRING )
   IF ( NINC < 1 )    CALL HANDLER( 1, STRING )
   IF ( NMAX < NMIN ) CALL HANDLER( 1, STRING )
!
!  Create list of dimensions using loop.   
!
   DO N = NMIN, NMAX, NINC
      LIST_OUT_LAST => LIST_OUT_CURRENT
      WRITE ( LIST_OUT_CURRENT%FIELD, FMT=IFORMAT ) N
      ALLOCATE( LIST_OUT_CURRENT%NEXT ); NULLIFY( LIST_OUT_CURRENT%NEXT%NEXT )
      LIST_OUT_CURRENT => LIST_OUT_CURRENT%NEXT; LIST_OUT_CURRENT%FIELD = ''
   END DO
!
   LIST_IN => LIST_IN%NEXT
!
END DO
!
IF ( LIST_OUT_CURRENT%FIELD == '' ) NULLIFY( LIST_OUT_LAST%NEXT )
ALLOCATE( OUTPUT_FROM_PARSERLIST )
OUTPUT_FROM_PARSERLIST = LIST_OUT
!
END FUNCTION PARSERLIST

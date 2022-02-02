! File transformer.f95
! A transformer program from FORTRAN 77 to Fortran 95 have been known 
! [Ref. (3)], but another program was made here. This program name is 
! transformer.  PROGRAM transformer is written in Fortran 95.
!
! PROGRAM transformer modifies the following terms.
! (1) The letters in the columns 73-80.
! (2) The intrinsic functions CDABS, CDSIN, CDCOS, CDEXP, CDLOG, CDSQRT,
!     DFOAT, IMAG and the statements COMPLEX*16, REAL*8, REAL*16.
! (3) Comment lines.
! (4) Continuation lines.
! (5) Unmatched DO constructs.
! (6) Type-statement CHARACTER*( ).
!
! The usage: first compile this file and execute the object program,
! and follow the messages written in the standard output unit.
!
  PROGRAM transformer
  IMPLICIT NONE
  INTEGER, PARAMETER:: i8=8
  INTEGER, PARAMETER:: il=200
  INTEGER:: i,io,i1,i2,i5,j,k,n
  CHARACTER(LEN=30):: charname
  CHARACTER(LEN=80):: char11, char12,char15(il)
  CHARACTER(LEN=5):: char5(i8),char6(i8),char1     ! Statement labels
  CHARACTER(LEN=5), PARAMETER:: char5num='00000'
!  CHARACTER(LEN=2), PARAMETER:: char2='kp'
  DATA char5/i8*' '/, char6/i8*' '/
  PRINT *,'The input file must store the original FORTRAN 77 program.'
  PRINT *,'The input file must be put in the present directory.'
  PRINT *,'The suffix of the input file name must be .f.'
  PRINT *,'Input the input file name without the suffix.'
  READ  *,charname
  OPEN(2,FILE=TRIM(charname)//'.f',ACTION='READ',STATUS='OLD')
  OPEN(3,STATUS='SCRATCH')
  OPEN(7,FILE=TRIM(charname)//'.f95')
  PRINT *,'The file '//TRIM(charname)//'.f95 was generated'
  PRINT *,'The transformed  Fortran 95 program is in the file ',&
      TRIM(charname)//'.f95.'
6 FORMAT(A)
! To delete the letters in the columns 73-80, and to fix comment lines.
  REWIND(2);  REWIND(3)
  DO
     READ(2,6,IOSTAT=io) char11
     IF(io <= -1) EXIT
     char11(73:) = '          '
     IF(SCAN(char11(1:1),'C*') /= 0) THEN
        char11(1:1)='!';   WRITE(3,6) char11;   CYCLE
     ENDIF
     IF(LEN_TRIM(char11)==0) THEN
        WRITE(3,6) char11;  CYCLE
     ENDIF
     DO i=1,72
        IF(IACHAR(char11(i:i))==9) THEN
           char11(i:i)=' ';  char11=char11(:i)//'       '//char11(i+1:)
        ENDIF
     ENDDO
     IF(VERIFY(char11(:5),' 0123456789') /= 0) char11=' '
     IF(VERIFY(char11(6:6),' 0')==1.AND.VERIFY(char11(:5),' ')>=1) char11=' '
     WRITE(3,6) char11
  END DO

! To modify the intrinsic functions CDABS, CDSIN, CDCOS, CDEXP, CDLOG, CDSQRT,
! DFOAT, IMAG, DCMPLX, and the statements COMPLEX*16, REAL*8, REAL*16.
!  COMPLEX FUNCTION XXX*16
  REWIND(3);  REWIND(7)
  DO
    READ(3,6,IOSTAT=io) char11
    IF(io <= -1) EXIT
    IF(char11(:1) == '!' .OR. LEN_TRIM(char11) == 0) THEN
       WRITE(7,6) char11;    CYCLE
    ENDIF
!  To transform CDABS to ABS.
    i=INDEX(char11,'CDABS(')
    IF(i >= 1) char11(i:i+4)='  ABS'
    i=INDEX(char11,'CDABS(')
    if(i >= 1) char11(i:i+4)='  ABS'
!  To transform DABS to ABS.
    i=INDEX(char11,'DABS(')
    IF(i >= 1) char11(i:i+3)=' ABS'
    i=INDEX(char11,'DABS(')
    IF(i >= 1) char11(i:i+3)=' ABS'
!  To transform CDSIN to SIN.
    i=INDEX(char11,'CDSIN(')
    IF(i >= 1) char11(i:i+4)='  SIN'
    i=INDEX(char11,'CDSIN,')
    IF(i >= 1) char11(i:i+5)='      '  ! To erase the declaration
!  To transform DSIN to SIN.
    i=INDEX(char11,'DSIN(')
    IF(i >= 1) char11(i:i+3)=' SIN'
!  To transform CDCOS to COS.
    i=INDEX(char11,'CDCOS(')
    IF(i >= 1) char11(i:i+4)='  COS'
    i=INDEX(char11,'CDCOS,')
    IF(i >= 1) char11(i:i+5)='      '  ! To erase the declaration
!  To transform DCOS to COS.
    i=INDEX(char11,'DCOS(')
    IF(i >= 1) char11(i:i+3)=' COS'
!  To transform DTAN to TAN.
    i=INDEX(char11,'DTAN(')
    IF(i >= 1) char11(i:i+3)=' TAN'
!  To transform CDEXP to EXP.
    i=INDEX(char11,'CDEXP(')
    IF(i >= 1) char11(i:i+4)='  EXP'
    i=INDEX(char11,'CDEXP,')
    IF(i >= 1) char11(i:i+5)='      '  ! To erase the declaration
!  To transform DEXP to EXP.
    i=INDEX(char11,'DEXP(')
    IF(i >= 1) char11(i:i+3)=' EXP'
!  To transform CDLOG to LOG.
    i=INDEX(char11,'CDLOG(')
    IF(i >= 1) char11(i:i+4)='  LOG'
    i=INDEX(char11,'CDLOG,')
    IF(i >= 1) char11(i:i+5)='      '  ! To erase the declaration
!  To transform DLOG to LOG.
    i=INDEX(char11,'DLOG(')
    IF(i >= 1) char11(i:i+3)=' LOG'
!  To transform CDEXP to EXP.
    i=INDEX(char11,'CDEXP(')
    IF(i >= 1) char11(i:i+4)='  EXP'
!  To transform CDSQRT to SQRT.
    i=INDEX(char11,'CDSQRT(')
    IF(i >= 1) char11(i:i+5)='  SQRT'
    i=INDEX(char11,'CDSQRT,')
    IF(i >= 1) char11(i:i+6)='       '
!  To transform DSQRT to SQRT.
    i=INDEX(char11,'DSQRT(')
    IF(i >= 1) char11(i:i+4)=' SQRT'
!  To transform DINT to AINT.
    i=INDEX(char11,'DINT(')
    IF(i >= 1) char11(i:i+3)='AINT'
!  To transform DATAN2 to ATAN2.
    i=INDEX(char11,'DATAN2(')
    IF(i >= 1) char11(i:i+5)=' ATAN2'
!  To transform DMAX1 to MAX.
    i=INDEX(char11,'DMAX1(')
    IF(i >= 1) char11(i:i+4)='  MAX'
!  To transform DMIN1 to MIN.
    i=INDEX(char11,'DMIN1(')
    IF(i >= 1) char11(i:i+4)='  MIN'
!  To transform CDDFOAT to DBLE.
    i=INDEX(char11,'DFLOAT(')
    IF(i >= 1) char11(i:i+5)='  DBLE'
    i=INDEX(char11,'DFLOAT,')
    IF(i >= 1) char11(i:i+6)='       '
!  To transform IMAG to AIMAG.
    i1=1
    DO
    i=INDEX(char11(i1:),'IMAG(')
    i2=i+i1-1
    IF(i >= 1 .and. char11(i2-1:i2-1) /= 'A') THEN
       char11=char11(:i2-1)//'AIMAG('//char11(i2+5:)
       i1=i2+2
    ELSE
       EXIT
    ENDIF
    ENDDO
!  To transform COMPLEX*16 to COMPLEX(//chk8//).
    i=INDEX(char11,'COMPLEX*16')
    IF(i >= 1) char11(i:)='COMPLEX(16)'//char11(i+10:)
!  To transform REAL*8 to REAL(8).
     i=INDEX(char11,' REAL*8')
     IF(i >= 1) THEN
        char11=char11(:i-1)//' REAL(8) '//char11(i+7:)
     ENDIF
!  To transform REAL*16 to REAL(16).
     i=INDEX(char11,' REAL*16')
     IF(i >= 1) THEN
        char11=char11(:i-1)//' REAL(16)'//char11(i+8:)
     ENDIF
     WRITE(7,6) char11
  ENDDO

! To modify continuation lines.
  REWIND(3);  REWIND(7)
  DO
     READ(7,6,IOSTAT=io) char11
     IF(io <= -1) EXIT
     WRITE(3,6) char11
  END DO
  REWIND(3);   REWIND(7)
  char15 = ' '
  READ(3,6,IOSTAT=io) char15(1)
  IF(char15(1)(6:6) == '0') char15(1)(6:6) = ' '
  aa: DO
     DO k=2,il
        READ(3,6,IOSTAT=io) char15(k)
        IF(io <= -1) EXIT aa
        IF(char15(k)(6:6) == '0') char15(k)(6:6) = ' '
        IF(char15(k)(1:1) == '!' .OR. LEN_TRIM(char15(k)) == 0) CYCLE
        IF(char15(k)(6:6) == ' ') EXIT
        char15(k)(6:6) = ' '
        j = LEN_TRIM(char15(1))
        char15(1)(j+1:j+1) = '&'
        j = VERIFY(char15(k),' ')
        char15(k)(j-1:j-1) = '&'
        EXIT
     ENDDO
     DO n=1,k-1
        WRITE(7,6) char15(n)
     ENDDO
     char15(1)=char15(k)
  END DO aa
  WRITE(7,6) char15(1)

! To modify unmatched DO constructs.
  REWIND(3);  REWIND(7)
  DO
     READ(7,6,IOSTAT=io) char11
     IF(io <= -1) EXIT
     WRITE(3,6) char11
  END DO
  REWIND(3);   REWIND(7)
  i5=0
  DO
     READ(3,6,IOSTAT=io) char11
     IF(io <= -1) EXIT
     IF(char11(1:1) == '!' .OR. LEN_TRIM(char11) == 0) THEN
        WRITE(7,6) char11
        CYCLE
     ENDIF
     char12 = char11(7:)
     char12 = ADJUSTL(char12)
     IF(char12(:3) == 'DO ') THEN
     char12(:2) = '  '
     char12 = ADJUSTL(char12)
     i = VERIFY(char12,'0123456789')
        IF(i >= 2) THEN
           char5(2:)=char5(:i8-1);  char6(2:)=char6(:i8-1)
           i5=i5+1
           char5(1) = char12(:i-1);  char6(1)=char5(1)
           IF(i5 >= 2 .AND. char5(1)==char5(2)) THEN
              char6(1) = char12(:i-1)//char5num(i:)
              WRITE(char6(1)(5:5),'(I1)') i5-2
              i1=VERIFY(char11(7:),' ')
              char11(6+i1:)='DO '//char6(1)//' '//char12(i+1:)
           ENDIF
           WRITE(7,6) char11;  CYCLE
        END IF
     END IF
     IF(i5 == 0 .OR. INDEX(char11(:5), TRIM(char5(1))) == 0) THEN
        WRITE(7,6) char11;   CYCLE
     END IF
     IF(INDEX(char11(:5), TRIM(char5(1))) >= 1) THEN
        DO i=1,i8
           IF(char5(i) /= char5(1)) EXIT
           IF(char6(i) == char5(1)) char6(i)=char11(:5)
        ENDDO
     ENDIF
     IF(ADJUSTL(char11(7:)) /= 'CONTINUE') THEN
        char11(:5)='     ';  WRITE(7,6) char11
     ENDIF
     DO
        char11(:5)=char6(1)
        i1=VERIFY(char11(7:),' ')
        char11(6+i1:)='CONTINUE'
        char1=char5(1)
        char5(:i8-1)=char5(2:);  char6(:i8-1)=char6(2:)
        i5=i5-1
        WRITE(7,6) char11
        IF(char1 /= char5(1)) EXIT
     ENDDO
  ENDDO

! To transform CHARACTER* to CHARACTER(LEN=  ), and
! to transform CHARACTER*(*) to CHARACTER(LEN=*).
  REWIND(3);   REWIND(7)
  DO
     READ(7,6,IOSTAT=io) char11
     IF(io <= -1) EXIT
     WRITE(3,6) char11
  END DO
  REWIND(3);   REWIND(7)
  DO
     READ(3,6,IOSTAT=io) char11
     IF(io <= -1) EXIT
     IF(char11(:1) == '!' .OR. LEN_TRIM(char11) == 0) THEN
        WRITE(7,6) TRIM(char11);   CYCLE
     ENDIF
     i=INDEX(char11,'CHARACTER*')
     IF(i >= 1) THEN
       char12=char11(i+10:)
       i1=verify(char12,' 0123456789')
       if(i1 >= 2) then
       char11=char11(:i-1)//'CHARACTER(LEN='//char12(:i1-1)//') '&
          //char11(i+8+i1:)
       ENDIF
     ENDIF
     i=INDEX(char11,'CHARACTER*(')
     IF(i >= 1) THEN
        char11 = char11(:i+8)//'(LEN='//char11(i+11:)
     ENDIF
     WRITE(7,6) TRIM(char11)
  ENDDO
  END PROGRAM transformer



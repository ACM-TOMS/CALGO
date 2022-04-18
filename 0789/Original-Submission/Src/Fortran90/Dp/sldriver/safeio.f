C***+****|****+****|* COPYRIGHT J D PRYCE 1998 **|****+****|****+****|**
      module SAFEIO
      logical:: PLAYBK=.FALSE.
      integer, parameter:: IPBK=19

      contains

C***+****|****+****|**SAFEIO package***|****+****|****+****|****+****|**
C Simple 'secure, abortable interactive input' routines
C
C The following are logical functions which return values through
C the argument list and also
C either
C   .TRUE. through the function name meaning 'successful input'
C   and the input value(s) through its (first) argument
C or
C    .FALSE. through the function name meaning 'aborted by user'
C The routines are
C GETI(X)  - get one Integer value X
C GETR(X)  - get one Real value X
C GETIR(X,XLO,XHI)  - get one Integer value X in a Range
C GETRR(X,XLO,XHI)  - get one Real value X in a Range
C GETIS(X,NX) - get NX integer values into array X
C GETRS(X,NX) - get NX real values into array X
C
C The user is to signal 'abort' by typing z or Z (followed by Enter)
C
C The following is a non-abortable routine requiring a yes/no answer
C YESNO() - character*1 function, returns either 'y' or 'n'
C The following non-abortable routine waits for ENTER to be pressed
C SPAUSE()
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine OPNPBK
      integer IOERR
      open(IPBK,file='playback.dat', status='REPLACE', action='WRITE',
     +  iostat=IOERR)
      if (IOERR.eq.0) then
        PLAYBK = .TRUE.
      else
        write(*,*) 'Unable to open playback-file unit',IPBK,
     +   ', file playback.out, IOSTAT=',IOERR
        write(*,*) 'It may be in use by another process?'
        PLAYBK = .FALSE.
        call SPAUSE
      end if
      end subroutine OPNPBK

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine CLOPBK
      close(IPBK)
      end subroutine CLOPBK

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function GETI(X)
      implicit none
C     .. Scalar Arguments ..
      integer X
      character*80 LINE
C
   10 read (*,FMT='(a)',ERR=30) LINE
      if (PLAYBK) write(IPBK,'(a)') TRIM(LINE)
      if (LINE(1:1)=='z' .or. LINE(1:1)=='Z') then
        GETI = .FALSE.
      else
        read (LINE,FMT=*,ERR=30) X
        GETI = .TRUE.
      end if
      return

   30 write (*,ADVANCE='NO',FMT=
     +  '(1x,''***Integer value required, please retype: '')')
      go to 10
      end function GETI

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function GETR(X)
      implicit none
C     .. Scalar Arguments ..
      double precision X
      character*80 LINE
C     ..
   10 read (*,FMT='(a)',ERR=30) LINE
      if (PLAYBK) write(IPBK,'(a)') TRIM(LINE)
      if (LINE(1:1)=='z' .or. LINE(1:1)=='Z') then
        GETR = .FALSE.
      else
        read (LINE,FMT=*,ERR=30) X
        GETR = .TRUE.
      end if
      return
   30 write (*,ADVANCE='NO',FMT=
     +  '(1x,''***Real value required, please retype: '')')
      go to 10
      end function GETR
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function GETIR(X,XLO,XHI)
      implicit none
C     .. Scalar Arguments ..
      integer X,XHI,XLO
C     ..
      character*80 LINE
C     ..
   10 read (*,FMT='(a)',ERR=30) LINE
      if (PLAYBK) write(IPBK,'(a)') TRIM(LINE)
      if (LINE(1:1)=='z' .or. LINE(1:1)=='Z') then
        GETIR = .FALSE.
      else
        read (LINE,FMT=*,ERR=30) X
        if (X.lt.XLO .or. X.gt.XHI) go to 40
        GETIR = .TRUE.
      end if
      return
   30 write (*,ADVANCE='NO',FMT=
     +  '(1x,''***Integer value required, please retype: '')')
      go to 10
   40 write (*,FMT=45) XLO,XHI
   45 format(1x,'***Not in range',i9,' :',i9,', please retype: ')
      go to 10
      end function GETIR

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function GETRR(X,XLO,XHI)
      implicit none
C     .. Scalar Arguments ..
      double precision X,XHI,XLO
C     ..
      character*80 LINE
C     ..
   10 read (*,FMT='(a)',ERR=30) LINE
      if (PLAYBK) write(IPBK,'(a)') TRIM(LINE)
      if (LINE(1:1)=='z' .or. LINE(1:1)=='Z') then
        GETRR = .FALSE.
      else
        read (LINE,FMT=*,ERR=30) X
        if (X.lt.XLO .or. X.gt.XHI) go to 40
        GETRR = .TRUE.
      end if
      return
   30 write (*,ADVANCE='NO',FMT=
     +  '(1x,''***Integer value required, please retype: '')')
      go to 10
   40 write (*,FMT=45) XLO,XHI
   45 format(1x,'***Not in range',i9,' :',i9,', please retype: ')
      go to 10
      end function GETRR

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function GETIS(X,NX)
      implicit none
C     .. Scalar Arguments ..
      integer NX
C     ..
C     .. Array Arguments ..
      integer X(NX)
C     ..
C     .. Local Scalars ..
      integer I
      character*80 LINE
C     ..
      I = 1
   10 read (*,FMT='(a)',ERR=30) LINE
      if (PLAYBK) write(IPBK,'(a)') TRIM(LINE)
      if (LINE(1:1)=='z' .or. LINE(1:1)=='Z') then
        GETIS = .FALSE.
      else
        read (LINE,FMT=*,ERR=30) (X(I),I=I,NX)
        GETIS = .TRUE.
      end if
      return

   30 write (*,ADVANCE='NO',FMT=35) I
   35 format
     +  ('***Integer values required, please retype from item',i2,': ')
      go to 10
      end function GETIS

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function GETRS(X,NX)
      implicit none
C     .. Scalar Arguments ..
      integer NX
C     ..
C     .. Array Arguments ..
      double precision X(NX)
C     ..
C     .. Local Scalars ..
      integer I
      character*80 LINE
C     ..
      I = 1
   10 read (*,FMT='(a)',ERR=30) LINE
      if (PLAYBK) write(IPBK,'(a)') TRIM(LINE)
      if (LINE(1:1)=='z' .or. LINE(1:1)=='Z') then
        GETRS = .FALSE.
      else
        read (LINE,FMT=*,ERR=30) (X(I),I=I,NX)
        GETRS = .TRUE.
      end if
      return

   30 write (*,ADVANCE='NO',FMT=35) I
   35 format
     +  ('***Real values required, please retype from item',i2,': ')
      go to 10
      end function GETRS

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      character*1 function YESNO(ASK)
      implicit none
      character*(*) ASK
C     .. Local Scalars ..
      character YN
C     ..
      write(*,'(a)',advance='NO') ASK
   10 read (*,FMT='(a)',END=40,ERR=30) YN
      if (PLAYBK) write(IPBK,'(a)') YN
      if (YN.eq.'y' .or. YN.eq.'Y') then
         YESNO = 'y'
      else if (YN.eq.'n' .or. YN.eq.'N') then
         YESNO = 'n'
      else
         go to 30
      end if
      return

   30 write (*,ADVANCE='NO',
     +  FMT='(1x,''***Please type one of y,Y,n,N: '')')
      go to 10
   40 write (*,*) 'Don''t use end-of-file,',
     + ' it means end of run in most F90 & some F77 systems'
      write (*,*) 'Sorry!'
      stop
      end function YESNO

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine SPAUSE

      write (*,ADVANCE='NO',FMT='(a)') 'Press ENTER to continue'
      read (*,FMT=*,END=100,ERR=100)
  100 return

      end subroutine SPAUSE
C***+****|****+****|****end of SAFEIO package****|****+****|****+****|**

      end module SAFEIO

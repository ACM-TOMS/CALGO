C***+****|****+****|* COPYRIGHT J D PRYCE 1998 **|****+****|****+****|**
C BATCHIO is a module with the same name, SAFEIO, and interface as the
C interactive SAFEIO package, but intended for use when running a
C program non-interactively. Differences:
C 1. SPAUSE doesn't wait for ENTER to be pressed but waits 1 second
C 2. 'ERR=' on input leads to a STOP.
C 3. They echo their input to the screen.
C The result is that the screen output in batch mode (using < to
C read data from a file) should be identical to that produced by
C running interactively, if the data has no errors!

      module SAFEIO
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
C     Dummy routine
      end subroutine OPNPBK

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine CLOPBK
C     Dummy routine
      end subroutine CLOPBK

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**

      logical function GETI(X)
      implicit none
C     .. Scalar Arguments ..
      integer X
      character*80 LINE
C     ..
   10 read (*,FMT='(a)') LINE
      call PRTRIM(LINE)
      if (LINE(1:1)=='z' .or. LINE(1:1)=='Z') then
        GETI = .FALSE.
        return
      else
        read (LINE,FMT=*,ERR=30) X
        GETI = .TRUE.
        return
   30   write (*,ADVANCE='NO',FMT=
     +  '(1x,''***Integer value required, please retype: '')')
        stop
      end if
      end function GETI

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function GETR(X)
      implicit none
C     .. Scalar Arguments ..
      double precision X
      character*80 LINE
C     ..
   10 read (*,FMT='(a)') LINE
      call PRTRIM(LINE)
      if (LINE(1:1)=='z' .or. LINE(1:1)=='Z') then
        GETR = .FALSE.
        return
      else
        read (LINE,FMT=*,ERR=30) X
        GETR = .TRUE.
        return
   30   write (*,ADVANCE='NO',FMT=
     +  '(1x,''***Real value required, please retype: '')')
        stop
      end if
      end function GETR
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function GETIR(X,XLO,XHI)
      implicit none
C     .. Scalar Arguments ..
      integer X,XHI,XLO
C     ..
      character*80 LINE
C     ..
   10 read (*,FMT='(a)') LINE
      call PRTRIM(LINE)
      if (LINE(1:1)=='z' .or. LINE(1:1)=='Z') then
        GETIR = .FALSE.
        return
      else
        read (LINE,FMT=*,ERR=30) X
        if (X.lt.XLO .or. X.gt.XHI) go to 40
        GETIR = .TRUE.
        return
   30   write (*,ADVANCE='NO',FMT=
     +  '(1x,''***Integer value required, please retype: '')')
        stop
   40   write (*,FMT=45) XLO,XHI
   45   format(1x,'***Not in range',i9,' :',i9,', please retype: ')
        stop
      end if
      end function GETIR

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function GETRR(X,XLO,XHI)
      implicit none
C     .. Scalar Arguments ..
      double precision X,XHI,XLO
C     ..
      character*80 LINE
C     ..
   10 read (*,FMT='(a)') LINE
      call PRTRIM(LINE)
      if (LINE(1:1)=='z' .or. LINE(1:1)=='Z') then
        GETRR = .FALSE.
        return
      else
        read (LINE,FMT=*,ERR=30) X
        if (X.lt.XLO .or. X.gt.XHI) go to 40
        GETRR = .TRUE.
        return
   30   write (*,ADVANCE='NO',FMT=
     +  '(1x,''***Integer value required, please retype: '')')
        stop
   40   write (*,FMT=45) XLO,XHI
   45   format(1x,'***Not in range',i9,' :',i9,', please retype: ')
        stop
      end if
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
   10 read (*,FMT='(a)') LINE
      call PRTRIM(LINE)
      if (LINE(1:1)=='z' .or. LINE(1:1)=='Z') then
        GETIS = .FALSE.
        return
      else
        read (LINE,FMT=*,ERR=30) (X(I),I=I,NX)
        GETIS = .TRUE.
        return

   30   write (*,ADVANCE='NO',FMT=35) I
   35   format
     +  ('***Integer values required, please retype from item',i2,': ')
        stop
      end if
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
   10 read (*,FMT='(a)') LINE
      call PRTRIM(LINE)
      if (LINE(1:1)=='z' .or. LINE(1:1)=='Z') then
        GETRS = .FALSE.
        return
      else
        read (LINE,FMT=*,ERR=30) (X(I),I=I,NX)
        GETRS = .TRUE.
        return

   30   write (*,ADVANCE='NO',FMT=35) I
   35   format
     +  ('***Real values required, please retype from item',i2,': ')
        stop
      end if
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
      write (*,'(a)') YN
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
      stop
   40 write (*,*) 'Don''t use end-of-file,',
     + ' it means end of run in most F90 & some F77 systems'
      write (*,*) 'Sorry!'
      stop
      end function YESNO

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine SPAUSE
      double precision DELAY
      integer CTRL_C
      parameter (DELAY=1d0, CTRL_C=3)
      integer ICURR,IRATE,ISTART
      integer(kind=2) K
      write (*,FMT='(a)') 'Press ENTER to continue'
C      read (*,FMT=*,END=100,ERR=100)
      call SYSTEM_CLOCK(ISTART,IRATE)
   50 continue
!C start of Salford FTN90-specific code
!C     Salford FTN90 routine to check if key has been pressed
!      call GET_KEY1@(K)
!      if (K.eq.CTRL_C) then
!        print*,'User-requested Break'
!        stop
!      end if
!C end of Salford FTN90-specific code
      call SYSTEM_CLOCK(ICURR,IRATE)
      if (dble(ICURR-ISTART)/dble(IRATE)<DELAY) go to 50
  100 return

      end subroutine SPAUSE
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine PRTRIM(LINE)
C This is superfluous in F90 but put in for F77 compatibility
      character*(*) LINE
      integer I
      I=LEN(LINE)+1
  10  I = I-1
      if (I.gt.1 .and. LINE(I:I).eq.' ') go to 10
      write (*,'(a)') LINE(1:I)
      end subroutine PRTRIM
C***+****|****+****|****end of SAFEIO package****|****+****|****+****|**

      end module SAFEIO

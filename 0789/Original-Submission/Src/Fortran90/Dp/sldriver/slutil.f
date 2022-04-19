C***+****|****+****|* COPYRIGHT J D PRYCE 1998 **|****+****|****+****|**
      module SLUTIL
      implicit none
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      contains
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine MKMESH(N,A,B,AREG,BREG,X)
      integer N
      double precision A,B,X(1:N)
      logical AREG,BREG
C Makes a mesh of N points in interval [A,B]. Points are more widely
C spaced for large (negative or positive) X in [A,B], a `natural scale'
C of O(10) is assumed so for small |A|,|B| the mesh is nearly equally
C spaced. If AREG (resp. BREG) is .FALSE. then A (resp. B)
C is `singular' and the mesh is formed by increasing N by 1 or 2 as
C needed, forming a mesh that includes A and/or B, and then discarding
C the singular points.
      double precision SCALE,XINFTY
      parameter (SCALE=10D0, XINFTY=1D35)
      integer I,IOFSET,NN
      double precision H,T,TA,TB
      NN = N-1
      IOFSET = 1
      if (.not. AREG) then
        NN = NN+1
        IOFSET = 0
      end if
      if (.not.BREG) NN = NN+1
      TA = A/(SCALE+abs(A))
      TB = B/(SCALE+abs(B))
      H = (TB-TA)/NN
      do 100 I=1,N
        T = TA + (I-IOFSET)*H
        if (T.ge.1D0) then
          X(I) = XINFTY
        else if (T.le.-1D0) then
          X(I) = -XINFTY
        else
          X(I) = SCALE*T/(1D0-abs(T))
        end if
  100 continue
C Ensure the endpoints if included are *exact*
      if (AREG) X(1) = A
      if (BREG) X(N) = B
c      print*,'exit MKMESH: NMESH(=N),MESH(=X):',N,X
      end subroutine MKMESH

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine SELMSH(NIN,XIN,MAXOUT,NOUT,XOUT)
C Selects at most MAXOUT elements from a mesh XIN, spreading the
C selection as evenly as possible over the indices of XIN and putting
C the selected elements in XOUT and the number of them in NOUT.
      integer NIN,MAXOUT,NOUT
      double precision XIN(NIN),XOUT(MAXOUT)
      integer I
      double precision C
      NOUT = min(NIN,MAXOUT)
      C = DBLE(NIN-1)/DBLE(NOUT-1)
      do 10 I=1,NOUT
        XOUT(I) = XIN(NINT(C*(I-1))+1)
c        write(*,'(i5)',advance='NO') NINT(C*(I-1))+1
   10 continue
c      write(*,*)
      end subroutine SELMSH

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      integer function LG2INT(Q)

C     .. Scalar Arguments ..
      logical Q
C     ..
      if (Q) then
         LG2INT = 1
      else
         LG2INT = 0
      end if

      end function LG2INT

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      integer function ITRIM(S)
C Returns position of last nonblank character of S
C     .. Scalar Arguments ..
      character*(*) S
C     ..
C     .. Local Scalars ..
      integer I
C     ..
C     .. Intrinsic Functions ..
      intrinsic LEN
C     ..
      I = LEN(S)
   10 if (I.eq.0) go to 20
      if (S(I:I).ne.' ') go to 20
      I = I - 1
      go to 10

   20 ITRIM = I
      end function ITRIM

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine STTRIM(STR,ISTR)
      character*(*) STR
      integer ISTR
C .. Local scalars ..
      integer I,J
C Does roughly the job of the F90 ITRIM function, in F77 terms.
C If I,J are the positions of 1st and last nonblank characters of STR
C it moves STR(I:J) to STR(1:) and sets ISTR=J-I+1
C it does NOT remove embedded blanks!
      J=LEN(STR)
      if (J.eq.0 .or. STR.eq.' ') then
        ISTR = 0
      else
C Now STR is nonblank so:
  10    if (STR(J:J).eq.' ') then
          J = J-1
          goto 10
        end if
        I = 1
  20    if (STR(I:I).eq.' ') then
          I = I+1
          goto 20
        end if
        STR = STR(I:J)
        ISTR = J-I+1
      end if
      end subroutine STTRIM

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine INT2ST(K,KSTR,IK)
      integer K,IK
      character*12 KSTR
C converts K to decimal character representation held in KSTR(1:IK)
      integer LEN
      parameter(LEN=12)
      integer I

      write(KSTR,'(i12)') K
      I=LEN
 100  I=I-1
      if (KSTR(I:I).ne.' ') goto 100
      IK=LEN-I
      KSTR=KSTR(I+1:LEN)
      end subroutine INT2ST

      end module SLUTIL

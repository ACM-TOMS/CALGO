C***+****|****+****|* COPYRIGHT J D PRYCE 1997 **|****+****|****+****|**
      module TESTMOD
      use SLTSTVAR
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
!Warning by JDP
! The ABINFO routine used by SLEIGN, SLEIGN2 may well be out of date
! since I kept trying new problems frequently
! GTHUHV is currently a dummy routine
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      contains
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine TSTSET(ISWTCH)

C     .. Parameters ..
C+------------------------------------------------------------------+
C| TSETNM is name given to Test Set by its implementer.             |
C| NPROBP is total no. of Problems.                                 |
C| There must be a Problem number I for every I from 1 to NPROBP.   |
C| *NOTE* When altering the Problem Set, re-set TSETNM & NPROBP !!! |
C+------------------------------------------------------------------+
      implicit none
      integer NPROBP
      character*8 TSETNM
      parameter(TSETNM='sample  ', NPROBP=10)
C     These are used in various problems. Some are given a value by the
C     TSTSET(1) call, which is held for later use, so they need to be
C     SAVEd.  Others are used as temporaries so don't need to be SAVEd:
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C     .. Scalar Arguments ..
      integer ISWTCH
C     .. Local variables ..
      double precision TMP,NU,C1,C2
C     .. Save statement ..
      save NU,C1,C2

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
cc      print*,'Calling SAMPLE''s TSTSET with ISWTCH=',ISWTCH
C Put ISWTCH.EQ.2 case first as this is the most frequently called:
      if (ISWTCH.eq.2) then
C     'COEFFN' call:
         go to (12,22,32,42,52,62,72,82,92,102) IPROB

      else if (ISWTCH.eq.-1) then
C     'TSTINI' call:
C        Copy title & total no. of Problems into COMMON variables:
         TITLE(1:8) = TSETNM
         NPROB = NPROBP
         return

      else if (ISWTCH.eq.0) then
C     'SETUP0' call:
C     Start of new problem so set default values of IPARM, NEPRM:
         IPARM = 0
         NEPRM = 0
         go to (10,20,30,40,50,60,70,80,90,100) IPROB

      else if (ISWTCH.eq.1) then
C     'SETUP1' call:
         go to (11,21,31,41,51,61,71,81,91,101) IPROB

      else if (ISWTCH.eq.3) then
         go to (13,23,33,43,53,63,73,83,93,103) IPROB

      else if (ISWTCH.eq.4) then
         go to (14,24,34,44,54,64,74,84,94,104) IPROB

      else
C       This signals a serious coding error in the Package so:
         print *,'Invalid ISWTCH',ISWTCH,' and/or IPROB',IPROB,
     +     ' input to TSTSET call'
         stop

      end if

C***+****|****+****|     START OF PROBLEM SET    |****+****|****+****|**
C Problem #1
   10 TITLE =  'Simple equation -u" = lambda u'
      NPARM = 0
      NEPRM = 0
      PARNM = ' '
      go to 1000

   11 A = 0D0
      B = PI
      ATYPE = 'R'
      BTYPE = 'R'
      SYM = .TRUE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

   12 P = 1D0
      Q = 0D0
      W = 1D0
      go to 1002

   13 U = A2
      PDU = A1
      go to 1003

   14 U = B2
      PDU = B1
      go to 1004

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #2
   20 TITLE = "transformed -u''=lambda u"
      NPARM = 1
      NEPRM = 0
      PARNM = 'NU, must be >0'
      go to 1000

   21 NU = PARM(1)
      C1 = (NU*NU-1D0)/4D0
      C2 = NU+NU-2D0
      A = 0D0
      B = PI**(1D0/NU)
      ATYPE = 'LCN'
      BTYPE = 'R'

      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

   22 P = 1D0
      Q = C1/(X*X)
      W = NU*NU*X**C2
      go to 1002

   23 U = A2
      PDU = A1
      go to 1003

   24 U = B2
      PDU = B1
      go to 1004

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #3
   30 TITLE = "Branko Curgus2: -y''=lambda x y  on [-1,1]"
      NPARM = 0
      NEPRM = 0
      PARNM = ' '
      go to 1000

   31 A = -1D0
      B = 1D0
      ATYPE = 'R'
      BTYPE = 'R'
      SYM = .TRUE.
      A2 = 1D0
      A1 = 01D0
      B2 = 1D0
      B1 = 01D0
      go to 1001

   32 P = 1D0
      Q = 0D0
      W = X
      go to 1002

   33 U = A2
      PDU = A1
      go to 1003

   34 U = B2
      PDU = B1
      go to 1004

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #4
   40 TITLE = "-u''=lambda u with a lambda-dependent BC"
      NPARM = 2
      NEPRM = 0
      PARNM = "A1', A2' where (A1+lam*A1')u(a)=(A2+lam*A2')pu'(a)"
      go to 1000

   41 A = 0d0
      B = PI
      ATYPE = 'R'
      BTYPE = 'R'
      SYM = .FALSE.
      A2 = 0d0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

   42 P = 1D0
      Q = 0D0
      W = 1D0
      go to 1002

   43 U = A2+EIG*PARM(2)
      PDU = A1+EIG*PARM(1)
      go to 1003

   44 U = B2
      PDU = B1
      go to 1004

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #5
   50 TITLE = "Branko Curgus4: -(sign(x)y')'=lambda sign(x) y on [-1,1]"
      NPARM = 0
      NEPRM = 0
      PARNM = ' '
      go to 1000

   51 A = -1d0
      B = 1d0
      ATYPE = 'R'
      BTYPE = 'R'
      SYM = .TRUE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

   52 P = sign(1D0,X)
      Q = 0D0
      W = P
      go to 1002

   53 U = A2
      PDU = A1
      go to 1003

   54 U = B2
      PDU = B1
      go to 1004

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #6
   60 TITLE =  "Klaus 1: -(y'/2)'+(x^2/8(x<0),x^2/2(x>=0)y = lambda y"
      NPARM = 0
      NEPRM = 0
      PARNM = ' '
      go to 1000

   61 A = -XINFTY
      B = XINFTY
      ATYPE = 'LPN'
      BTYPE = 'LPN'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

   62 P = 0.5D0
      if (X.lt.0d0) then
        Q = 0.125d0*X*X
      else
        Q = 0.5d0*X*X
      end if
      W = 1D0
      go to 1002

   63 U = A2
      PDU = A1
      go to 1003

   64 U = B2
      PDU = B1
      go to 1004

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #7
   70 TITLE =  'Simple equation -u" = lambda u with u = A1*x+A2 near 0'
      NPARM = 0
      NEPRM = 0
      PARNM = ' '
      go to 1000

   71 A = 0D0
      B = PI
      ATYPE = 'R'
      BTYPE = 'R'
      SYM = .TRUE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

   72 P = 1D0
      Q = 0D0
      W = 1D0
      go to 1002

   73 U = A1*X+A2
      PDU = A1
      go to 1003

   74 U = B2
      PDU = B1
      go to 1004

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #8
   80 TITLE =  'Another transformed version of Hinz problem'
      NPARM = 0
      NEPRM = 0
      PARNM = ' '
      go to 1000

   81 A = 0D0
      B = 1D0
      ATYPE = 'LCN'
      BTYPE = 'LPNO'
      SYM = .FALSE.
      A2 = 1D0
      A1 = 0D0
      B2 = 1D0
      B1 = 0D0
      go to 1001

   82 TMP = 1D0-X
      P = X*TMP
      Q = X/(TMP**3)*COS(X/TMP)
      W = X/(TMP**3)
      go to 1002

   83 U = A2
      PDU = A1
      go to 1003

   84 U = B2
      PDU = B1
      go to 1004

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #9
   90 TITLE =
     +'Transformed Hinz prob, -(xu'')'' + x cos(x) u = lambda x u'
      NPARM = 0
      NEPRM = 0
      PARNM = ' '
      go to 1000

   91 A = 0D0
      B = XINFTY
      ATYPE = 'LCN'
      BTYPE = 'LPNO'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

   92 P = X
      Q = X*COS(X)
      W = X
      go to 1002

   93 U = A2
      PDU = A1
      go to 1003

   94 U = B2
      PDU = B1
      go to 1004

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Problem #10
  100 TITLE =  'Hinz problem -u" + (cos(x) - 1/(4x^2))u = lambda u'
      NPARM = 0
      NEPRM = 0
      PARNM = ' '
      go to 1000

  101 A = 0D0
      B = XINFTY
      ATYPE = 'LCN'
      BTYPE = 'LPNO'
      SYM = .FALSE.
      A2 = 0D0
      A1 = 1D0
      B2 = 0D0
      B1 = 1D0
      go to 1001

  102 P = 1D0
      Q = COS(X) - 0.25D0/(X*X)
      W = 1D0
      go to 1002

  103 U = A2
      PDU = A1
      go to 1003

  104 U = B2
      PDU = B1
      go to 1004

C***+****|****+****|      END OF PROBLEM SET     |****+****|****+****|**

 1000 continue
 1001 continue
 1002 continue
 1003 continue
 1004 continue
      end subroutine TSTSET

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine ABINFO(P0ATA,QFATA,P0ATB,QFATB)
C This is a special routine for SLEIGN, SLEIGN2
C******* THIS ROUTINE APPLIES TO A SPECIFIC TEST SET ONLY ************
C User beware. It may be out of date.
      implicit none
      integer NPROB
      parameter (NPROB=10)
C     .. Scalar Arguments ..
      double precision P0ATA,QFATA,P0ATB,QFATB
C     .. Local Arrays ..
      integer P0A(1:NPROB),QFA(1:NPROB),P0B(1:NPROB),QFB(1:NPROB)

C      1  2  3  4  5  6  7  8  9 10
      data P0A/0,0,0,0,0,0,0,0,0,0/
      data QFA/1,0,1,1,1,1,1,1,1,0/
      data P0B/0,0,0,0,0,0,0,0,0,1/
      data QFB/1,1,1,1,1,1,1,1,1,1/

      P0ATA = 2*P0A(IPROB) - 1
      QFATA = 2*QFA(IPROB) - 1
      P0ATB = 2*P0B(IPROB) - 1
      QFATB = 2*QFB(IPROB) - 1

C Now special cases need to be handled!
C In some problems this info may depend on values of parameters
C ...

      end subroutine ABINFO

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine GTHUHV(IEND,XA,UA,VA,HU,HV)
C This routine is used only by SOLVRS in its interface to SLEIGN2.
C Its functionality should be incorporated into TSTSET eventually.
C******* THIS ROUTINE APPLIES TO THE SAMPLE TEST SET ONLY ************
C When changing the test set, you must put code in here if you wish
C the limit-circle facility of SLEIGN2 to work correctly!

      implicit none
      integer IEND
      double precision XA,UA,VA,HU,HV
      double precision AA,BB,ALPHA,NU,OM,OMSQ,TEMP

C Returns the HU,HV info needed by SLEIGN2 as part of its UV routine.
C (May be included in main TSTSET eventually)
C Arguments:
C IEND: (input) 0 or 1 for left or right endpoint
C XA:   (input) point where boundary condition info is evaluated
C For the others, see SLEIGN2's documentation.

C     Kluge, to set valid fl.pt. values: SLEIGN2 does funny things..
      HU = 0D0
      HV = 0D0

C Change the labels below to execute appropriate code for any problem
C with a LC endpoint.
      if (ATYPE(1:2).eq.'LC' .or. BTYPE(1:2).eq.'LC') then
        go to(10,10,10,10,10,10,10,10,10,10) IPROB
      else
C       If neither endpoint is LC, UV is not needed so do nothing
        return
      end if

   10 write(*,*) 'Error in GTHUHV, Problem',IPROB,' has LC endpoint',
     +  ' but no HU,HV info implemented for it'
      stop

C Insert code on these lines:
C  120 if (IEND.eq.0) then
C        HU = ...
C        HV = ...
C      else ...
C      end if
C      go to 1000
C
 1000 end subroutine GTHUHV

      end module TESTMOD

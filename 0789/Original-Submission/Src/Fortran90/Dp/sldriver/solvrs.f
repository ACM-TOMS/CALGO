C***+****|****+****|* COPYRIGHT J D PRYCE 1998 **|****+****|****+****|**
C This file contains
C - main module SOLVRS
C - small modules
C     PQWCOM, SL02CM
C   which act as COMMONs for inter-routine communication.
C - auxiliary routines
C     used by SLEDGE:         COEFF
C     used by SLEIGN,SLEIGN2: P,Q,W,R
C     used by SLEIGN2:        UV
C     used by SL02F:          SETUP
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      module PQWCOM
C PQWCOM acts as a COMMON area. It is used by SSLEIG, SSLEI2 because of
C FLIPQ which copes with (q(x) in SLEIGN) being -(q(x) in SLEIGN2 &
C elsewhere)
      double precision XX,PP,QQ,WW
      logical FLIPQ
      end module PQWCOM

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C SL02CM is used to communicate ISING data to SL02F's SETUP routine
      module SL02CM
      character*1 AINFO,BINFO
      end module SL02CM

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      module SOLVRS
      use SLCONSTS
      use SLTSTPAK,only:GETBCS,COEFFN,CPU
      use SAFEIO,only:SPAUSE,YESNO
      implicit none
C***+****|****+****|****+*** Global Data ***+****|****+****|****+****|**
      integer,parameter:: NSOLVR=4
     +                   ,ISLED=1,ISLEIG=2,ISL02=3,ISLEI2=4
     +                   ,OKCALC=0, DEPREC=1, FORBID=2
      character*6,parameter:: SLVNAM(NSOLVR)=
     +  (/'sledge','sleign','sl02f ','sleig2'/)
      integer,parameter:: OKEXTS(NSOLVR)=(/0,1,0,1/)

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      contains
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine SOLVIT(ISOLVR,QCALC,A,B,A1,A2,B1,B2,ATYPE,BTYPE,
     +           KLO,KHI,TOLER,
     +           NXMESH,XMESH,NLMESH,LMESH,
     +           IFAIL,EV,EF,PDEF,RHO,
     +           ELAPSE,FLAG)
C     ..
C     .. Scalar Arguments ..
      double precision A,B,A1,A2,B1,B2,TOLER
      double precision,intent(out)::ELAPSE
      integer ISOLVR,KHI,KLO,NLMESH,NXMESH,QCALC,FLAG
      character*4 ATYPE,BTYPE
C     ..
C     .. Array Arguments ..
      integer IFAIL(KLO:KHI)
      double precision LMESH(1:MAXMSH),XMESH(1:MAXMSH),
     +  EV(KLO:KHI,1:2),EF(1:MAXMSH*MAXEVS),PDEF(1:MAXMSH*MAXEVS),
     +  RHO(1:MAXMSH)
C SOLVIT oversees the passing of a given Problem to a given Solver with
C the extra defining info: type of calculation (shown by QCALC);e-v
C index range, tolerance etc. Inside it are Internal Procedures, one for
C each Solver, which set up the call to the solver in question.
C
C SOLVIT has a 'table of capabilities' by which it weeds out
C calculations that are deprecated (e.g. giving a non-AUTO x-mesh to
C SLEDGE) or impossible (e.g. asking any code except SLEDGE to do SDF
C calculation).
C - If an impossible calculation is requested it does not attempt it
C   and returns with FLAG=FORBID.
C - If a deprecated one is requested it asks if user wants to go ahead
C   anyway. If no it returns with FLAG=FORBID, if yes it attempts it and
C   returns with FLAG=DEPREC.
C - Otherwise it attempts it and returns with FLAG=OKCALC (=0).
C In this table blank entries mean impossible, (1,2) are deprecated:
C             +-evs only
C             |    +-evs, & efns on auto mesh
C             |    |    +-evs, & efns on uniform or user mesh
C             |    |    |    +-spectral density fn on unif/user mesh
C       QCALC=1    2    3    4
C     -------------------------
C     sledge  Y    Y   (1)   Y
C     sleign  Y         Y
C     sl02f   Y    Y    Y
C     sleig2  Y        (2)
C  (1) Pruess says *don't* allow QCALC=3)
C  (2) SLEIGN2 *can* do this but the code is in the SLEIGN2 driver & I
C      haven't worked out what to extract yet
C NOTE:The sort of thing that leads to FLAG>0 is asking SLEIGN to make
C      an AUTO-mesh. Maybe there is a way to make SLEIGN do this but I
C      haven't found it. On the other hand it was easy to form an
C      AUTO-mesh from SL02F's output, while SLEDGE can provide it
C      explicitly. That is, FLAG>0 may be a comment on the competence of
C      the person who wrote the interface code as much as on the solver.
C     ..
C     .. Array Parameters ..
      integer, parameter:: HOWOK(1:4,1:NSOLVR)=reshape(
     +  (/OKCALC, OKCALC, DEPREC, OKCALC,
     +    OKCALC, FORBID, OKCALC, FORBID,
     +    OKCALC, OKCALC, OKCALC, FORBID,
     +    OKCALC, FORBID, DEPREC, FORBID/)
     +  ,(/4,NSOLVR/))
C     ..
C     .. Local Scalars ..
      integer I,NXMSH0
      double precision UA,UB,PDUA,PDUB,X
C     ..
C     .. Local Arrays ..
      double precision XMESH0(1:MAXMSH)

C Input Arguments:
C ISOLVR     identifies which solver to use, in range 1..NSOLVR
C QCALC      =1 compute just evs, no efns
C            =2 compute evs, and efns on mesh automatically generated
C               by solver
C            =3 compute evs, and efns on user-supplied mesh
C            =4 compute spectral density fn (on user-supplied mesh)
C A,B,A1,A2,B1,B2,ATYPE,BTYPE
C            Endpoint & (regular) bdry cond info.
C            ATYPE,BTYPE are as returned by SLTSTPAK
C KLO,KHI    range of eigenvalues to find
C TOLER      tolerance (both abs & rel) passed to Solver
C NXMESH,XMESH (Only when QCALC.eq.3)
C            NXMESH=no. of points in user-supplied x-mesh.
C            NB! NXMESH may be as large as NXMESH+2 on output so must
C            have NXMESH<=MAXMSH-2 on input.
C            User-provided x-mesh is in XMESH(1:NXMESH)
C NLMESH,LMESH (Only when QCALC.eq.4)
C            NLMESH=no. of points in user-supplied lambda-mesh.
C            User-provided lambda-mesh is in LMESH(1:NLMESH)
C
C Output Arguments:
C   Note IFAIL is indexed KLO:KHI in this routine and in REPORT but
C   indexed 1:MAXEVS in main program where actually declared.
C   Similarly EV is indexed KLO:KHI,1:2 here and in REPORT
C   but indexed 1:2*MAXEVS in main program.
C   Similarly EF,PDEF are indexed 1:MAXMSH,KLO:KHI in REPORT
C   but indexed 1:(MAXMSH*MAXEVS) here and in main program.
C   !!NOTE!!
C   Since there is the possibility for NXMESH to be set by the Solver,
C   EF,PDEF can't have NXMESH in their dimension definition (within F77
C   rules) convenient tho' this would be.
C
C IFAIL      IFAIL(K) must hold Solver's error flag for calc of ev of
C            index K.
C          - There is one exception to this: if the input arguments make
C            an impossible request so that the Solver cannot even be
C            called (e.g. asking SLEIGN to produce an AUTO-mesh) then
C            set IFAIL(KLO:KHI) to 999. In this case please set all
C            EV(K,1) to 0.d0 and all EV(K,2) to 1.d20, and similar
C            'silly' values for other output arguments in case REPORT
C            has a bug & fails to report this fact.
C          - Routine SABORT can be used to do this.
C EV         EV(K,1) must hold computed ev of index K.
C            EV(K,2) must hold an absolute error estimate of the ev.
C            Depending on the solver, this may be an estimate of
C            (computed-true) or of a bound on abs(computed-true).
C NXMESH,XMESH
C            (only when QCALC.eq.2)
C            XMESH(1:NXMESH) is to hold solver-defined x-mesh.
C            (only when QCALC.eq.3)
C            XMESH(1:NXMESH) is to hold mesh after SOLVIT has
C            'truncated' it if necessary by removing points outside
C            current endpoints A,B and possibly re-inserting A,B.
C EF,PDEF    (only when QCALC.eq. 2 or 3)
C            Indexed 1:MAXMSH,KLO:KHI then EF(I,K),PDEF(I,K)
C            Indexed 1:...  then   EF((K-KLO)*MAXMSH+I),
C                                PDEF((K-KLO)*MAXMSH+I)
C            must hold value of
C            u(x),pu'(x) at I-th point of XMESH, for K-th efn
C ELAPSE     Total CPU time for Solver call
C FLAG       General error flag. Values on exit:
C =OKCALC(=0)  The solver was able to attempt the requested calculation
C              even if it failed to produce eigenvalues for some or
C              all K in KLO:KHI as shown by IFAIL(K) values.
C =DEPREC(=1)  The calculation was 'deprecated' but user chose to go
C              ahead anyway
C =FORBID(=2)  SOLVIT couldn't even attempt the requested calculation
C              with the current solver.

c      print*,'enter SOLVIT: NXMESH,XMESH=',NXMESH,XMESH


C Test the requested calculation against the 'capabilities table' and
C set FLAG accordingly:
      FLAG = HOWOK(QCALC,ISOLVR)
      if (FLAG.eq.FORBID) then
        write(*,*) '  ***Sorry: ',SLVNAM(ISOLVR),
     +   ' cannot provide the requested calculation'
      else if (FLAG.eq.DEPREC) then
        write(*,*)'  ***Note: this calculation is not recommended for ',
     +    SLVNAM(ISOLVR),' at present'
        if (YESNO(
     +     '     Results may be suspect, go ahead anyway?').eq.'n') then
          FLAG = FORBID
        end if
      end if
      if (FLAG.eq.FORBID) return
C ..else continue with normal processing ..

C If current endpoints are regular, either originally or by truncation
C of singular point
C   GETBCS will find coefficients (PDUA,UA), (PDUB,UB) which define
C   the regular BCs to be imposed at current endpoints, in case of
C   solvers that need this in their argument list.
C else
C   set (PDUA,UA), (PDUB,UB) to the values (A1,A2), (B1,B2) that
C   appear on the Main Menu
C end if
C The EIG argument is arbitrarily taken as 0.0. SLEDGE will extract
C what lambda-dependence it can support by making a further call to
C GETBCS with EIG=1.0.
      if (ATYPE.eq.'R') then
        call GETBCS(0,A,0d0,PDUA,UA)
      else
        PDUA = A1
        UA = A2
      end if
      if (BTYPE.eq.'R') then
        call GETBCS(1,B,0d0,PDUB,UB)
      else
        PDUB = B1
        UB = B2
      end if

C If QCALC=3, convert NXMESH,XMESH by removing points outside (A,B) and
C re-inserting A (sim. B) if it is regular or if solver is SLEDGE.
C NB! No. of points NXMESH may increase by up to 2 as a result!
C When checking whether a point is 'inside (A,B)', we require it to be
C actually in ((A+1e-6*|A|), (B-1e-6*|B|)) This avoids including end
C values twice if they have been read from the database and are a bit
C inaccurate. The factor 1e-6 seems sensible: the user must give >6 sig
C figs when he types the x values into the database.

      NXMSH0 = NXMESH
      do 90 I=1,NXMESH
        XMESH0(I) = XMESH(I)
   90 continue

      if (QCALC.eq.3) then
        NXMESH = 0
        if (ATYPE.eq.'R' .or. ISOLVR.eq.ISLED) then
          NXMESH = NXMESH+1
          XMESH(NXMESH) = A
        end if
        do 100 I=1,NXMSH0
          X = XMESH0(I)
          if (X.gt.(A+1d-6*abs(A)) .and. X.lt.(B-1d-6*abs(B))) then
            NXMESH = NXMESH+1
            XMESH(NXMESH) = X
          end if
  100   continue
        if (BTYPE.eq.'R' .or. ISOLVR.eq.ISLED) then
          NXMESH = NXMESH+1
          XMESH(NXMESH) = B
        end if
      end if

      write (*,*)
C See at bottom of TSTSET routine where dots are generated:
      write(*,*) 'Each dot represents 1000 function evaluations:'

      if (ISOLVR.eq.ISLED) then
        call SSLED(QCALC,A,B,ATYPE,BTYPE,
     +           KLO,KHI,TOLER,
     +           NXMESH,XMESH,NLMESH,LMESH,
     +           IFAIL,EV,EF,PDEF,RHO,
     +           ELAPSE,
     +           UA,PDUA,UB,PDUB)
      else if (ISOLVR.eq.ISLEIG) then
        call SSLEIG(QCALC,A,B,ATYPE,BTYPE,
     +           KLO,KHI,TOLER,
     +           NXMESH,XMESH,NLMESH,LMESH,
     +           IFAIL,EV,EF,PDEF,RHO,
     +           ELAPSE,
     +           UA,PDUA,UB,PDUB)
      else if (ISOLVR.eq.ISL02) then
        call SSL02(QCALC,A,B,ATYPE,BTYPE,
     +           KLO,KHI,TOLER,
     +           NXMESH,XMESH,NLMESH,LMESH,
     +           IFAIL,EV,EF,PDEF,RHO,
     +           ELAPSE,
     +           UA,PDUA,UB,PDUB)
      else if (ISOLVR.eq.ISLEI2) then
        call SSLEI2(QCALC,A,B,ATYPE,BTYPE,
     +           KLO,KHI,TOLER,
     +           NXMESH,XMESH,NLMESH,LMESH,
     +           IFAIL,EV,EF,PDEF,RHO,
     +           ELAPSE,
     +           UA,PDUA,UB,PDUB)
      else
        write (*,*) 'Invalid Solver-ID number requested: ',ISOLVR
        stop
      end if
C     To terminate the ...'s in TSTSET counting 1000 fn evals:
      write(*,*)
      end subroutine SOLVIT

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
      subroutine SSLED(QCALC,A,B,ATYPE,BTYPE,
     +           KLO,KHI,TOLER,
     +           NXMESH,XMESH,NLMESH,LMESH,
     +           IFAIL,EV,EF,PDEF,RHO,
     +           ELAPSE,
     +           UA,PDUA,UB,PDUB)
      use SLEDGEMD
      implicit none
C     ..
C     .. Scalar Arguments ..
      double precision A,B,TOLER,
     +           UA,PDUA,UB,PDUB
      double precision,intent(out)::ELAPSE
      integer KHI,KLO,NLMESH,NXMESH,QCALC
      character*4 ATYPE,BTYPE
C     ..
C     .. Array Arguments ..
      integer IFAIL(KLO:KHI)
      double precision LMESH(1:MAXMSH),XMESH(1:MAXMSH),
     +  EV(KLO:KHI,1:2),EF(1:MAXMSH*MAXEVS),PDEF(1:MAXMSH*MAXEVS),
     +  RHO(1:MAXMSH)

C     Interface to SLEDGE whose argument list is
C      subroutine SLEDGE(JOB,CONS,ENDFIN,INVEC,TOL,TYPE,EV,NUMX,XEF,EF,
C     +                  PDEF,T,RHO,IFLAG,STORE)
C     In the call, actual arguments have the same names as dummy ones
C     except as in the table:
C       Dummy    Actual
C       EV       EV(:,1)
C       XEF      XMESH
C       T        RHO
C       IFLAG    IFAIL

C     .. Parameters ..
C     ISTORE is workspace size, see description of STORE in SLEDGE code
      integer ISTORE
      parameter (ISTORE=26*MAXMSH+16)
C     ..
C     .. Local Scalars ..
      double precision PDUA1,UA1,CURREN
      integer I,IK,K,NUMX
C     ..
C     .. Local Arrays ..
      double precision CONS(1:8),STORE(1:ISTORE),TOL(1:6)
      integer INVEC(1:3+MAXEVS)
      logical ENDFIN(1:2),JOB(1:5),TYPE(1:4,1:2)
C     ..
C     .. External Subroutines ..
cc      external SLEDGE
C     ..

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**

C Set up endpoint and regular BC info for SLEDGE:
C     Note SLEDGE's B1,B2 is my PDUB,-UB
      CONS(1) = PDUA
      CONS(3) = UA
      CONS(5) = PDUB
      CONS(6) = -UB
C     At x=a, SLEDGE allows lambda-dependent BC of form
C       (A1-EIG*A1')u(a) = (A2-EIG*A2')pu'(a)
C     Extract A1',A2' by a further call to GETBCS with EIG=1.0, and
C     linearizing. But only if of type 'R', to avoid zerodivide etc in
C     evaluating BC fns in TSTSET
      if (ATYPE.eq.'R') then
        call GETBCS(0,A,1d0,PDUA1,UA1)
        CONS(2) = PDUA-PDUA1
        CONS(4) = UA-UA1
      else
        CONS(2) = 0D0
        CONS(4) = 0D0
      end if
      CONS(7) = A
      CONS(8) = B
      ENDFIN(1) = A .ne. -XINFTY
      ENDFIN(2) = B .ne. XINFTY
C Set correct pattern of JOB(1:5) settings for allowed QCALC values:
C (See SLEDGE documentation for details. Certain SLEDGE options can not
C be exercised, e.g. to compute evs & SDF at a single call.)
C QCALC                    \ JOB(I),I= 1  2  3  4  5
C----------------------------------------------------
C  1 evs, no efns                    | T  F  F  F  T
C  2 evs + efns AUTO x-mesh          | F  T  F  F  T
C  3 evs + efns USER x-mesh          | F  T  F  F  F
C  4 sp dens fn USER lambda-mesh     | F  F  T  F  T
C
      JOB(1) = QCALC.eq.1
      JOB(2) = QCALC.eq.2 .or. QCALC.eq.3
      JOB(3) = QCALC.eq.4
C JOB(4): automatic end-classif'n wanted
      JOB(4) = .FALSE.
C Choose SLEDGE-supplied (JOB(5)= .TRUE.) or user-supplied x-mesh:
      if (QCALC.ne.3) then
        JOB(5) = .TRUE.
        NUMX = 0
      else
C       Supplied mesh MUST contain A,B: I think I've achieved this!
        JOB(5) = .FALSE.
        NUMX = NXMESH
      end if

C Set printing level (0 unless debugging):
      INVEC(1) = 0
C Set no. of lambda-points for SDF calculation:
      INVEC(2) = NLMESH
C Number of eigenvalues to be computed:
      if (KHI-KLO+1.le.MAXEVS) then
         INVEC(3) = KHI - KLO + 1
      else
         write (*,FMT=*)
     +     '  Too many eigenvalues requested, reduced to ',MAXEVS
         INVEC(3) = MAXEVS
      end if
C Indices of eigenvalues sought:
C SLEDGE allows them to be non-contiguous but SOLVIT interface doesn't
      do 300 I = 1,INVEC(3)
         INVEC(3+I) = KLO + I - 1
  300 continue
C Set vector of tolerances:
      TOL(1) = TOLER
      TOL(2) = TOLER
      TOL(3) = TOLER
      TOL(4) = TOL(2)
      TOL(5) = TOLER
      TOL(6) = TOL(2)

c      print*,'before SLEDGE:NXMESH,XMESH(1:NXMESH)=',
c     +  NXMESH,XMESH(1:NXMESH)
C Reset the CPU clock
      call CPU(0,CURREN,ELAPSE)
C In this call, note the ev's are put in the 1st column of array EV:
      call SLEDGE(JOB,CONS,ENDFIN,INVEC,TOL,TYPE,EV,NUMX,XMESH,EF,PDEF,
     +            LMESH,RHO,IFAIL,STORE)

C Read the CPU clock
      call CPU(1,CURREN,ELAPSE)
C Update NXMESH since for QCALC.ne.3 NUMX is set by SLEDGE
      NXMESH = NUMX
c      print*,'after SLEDGE:NXMESH,XMESH(1:NXMESH)=',
c     +  NXMESH,XMESH(1:NXMESH)

C SLEDGE doesn't give error estimates as such so set the error-estimate
C entries of EV to indicate estimated mixed rel-abs bound of TOLER, or
C dummy value if failure exit.
      do 380 K=KLO,KHI
        if (IFAIL(K).eq.0) then
          EV(K,2) = TOLER*max(abs(EV(K,1)),1D0)
        else
          EV(K,2) = 1D20
        end if
  380 continue
C Output SLEDGE's endpoint-classification data (as REPORT doesn't know
C about this). It's supposed to be independent of K so I assume checking
C IFAIL(KLO) tells me whether SLEDGE has failed to produce it (IFAIL<0)
C or believes it is not to be trusted (IFAIL contains the digit 2)
        write(*,FMT=9999) 'A',(TYPE(I,1),I=1,4), 'B',(TYPE(I,2),I=1,4)
 9999   format(/,'Classification of endpoint ',a,': Regular=',L1,
     +  ', LC =',L1,', Nonosc all EV=',L1,', Osc all EV=',L1)
        write(20,FMT=9998) 'A',(TYPE(I,1),I=1,4), 'B',(TYPE(I,2),I=1,4)
 9998   format(' %Classification of endpoints: ',2(1x,a,' =',4L2))
        if (HASDIG(IFAIL(KLO),2)) then
          write(*,*) 'There is doubt about this (IFAIL=2)'
          write(20,*) ' %There is doubt about this (IFAIL=2)'
        end if
        call SPAUSE

C If efns were computed (JOB(2) is true) relocate their values from
C SLEDGE's storage (consecutive blocks of length NUMX) to REPORT's
C storage (EF,PDEF treated as dimensioned (1:MAXMSH,KLO:KHI)).
C Work backwards as we are rearranging storage on top of itself.
C However if IFAIL(K)<0 (SLEDGE fatal error) set all the values
C to 0.
      if (JOB(2)) then
        do 410 K=KHI,KLO,-1
          IK = K-KLO
          if (IFAIL(K).ge.0) then
            do 400 I=NUMX,1,-1
              EF(IK*MAXMSH+I) = EF(IK*NUMX+I)
              PDEF(IK*MAXMSH+I) = PDEF(IK*NUMX+I)
  400       continue
          else
            do 405 I=NUMX,1,-1
              EF(IK*MAXMSH+I) = 0D0
              PDEF(IK*MAXMSH+I) = 0D0
  405       continue
          end if
  410   continue
      end if
      end subroutine SSLED

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function HASDIG(N,D)
      implicit none
      integer N,D
C To help interpret SLEDGE's positive IFAIL values.
C Gives TRUE iff D is one of the decimal digits of N
C (1041 has digits 0,1,4; so does -1041; 0 has digit 0)
      integer NN
      NN = abs(N)
   10 if (mod(NN,10).eq.D) goto 20
      NN = NN/10
      if (NN.ne.0) go to 10
      HASDIG = .FALSE.
      return
   20 HASDIG = .TRUE.
      end function HASDIG
      
C---+----|----+----|----+----|----+----|----+----|----+----|----+----|--
      subroutine SSLEIG(QCALC,A,B,ATYPE,BTYPE,
     +           KLO,KHI,TOLER,
     +           NXMESH,XMESH,NLMESH,LMESH,
     +           IFAIL,EV,EF,PDEF,RHO,
     +           ELAPSE,
     +           UA,PDUA,UB,PDUB)
      use SLEIGNMD
      use PQWCOM
      use TESTMOD,only:ABINFO
      implicit none
C     ..
C     .. Scalar Arguments ..
      double precision A,B,TOLER,
     +           UA,PDUA,UB,PDUB
      double precision,intent(out)::ELAPSE
      integer KHI,KLO,NLMESH,NXMESH,QCALC
      character*4 ATYPE,BTYPE
C     ..
C     .. Array Arguments ..
      integer IFAIL(KLO:KHI)
      double precision LMESH(1:MAXMSH),XMESH(1:MAXMSH),
     +  EV(KLO:KHI,1:2),EF(1:MAXMSH*MAXEVS),PDEF(1:MAXMSH*MAXEVS),
     +  RHO(1:MAXMSH)

C     Interface to SLEIGN whose argument list is
C      SUBROUTINE SLEIGN(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,
C     1                  NUMEIG,EIG,TOL,IFLAG,ISLFUN,SLFUN)
C     In the call, actual arguments have the same names as dummy ones
C     except as in the table:
C       Dummy    Actual
C       A1       PDUA
C       A2       -UA
C       B1       PDUB
C       B2       -UB

C     ..
C     .. Local Scalars ..
      integer I,IFLAG,IK,INTAB,ISLFUN,J,K,NUMEIG,AOFSET,BOFSET
      double precision CURREN
      double precision DUM,EIG,TOL,P0ATA,QFATA,P0ATB,QFATB
C     ..
C     .. Local Arrays ..
      double precision SLFUN(1:9+MAXMSH)
C     ..
C     .. External Functions ..
      double precision P
      external P
C     ..
C     .. External Subroutines ..
cc      external SLEIGN

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C NOTE: SLEIGN's notation differs from my standard in these ways:
C     Its NUMEIG is (my K)+1
C     Its A1,A2 are my PDUA,-UA & same for B1,B2
C     Its q(x) is -(my q(x))
C     His coeff functions (my p,-q,w) MUST be called p,q,r
C
C SLEIGN has no way (that I know) to generate AUTO mesh for efns:
      if (QCALC.eq.2) then
        write (*,9999)
 9999   format('SLEIGN interface can''t provide AUTO-mesh, please '/
     +    'select a different option in EFMENU or another solver')
C!        call SABORT
        return
      end if
C Make Q(X) flip the sign of q:
      FLIPQ = .TRUE.
C Set up endpoint info for SLEIGN
C     Say if finite or infinite endpoints by INTAB=1,2,3 or 4:
      INTAB = 1
      if (B.eq.XINFTY) INTAB = INTAB + 1
      if (A.eq.-XINFTY) INTAB = INTAB + 2

C Say if p(x)=0 and q(x)=finite at endpoints.
C     The ABINFO routine gives the data for *untruncated* endpoints:
      call ABINFO(P0ATA,QFATA,P0ATB,QFATB)
C     .. but if end is regular, originally or by *truncation* then the
C     following apply:
      if (ATYPE.eq.'R') then
        P0ATA = -1D0
        QFATA = 1D0
      end if
      if (BTYPE.eq.'R') then
        P0ATB = -1D0
        QFATB = 1D0
      end if

C Dummy evaluation of coefficient function, see code of P,Q,R.
C     Use the x/(1+|x|) map & its inverse t/(1-|t|)
C     to safely find an interior point of (A,B)
      DUM = 0.5D0*(A/(1+abs(A)) + B/(1+abs(B)))
      DUM = P(DUM/(1-abs(DUM)))

C     Reset the CPU clock:
      call CPU(0,CURREN,ELAPSE)

C Loop over eigenvalues:
      do 300 K = KLO,KHI
C       Make SLEIGN choose own guess of EIG:
        EIG = 0.0D0
C       Make local copy of TOLER (since changed by SLEIGN):
        TOL = TOLER
C       Make local copy of K (since may be changed by SLEIGN)
C       Also SLEIGN indexes evs from 1, not 0
        NUMEIG = K+1

C       Control eigenfunction calculation.
C       If ISLFUN>0 then:
C       SLFUN(9+1:9+ISLFUN) holds x(i) on entry, u(x(i)) on exit
        if (QCALC.eq.1) then
C         Evs only
          ISLFUN = 0
        else if (QCALC.eq.2) then
C         We should have excluded this so:
          write(*,*) 'ERROR in SSLEIG: QCALC=2 shouldn''t get here!'
          stop
        else
C         Efns computed on USER-mesh
C         Careful! All my meshes contain A,B. But SLEIGN bombs out if
C         given a singular endpoint as part of its mesh. So discard
C         XMESH(1) and/or XMESH(NXMESH) if singular, and insert dummy
C         e-fn values (zero) in the output array, after the solver call.
C         NB: moving the end meshpoints inwards a bit is not an option,
C         as the results would not be comparable with those produced by
C         other solvers
          AOFSET = 0
          BOFSET = 0
          if (ATYPE.ne.'R') AOFSET = 1
          if (BTYPE.ne.'R') BOFSET = 1
          ISLFUN = NXMESH-AOFSET-BOFSET
          do 110 I=1,ISLFUN
            SLFUN(9+I) = XMESH(I+AOFSET)
  110     continue
        end if

        call SLEIGN(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,PDUA,-UA,
     +               PDUB,-UB,NUMEIG,EIG,TOL,IFLAG,ISLFUN,SLFUN)

C Extract SLEIGN's 'out' or 'inout' arguments:
        IFAIL(K) = IFLAG
C       Get ev estimate, & estimate of absolute error bound:
        EV(K,1) = EIG
        EV(K,2) = max(1D0,abs(EIG))*TOL
C       Get efn values u(x(i)); note SLEIGN doesn't provide pu' values
C       Insert dummy u=0 at singular ends, see "Careful!" above
        if (QCALC.ge.2) then
c          print*,'SLFUN(1:9+NXMESH)=',SLFUN(1:9+NXMESH)
          IK = (K-KLO)*MAXMSH
          do 120 J=1,NXMESH
            I = J-AOFSET
            if (I.ge.1 .and. I.le.ISLFUN) then
              EF(IK+J) = SLFUN(9+I)
            else
              EF(IK+J) = 0D0
            end if
            PDEF(IK+J) = 0D0
  120     continue
        end if
C       Report if SLEIGN altered NUMEIG:
        if (NUMEIG-1.ne.K) then
          write (*,*) 'It appears eigenvalue',K,' does not exist'
          write (*,*) 'Estimated index of last eigenvalue is',NUMEIG-1
        end if

  300 continue

C     Read the CPU clock:
      call CPU(1,CURREN,ELAPSE)

      end subroutine SSLEIG

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine SSLEI2(QCALC,A,B,ATYPE,BTYPE,
     +           KLO,KHI,TOLER,
     +           NXMESH,XMESH,NLMESH,LMESH,
     +           IFAIL,EV,EF,PDEF,RHO,
     +           ELAPSE,
     +           UA,PDUA,UB,PDUB)
      use SLEIG2MD
      use PQWCOM
      use TESTMOD,only:ABINFO
      implicit none
C     ..
C     .. Scalar Arguments ..
      double precision A,B,TOLER,
     +           UA,PDUA,UB,PDUB
      double precision,intent(out)::ELAPSE
      integer KHI,KLO,NLMESH,NXMESH,QCALC
      character*4 ATYPE,BTYPE
C     ..
C     .. Array Arguments ..
      integer IFAIL(KLO:KHI)
      double precision LMESH(1:MAXMSH),XMESH(1:MAXMSH),
     +  EV(KLO:KHI,1:2),EF(1:MAXMSH*MAXEVS),PDEF(1:MAXMSH*MAXEVS),
     +  RHO(1:MAXMSH)

C     Interface to SLEIGN2 whose argument list is
C      subroutine SLEIGN2(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,A1,A2,B1,B2,
C     +                   NUMEIG,EIG,TOL,IFLAG,ISLFUN,SLFUN,SINGATA,
C     +                   SINGATB,CIRCLA,CIRCLB,OSCILA,OSCILB)
C     In the call, actual arguments have the same names as dummy ones
C     except as in the table:
C       Dummy    Actual
C       A1       PDUA
C       A2       -UA
C       B1       PDUB
C       B2       -UB

C     ..
C     .. Local Scalars ..
      integer I,IFLAG,IK,INTAB,ISLFUN,J,K,NUMEIG,AOFSET,BOFSET
      double precision CURREN
      double precision DUM,EIG,TOL,P0ATA,QFATA,P0ATB,
     +                 QFATB,SINGATA,SINGATB,CIRCLA,CIRCLB,OSCILA,OSCILB
C     ..
C     .. Local Arrays ..
      double precision SLFUN(1:9+MAXMSH)
C     ..
C     .. External Functions ..
      double precision P
      external P
C     ..
C     .. External Subroutines ..
cc      external SLEIGN2

C Special COMMON for counting no. of calls to SLEIGN2's UV routine
      common/SLEI2C/ NUVEVL
      integer NUVEVL

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C NOTE: SLEIGN2's notation differs from my standard in these ways:
C     Its A1,A2 are my PDUA,-UA & same for B1,B2
C
C Make Q(X) NOT flip the sign of q:
      FLIPQ = .FALSE.
C SLEIGN2 has no way (that I know) to generate AUTO mesh for efns:
      if (QCALC.eq.2) then
        write (*,9999)
 9999   format('SLEIGN2 interface can''t provide AUTO-mesh, please '/
     +    'select a different option in EFMENU or another solver')
C!        call SABORT
        return
      end if
C Set up endpoint info for SLEIGN2
C     Say if finite or infinite endpoints by INTAB=1,2,3 or 4:
      INTAB = 1
      if (B.eq.XINFTY) INTAB = INTAB + 1
      if (A.eq.-XINFTY) INTAB = INTAB + 2

C Say if p(x)=0 and q(x)=finite at endpoints.
C     The ABINFO routine gives the data for *untruncated* endpoints:
      call ABINFO(P0ATA,QFATA,P0ATB,QFATB)
C     .. but if end is regular, originally or by *truncation* then the
C     following apply:
      if (ATYPE.eq.'R') then
        P0ATA = -1D0
        QFATA = 1D0
      end if
      if (BTYPE.eq.'R') then
        P0ATB = -1D0
        QFATB = 1D0
      end if

C     Classify type of singularity if any:
      if (ATYPE.eq.'R' .or. ATYPE.eq.'WR') then
         SINGATA = -1.0D0
         CIRCLA  = -1.0D0
         OSCILA  = -1.0D0
      else if (ATYPE.eq.'LCN') then
         SINGATA =  1.0D0
         CIRCLA  =  1.0D0
         OSCILA  = -1.0D0
      else if (ATYPE.eq.'LCO') then
         SINGATA =  1.0D0
         CIRCLA  =  1.0D0
         OSCILA  =  1.0D0
      else
C     for LP, LPNO it doesn't matter but:
         SINGATA =  1.0D0
         CIRCLA  = -1.0D0
         OSCILA  = -1.0D0
      end if

      if (BTYPE.eq.'R' .or. BTYPE.eq.'WR') then
         SINGATB = -1.0D0
         CIRCLB  = -1.0D0
         OSCILB  = -1.0D0
      else if (BTYPE.eq.'LCN') then
         SINGATB =  1.0D0
         CIRCLB  =  1.0D0
         OSCILB  = -1.0D0
      else if (BTYPE.eq.'LCO') then
         SINGATB =  1.0D0
         CIRCLB  =  1.0D0
         OSCILB  =  1.0D0
      else
C     for LP, LPNO it doesn't matter but:
         SINGATB =  1.0D0
         CIRCLB  = -1.0D0
         OSCILB  = -1.0D0
      end if

C Dummy evaluation of coefficient function, see code of P,Q,R.
C     Use the x/(1+|x|) map & its inverse t/(1-|t|)
C     to safely find an interior point of (A,B)
      DUM = 0.5D0*(A/(1+abs(A)) + B/(1+abs(B)))
      DUM = P(DUM/(1-abs(DUM)))

C     Reset the CPU clock:
      call CPU(0,CURREN,ELAPSE)
C     Zero the counter of calls to SLEIGN2's UV routine:
      NUVEVL = 0
C See at bottom of UV routine in problem set, where dots are generated:
      write(*,*)
     + 'Each + represents 1000 calls to SLEIGN2''s UV routine:'

C Loop over eigenvalues:
      do 300 K = KLO,KHI
C       Make SLEIGN2 choose own guess of EIG:
        EIG = 0.0D0
C       Make local copy of TOLER (since changed by SLEIGN2):
C       Minus sign requests trace-output
        TOL = -TOLER
C       Make local copy of K (since may be changed by SLEIGN2)
        NUMEIG = K

C       Control eigenfunction calculation.
C       If ISLFUN>0 then:
C       SLFUN(9+1:9+ISLFUN) holds x(i) on entry, u(x(i)) on exit
        if (QCALC.eq.1) then
C         Evs only
          ISLFUN = 0
        else if (QCALC.eq.2) then
C         We should have excluded this so:
          write(*,*) 'ERROR in SSLEI2: QCALC=2 shouldn''t get here!'
          stop
        else
C         Efns computed on USER-mesh
C         See "Careful!" comments in SLEIGN interface.
C         Same applies here except we shouldn't be using this way
C         to compute e-fns, which seems faulty for SLEIGN2
          AOFSET = 0
          BOFSET = 0
          if (ATYPE.ne.'R') AOFSET = 1
          if (BTYPE.ne.'R') BOFSET = 1
          ISLFUN = NXMESH-AOFSET-BOFSET
          do 110 I=1,ISLFUN
            SLFUN(9+I) = XMESH(I+AOFSET)
  110     continue
        end if

c Testprints:
c      print*,'SSLEI2: A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,PDUA,-UA=',
c     + A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,PDUA,-UA
c      print*,'SSLEI2: PDUB,-UB,NUMEIG,EIG,TOL,IFLAG,ISLFUN=',
c     + PDUB,-UB,NUMEIG,EIG,TOL,IFLAG,ISLFUN
c      print*,'SSLEI2: SINGATA,CIRCLA,OSCILA,SINGATB,CIRCLB,OSCILB=',
c     + SINGATA,CIRCLA,OSCILA,SINGATB,CIRCLB,OSCILB
c      call SPAUSE
c
        call SLEIGN2(A,B,INTAB,P0ATA,QFATA,P0ATB,QFATB,PDUA,-UA,
     +               PDUB,-UB,NUMEIG,EIG,TOL,IFLAG,ISLFUN,SLFUN,
     +               SINGATA,SINGATB,CIRCLA,CIRCLB,OSCILA,OSCILB)

C Extract SLEIGN2's 'out' or 'inout' arguments:
C       It seems a multiple of 100 is added to the 'normal' flag value
C       to give extra info to SLEIGN2's driver. Remove this:
        IFAIL(K) = modulo(IFLAG,100)
C       Get ev estimate, & estimate of absolute error bound:
        EV(K,1) = EIG
        EV(K,2) = max(1D0,abs(EIG))*TOL

C       Get efn values u(x(i));
C       Insert dummy u=0 at singular ends
C       .. see "Careful!" above in this routine & in SLEIGN interface
        if (QCALC.ge.2) then
c          print*,'SLFUN(1:9+NXMESH)=',SLFUN(1:9+NXMESH)
          IK = (K-KLO)*MAXMSH
          do 120 J=1,NXMESH
            I = J-AOFSET
            if (I.ge.1 .and. I.le.ISLFUN) then
              EF(IK+J) = SLFUN(9+I)
            else
              EF(IK+J) = 0D0
            end if
            PDEF(IK+J) = 0D0
  120     continue
        end if
C       Report if SLEIGN2 altered NUMEIG:
        if (NUMEIG.ne.K) then
          write(*,*)
          write (*,*) 'It appears eigenvalue',K,' does not exist'
          write (*,*) 'Estimated index of last eigenvalue is',NUMEIG
        end if
        write (*,*) 'K = ',K,' done'

  300 continue

C     Read the CPU clock:
      call CPU(1,CURREN,ELAPSE)
      write (*,fmt=
     + '(1X,''No. of calls to UV routine: '',i8)') NUVEVL

      end subroutine SSLEI2

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine SSL02(QCALC,A,B,ATYPE,BTYPE,
     +           KLO,KHI,TOLER,
     +           NXMESH,XMESH,NLMESH,LMESH,
     +           IFAIL,EV,EF,PDEF,RHO,
     +           ELAPSE,
     +           UA,PDUA,UB,PDUB)
      use MARCOMOD
C     Needed for SELMSH routine:
      use SLUTIL
      use SL02CM
      implicit none
C     ..
C     .. Scalar Arguments ..
      double precision A,B,TOLER,
     +           UA,PDUA,UB,PDUB
      double precision,intent(out)::ELAPSE
      integer KHI,KLO,NLMESH,NXMESH,QCALC
      character*4 ATYPE,BTYPE
C     ..
C     .. Array Arguments ..
      integer IFAIL(KLO:KHI)
      double precision LMESH(1:MAXMSH),XMESH(1:MAXMSH),
     +  EV(KLO:KHI,1:2),EF(1:MAXMSH*MAXEVS),PDEF(1:MAXMSH*MAXEVS),
     +  RHO(1:MAXMSH)

C     Interface to SL02F whose argument list is
C      subroutine SL02F(ELAM,A,B,K,AINFO,BINFO,SYM,TOL,COEFFN,SETUP,N,WK,
C     &                 IWK,WKSMAL,ISMAL,KNOBS,IFAIL)
C     In the call, actual arguments have the same names as dummy ones
C     except as in the table:
C       Dummy    Actual
C       A        ALOC
C       B        BLOC
C       TOL      TOLER
C       IFAIL    IFAILA

C     .. Parameters ..
      integer IWK,ISMAL
      parameter (IWK=8000,ISMAL=130)
C     ..
C     .. Local Scalars ..
      double precision ALOC,BLOC,CURREN
      integer I,IFAILA,IK,K,N
      logical SYM,HAVEIG
C     ..
C     .. Local Arrays ..
      double precision EIGFNS(0:1,1:3),ELAM(1:2),
     +        WK(0:IWK,1:4),WKSMAL(0:ISMAL,1:7)
      integer KNOBS(1:2)
C     ..
C     .. External Subroutines ..
cc      external SETUP,SL02F
      external SETUP
C     ..
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Set up endpoint info for SL02F
      if (ATYPE.eq.'R') then
         AINFO = 'R'
      else if (ATYPE.eq.'LCN' .and. A.gt.-XINFTY) then
         AINFO = 'S'
      else if (ATYPE.eq.'LCO') then
         AINFO = 'S'
         write (6,FMT=*)
     +     ' ** Endpoint b is of type LCO which SL02F cannot handle',
     +     ' but will have a go'
      else if (ATYPE(1:2).eq.'LP' .and. A.gt.-XINFTY) then
         AINFO = 'S'
      else if (ATYPE.eq.'WR' .and. A.gt.-XINFTY) then
         AINFO = 'S'
      else
         AINFO = 'I'
      end if

      if (BTYPE.eq.'R') then
         BINFO = 'R'
      else if (BTYPE.eq.'LCN' .and. B.lt.XINFTY) then
         BINFO = 'S'
      else if (BTYPE.eq.'LCO') then
         BINFO = 'S'
         write (6,FMT=*)
     +     ' ** Endpoint b is of type LCO which SL02F cannot handle',
     +     ' but will have a go'
      else if (BTYPE(1:2).eq.'LP' .and. B.lt.XINFTY) then
         BINFO = 'S'
      else if (BTYPE.eq.'WR' .and. B.lt.XINFTY) then
         BINFO = 'S'
      else
         BINFO = 'I'
      end if
C     At present SOLVIT isn't passing symmetry info across so:
      SYM = .false.

C     Reset the CPU clock and the function evaluation counter:
      call CPU(0,CURREN,ELAPSE)

C Loop making separate SL02F call for each eigenvalue index:
      do 200 K = KLO,KHI
         IK = K-KLO

C        Set size of storage given to SL02F in WK array
C        (normally reset to what was used, on exit)
         N = IWK

         KNOBS(1) = 0
         KNOBS(2) = 0
C Always use 0 as initial guess of ev, and 1 as initial search step:
         ELAM(1) = 0.D0
         ELAM(2) = 1.D0
C Soft fail
         IFAILA = -1
C Use copies of endpoints as SL02F alters them at singular points:
         ALOC = A
         BLOC = B

        call SL02F(ELAM,ALOC,BLOC,K,AINFO,BINFO,SYM,TOLER,COEFFN,
     +             SETUP,N,WK,IWK,WKSMAL,ISMAL,KNOBS,IFAILA)

C Copy output into arrays for output:
        EV(K,1) = ELAM(1)
        EV(K,2) = ELAM(2)
        IFAIL(K) = IFAILA
C Report this for general interest:
        if (IFAILA.eq.0 .or. N.ne.IWK) then
          write(*,*)'K=',K,', error flag =',IFAILA,
     +      ', no. of mesh intervals used=',N
        else
          write(*,*)'K=',K,', error flag =',IFAILA,
     +      ', SL02F was unable to form a satisfactory mesh'
        end if

C If SL02F returned anything other than IFAIL=0 or 12 then no sensible
C ev estimate was found so no point attempting eigenfunctions.
        HAVEIG = (IFAILA.eq.0 .or. IFAILA.eq.12)

        if (QCALC.eq.2 .and. K.eq.KLO) then
          if (HAVEIG) then
C***Do this once, for K=KLO only*** (irrelevant when SL02FM used)
C In case of AUTO-mesh, return, in NXMESH & XMESH, a mesh derived from
C that used by SL02F. If the latter has NN *interior* meshpoints and
C NN<=MAXMSH, use all of them, else use a selection got by linear
C map from indices 1:NXMESH to 1:NN, with A,B added.
C   Note that the *interior* meshpoints are WK(ia:ib,1) where
C   ia is 0 if A is singular, 1 if regular
C   ib is N if B is singular, N-1 if regular
            call SELMSH(N,WK(0,1),MAXMSH,NXMESH,XMESH)
            if (AINFO.eq.'R') XMESH(1) = A
            if (BINFO.eq.'R') XMESH(NXMESH) = B
          else
C If no sensible mesh around, return a trivial mesh:
C ???faulty!!!
            NXMESH = 2
            XMESH(1) = A
            XMESH(2) = B
          end if
        end if

        if (QCALC.ge.2) then
C If efns are wanted, then either XMESH held USER-mesh on entry, or
C AUTO-mesh has just been created.
C Call MARCOPAK routines to evaluate normalized eigenfunction on XMESH.

          if (HAVEIG) then
C           ***Informative fail while debugging***
            IFAILA = -1
            call SL03F(ELAM,K,EIGFNS,WK,IWK,N,SETUP,IFAILA)
            IFAILA = -1
            call SL04F(XMESH,EF(IK*MAXMSH+1),PDEF(IK*MAXMSH+1),
     +        MAXMSH-1,ELAM,EIGFNS,WK,IWK,NXMESH-1,N,IFAILA)
          else
C           If no sensible mesh, set all u & pu' values to 0
            write (*,*)
     +      'No eigenfunctions available, so zeros will be displayed'
            do 180 I=1,NXMESH
              EF(IK*MAXMSH+I) = 0D0
              PDEF(IK*MAXMSH+I) = 0D0
  180       continue
          end if
        end if

  200 continue

C Read the CPU clock:
      call CPU(1,CURREN,ELAPSE)

      end subroutine SSL02

      end module SOLVRS

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C             Auxiliary routines for some of the solvers
C These can't be in a module until the solvers themselves are made into
C modules since they must be F77 external names to the linker.
C for SL02F **+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine SETUP(Y,PDY,EIG,X,IEND,ISING)
C SETUP interfaces between SL02F and GETBCS, putting the arguments in a
C different order. It also extracts the ISING data from SSL02 via
C SL02CM.
      use SL02CM
      use SLTSTPAK,only:GETBCS
      implicit none
C     .. Scalar Arguments ..                                            
      double precision EIG,PDY,X,Y
      integer IEND
      logical ISING
C     ..
      call GETBCS(IEND,X,EIG,PDY,Y)
      if (IEND.eq.0) then
        ISING = AINFO.ne.'R'
      else
        ISING = BINFO.ne.'R'
      end if
c      write(*,FMT=9999)IEND,X,EIG,PDY,Y,ISING
c 9999 format('SETUP iend=',i1,' x=',1pg12.6,' eig=',g12.6,' pdy,y=',
c     + 2g12.6,' ising=',l1)

      end subroutine SETUP

C for SLEDGE: +****|****+****|****+****|****+****|****+****|****+****|**
      subroutine COEFF(X,PX,QX,RX)
      use SLTSTPAK,only:COEFFN
      implicit none
C This just interfaces to COEFFN to give the name required by SLEDGE
C     .. Scalar Arguments ..
      double precision PX,QX,RX,X
C     ..
      call COEFFN(X,PX,QX,RX)
      end subroutine COEFF

C for SLEIGN, SLEIGN2: *+****|****+****|****+****|****+****|****+****|**
C Functions for p(x), q(x), w(x) saving redundant calls to COEFFN:
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      double precision function P(X)
      use PQWCOM
      use SLTSTPAK,only:COEFFN
      implicit none
      double precision X

      if (X.ne.XX) then
         XX = X
         call COEFFN(XX,PP,QQ,WW)
      end if

      P = PP
      end
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C This is my q(x) in SLEIGN2, but minus (my q(x)) in SLEIGN:
      double precision function Q(X)
      use PQWCOM
      use SLTSTPAK,only:COEFFN
      implicit none
      double precision X

      if (X.ne.XX) then
         XX = X
         call COEFFN(XX,PP,QQ,WW)
      end if

      Q = QQ
      if (FLIPQ) Q=-Q
      end
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      double precision function W(X)
      use PQWCOM
      use SLTSTPAK,only:COEFFN
      implicit none
      double precision X

      if (X.ne.XX) then
         XX = X
         call COEFFN(XX,PP,QQ,WW)
      end if

      W = WW
      end

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C SLEIGN's name for w(x) is r(x):
      double precision function R(X)
      use PQWCOM
      use SLTSTPAK,only:COEFFN
      implicit none
      double precision X

      if (X.ne.XX) then
         XX = X
         call COEFFN(XX,PP,QQ,WW)
      end if

      R = WW
      end

C for SLEIGN2: ****|****+****|****+****|****+****|****+****|****+****|**
      subroutine UV(XA,UA,PUP,VA,PVP,HU,HV)
      use TESTMOD,only: A1,A2,A,B1,B2,B,EIG,GTHUHV
      use SLTSTPAK,only:GETBCS
      implicit none
      double precision XA,UA,PUP,VA,PVP,HU,HV
C UV computes the boundary condition information for SLEIGN2 to use at
C a LCN or LCO endpoint.
C Input argument
C   XA     Value of independent variable x
C Output arguments
C   UA,PUP Values of u(x), pu'(x) where u is first maximal domain function
C         defining BCs
C   VA,PVP Values of v(x), pv'(x) where v is 2nd maximal domain function
C         defining BCs
C   HU,HV Values of -(pu')'+qu, -(pv')'+qv at x
C If functions used at the left end are different from those at the
C right end, this must be explicitly coded by choosing a breakpoint c
C and writing code of the form
C   if (XA.lt.C) then
C     functions for left end
C   else
C     functions for right end
C   end if

C     .. Local Variables ..
      integer IEND
      double precision A1SAV,A2SAV,B1SAV,B2SAV,C,TA,TB,TC

C Special COMMON for counting no. of calls to this routine
      common/SLEI2C/ NUVEVL
      integer NUVEVL

C     .. Executable Statements ..
C Form arbitrary 'mid-point' of interval
      TA = A/(1D0+abs(A))
      TB = B/(1D0+abs(B))
      TC = 0.5D0*(TA+TB)
      C = TC/(1D0-abs(TC))

C Save SLTSTPAK state info, mess around with it, then restore it:
      if (XA.le.C) then
        IEND = 0
        A1SAV = A1
        A2SAV = A2
C      Set coeffs for 1st BC function:
        A1 = 1D0
        A2 = 0D0
        call GETBCS(IEND,XA,EIG,PUP,UA)
C      Set coeffs for 2nd BC function:
        A1 = 0D0
        A2 = 1D0
        call GETBCS(IEND,XA,EIG,PVP,VA)
        A1 = A1SAV
        A2 = A2SAV
      else
        IEND = 1
        B1SAV = B1
        B2SAV = B2
C      Set coeffs for 1st BC function:
        B1 = 1D0
        B2 = 0D0
        call GETBCS(IEND,XA,EIG,PUP,UA)
C      Set coeffs for 2nd BC function:
        B1 = 0D0
        B2 = 1D0
        call GETBCS(IEND,XA,EIG,PVP,VA)
        B1 = B1SAV
        B2 = B2SAV
      end if
c      print*,'UV: C,IEND,XA,EIG,PUP,UA,PVP,VA=',C,IEND,XA,EIG,PUP,UA,PVP,VA

C Now the HU,HV info from Problem Set module.
      call GTHUHV(IEND,XA,UA,VA,HU,HV)

      NUVEVL = NUVEVL + 1
      if (mod(NUVEVL,1000).eq.0) write(*,fmt='("+")',advance='NO')
c      print*,'XA,UA,PUP,VA,PVP,HU,HV=', XA,UA,PUP,VA,PVP,HU,HV
      return
      end subroutine UV

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C!      subroutine SABORT
C!      implicit none
C!C To be called by one of the Solver interface routines in case of an
C!C impossible request, see description of IFAIL in SOLVIT
C!C Probably superseded at version 3.0 by SOLVIT's 'capabilities table'
C!C SABORT calls left in, commented out., as defensive measure.
C!      integer K
C!      do 10 K=KLO,KHI
C!        IFAIL(K) = 999
C!        EV(K,1) = 0d0
C!        EV(K,2) = 1d20
C!   10 continue
C!      FLAG = 1
C!      end subroutine SABORT

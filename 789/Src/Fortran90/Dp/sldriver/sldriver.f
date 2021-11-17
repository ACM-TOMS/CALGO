C***+****|****+****|* COPYRIGHT J D PRYCE 1998 **|****+****|****+****|**
C          SLDRIVER Version 4.1 by John D Pryce, June 1998
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C
C   To go with SLTSTPAK Version 3.2 Jun 98.
C   Driver Program for SL-Solvers and my design of Test Set implemented
C   in SLTSTPAK package. Includes facilities for computing both eigen-
C   values (e-vs) and eigenfunctions (e-fns) and for comparing computed
C   e-vs and e-fns with "true" values held in a database file.
C
C Revision history:
C   Apr 1994: Version 3
C   Aug 94: Revised by Steve Pruess
C   Sep-Dec 94: Improved screen appearance and summary output
C   Nov-Dec 95: Make compatible with Salford FTN90 compiler
C   Jan-Apr 96: Major revision involving:
C   1.Add the database file facilities.
C     This involved redefining the state-variables, which describe what
C     data is/is not present and cause appropriate action to be taken
C     when various options are invoked from the menus.
C   2.Divide the code into
C     - Driver (this module).
C     - DBMOD module to read "true" e-v and e-fn data
C     - The SAFEIO module
C     - SLTSTPAK
C     - The current Test Set of problems
C       (JDPTSET is the standard set of 60 problems)
C     - SOLVRS, a general interface module for several Solvers. This
C       contains the names of the solvers and other information as
C       PARAMETER data, and a general interface routine SOLVIT which
C       CONTAINS specific interface routines:
C       - SSLEIG which calls SLEIGN
C       - SSLED which calls SLEDGE
C       - etc
C     Thus adding a new Solver involves adding the data items to SOLVRS
C     and adding a new specific interface routine inside SOLVIT. Of
C     course, this interface routine is the key!! It must connect
C     whatever special facilities the Solver offers, to the information
C     available from SLTSTPAK routines.
C
C   Mar 1997. In response to comments from TOMS referee.
C     1.Add Sp Dens Fn facility offered by SLEDGE. This involved a
C       revision of the main menu & converting 'e-fn choices' submenu to
C       'calculation choices' submenu.
C     2.Make database interface simpler & more maintainable (more object
C       oriented).
C     3.Tidy up the routine documentation.
C     4.Add more info about the problem being solved, to the REPORT
C       output
C
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C                        GENERAL DOCUMENTATION
C!..yet to be written


C***+****|****+****|***MAIN ROUTINES USED BY SLDRIVER*+****|****+****|**
      module SLMOD
      use SLTSTPAK
      use SLCONSTS
      use SLPSET
      use SLUTIL
      use SAFEIO

      implicit none
      character*8 TSETNM

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      contains
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function GTSLVR(ISOLVR,SOLVER,OKEXIT)
      use SOLVRS, only: NSOLVR, SLVNAM, OKEXTS
C GTSLVR asks the user to choose the SL-Solver to be used by number in
C the range 1 to NSOLVR.
C It returns TRUE if this is done. If user aborts it returns FALSE and
C the in/out arguments are left unchanged.
C Global data used from SOLVRS (input):
C   NSOLVRS, SLVNAM, OKEXTS
C In/out arguments:
C   ISOLVER,SOLVER,OKEXIT Index within list of solver chosen; its
C           name;its 'success exit' flag value.
C     .. Scalar Arguments ..
      integer ISOLVR,OKEXIT
      character*6 SOLVER
C     .. Local Scalars ..
      integer I
C
      write(*,fmt=9999) (I,SLVNAM(I),I=1,NSOLVR)
 9999 format('Choose Solver from: ',5(i2,'=',a6,2x),
     +      (/20x,5(i2,'=',a6,2x)))
      if (.not. GETIR(ISOLVR,1,NSOLVR)) go to 100
C Global data from SOLVRS module:
      SOLVER = SLVNAM(ISOLVR)
      OKEXIT = OKEXTS(ISOLVR)
      GTSLVR = .TRUE.
      return

C Respond to user aborting input:
  100 GTSLVR = .FALSE.
      end function GTSLVR

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function GTPROB(NPROB,IPROB,TITLE,NPARM,NEPRM,PARNM,PARM,
     +                        A,B,ATYPE,BTYPE,A1,A2,B1,B2,SYM)
C GTPROB asks the user to choose the SL Problem Number IPROB in the
C range 1 to NPROB.
C It returns TRUE if this is done. If user aborts input of Problem No. or
C subsequent input of parameter values it returns FALSE and the in/out
C arguments are left unchanged.
C It sets the remaining arguments TITLE to SYM by calls to the SLTSTPAK
C routine SETUP0 and the SLDRIVER routine GTPARM (which calls the
C SLTSTPAK routine SETUP1).
C!Still a kluge!! as I haven't mastered the logic of aborting one
C state-setting when within another. Will be better when I go more
C object-oriented & only keep *one* copy of SLTSTPAK state variables.

C     .. Scalar Arguments ..
      double precision A,A1,A2,B,B1,B2
      integer IPROB,NEPRM,NPARM,NPROB
      logical SYM
      character*4 ATYPE,BTYPE
      character*72 PARNM,TITLE
C     ..
C     .. Array Arguments ..
      double precision PARM(0:10)
C     .. Local Scalars ..
      integer IPRSAV

      IPRSAV = IPROB
      write(*,advance='NO',fmt=
     +  '(/''Give Problem no. in range 1 to'',i3,'': '')') NPROB
      if (.not.GETIR(IPROB,1,NPROB)) go to 100

C Call the version of SETUP0 that doesn't change SLTSTPAK state:
      call GTDAT0(IPROB,TITLE,NPARM,NEPRM,PARNM)
      write(*,fmt='(/1x,a)') TITLE
      if (NPARM.gt.0) then
        write(*,fmt=9999,advance='NO') PARNM(1:ITRIM(PARNM))
 9999   format(/'Give values of the following parameter(s)',/,1x,a,':')
C       get entries to PARM starting with PARM(1):
        if (.not.GETRS(PARM(1),NPARM)) go to 100
      end if

C    else read in OK so:
C    Update SLTSTPAK state for new IPROB
      call SETUP0(IPROB,TITLE,NPARM,NEPRM,PARNM)
C    Update SLTSTPAK state for new PARM
      call SETUP1(PARM,A,B,ATYPE,BTYPE,A1,A2,B1,B2,SYM)
      GTPROB = .TRUE.
      return

C Respond to user aborting input:
  100 GTPROB = .FALSE.
      IPROB = IPRSAV
      call GTDAT0(IPROB,TITLE,NPARM,NEPRM,PARNM)

      end function GTPROB

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function GTPARM(NPARM,PARNM,PARM,A,B,ATYPE,BTYPE,A1,A2,B1,
     +                        B2,SYM)
C GTPARM asks the user to choose the parameter values if any for the
C current SL Problem in PARM(1:NPARM).
C It returns TRUE if this is done. If user aborts it returns FALSE and
C the in/out arguments are left unchanged.
C After the user does this it calls the SLTSTPAK routine SETUP1 to
C install the PARM values in SLTSTPAK's private memory and thus reset the
C values of A,B, ..., SYM inSLTSTPAK's memory. The reset values are returned
C through the arguments of GTPARM.
C Input arguments: NPARM, PARNM
C In/out arguments: the others
C     .. Scalar Arguments ..
      double precision A,A1,A2,B,B1,B2
      integer NPARM
      logical SYM
      character*4 ATYPE,BTYPE
      character*72 PARNM
C     ..
C     .. Array Arguments ..
      double precision PARM(0:10)
C     ..
      if (NPARM.gt.0) then
        write(*,fmt=9999,advance='NO') PARNM(1:ITRIM(PARNM))
 9999   format(/'Give values of the following parameter(s)',/,1x,a,':')
C       get entries to PARM starting with PARM(1):
        if (.not.GETRS(PARM(1),NPARM)) go to 100
      end if

C    else read in OK so:
      call SETUP1(PARM,A,B,ATYPE,BTYPE,A1,A2,B1,B2,SYM)
      GTPARM = .TRUE.
      return

C Respond to user aborting input:
  100 GTPARM = .FALSE.
      end function GTPARM

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine GTENDS(A,B,ATYPE,BTYPE,AORIG,BORIG,ATORIG,BTORIG)
C (version for drivers that handle singular ends automatically.)
C GTENDS asks the user to choose the current endpoints for the
C current SL Problem.
C It returns TRUE if this is done. If user aborts it returns FALSE and
C the in/out arguments are left unchanged.
C Input arguments:
C   AORIG,BORIG = standard endpoints
C   ATORIG,BTORIG = their type
C In/out arguments:
C   A,B   endpoints to be used for a given solver run. GTENDS checks
C         they are compatible with the standard endpoints: if AORIG is
C         singular (including of type 'WR'), A must be >= AORIG, if
C         regular no restriction is imposed. Similarly at BORIG. For
C         D02KEF which cannot handle singular endpoints automatically,
C         '>=' is replaced by '>'
C   ATYPE is set to ATORIG if A=AORIG, to 'R' if A>AORIG
C   BTYP, is set to BTORIG if B=BORIG, to 'R' if B<BORIG
C
C     .. Scalar Arguments ..
      double precision AORIG,A,BORIG,B
      character*4 ATORIG,BTORIG,ATYPE,BTYPE
C     ..
C     .. Local Scalars ..
      double precision ATRY,BTRY

C     ..
   10 write(*,ADVANCE='NO',fmt=
     +  '(/1x,''New value of a (type z to leave unchanged): '')')
      if (.not.GETR(ATRY)) ATRY = A
      write(*,ADVANCE='NO',fmt=
     +  '(1x,''New value of b (type z to leave unchanged): '')')
      if (.not.GETR(BTRY)) BTRY = B

C    Check for nontrivial interval:
      if (.not. (ATRY.lt.BTRY)) then
         write(*,ADVANCE='NO',fmt='(a)')
     +          'Supplied B is <= A, please retype: '
         go to 10
      end if

      if ((ATORIG.ne.'R'.and.ATRY.lt.AORIG) .or.
     +    (BTORIG.ne.'R'.and.BTRY.gt.BORIG)) then
         write(*,fmt=*) 'At singular end, truncated endpt must not be',
     +     'outside (A,B), please retype'
         go to 10
      end if

C     Got valid ATRY,BTRY:
      A = ATRY
      B = BTRY
      ATYPE = ATORIG
      BTYPE = BTORIG
      if (A.gt.AORIG) ATYPE = 'R'
      if (B.lt.BORIG) BTYPE = 'R'

      end subroutine GTENDS

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function GT1BC(C1,C2)
C GT1BC gets one pair of regular-type BC coefficients.
C It returns TRUE if this is done. If user aborts it returns FALSE and
C the in/out arguments are left unchanged.
C In/out arguments:
C   C1,C2  Set by user to values with the restriction: not both 0.

C     .. Scalar Arguments ..
      double precision C1,C2
C     ..
C     .. Local Arrays ..
      double precision CC(1:2)
C     ..
C    Get a pair of BC coefficients:
      if (.not.GETRS(CC,2)) go to 100
C    Check for nontrivial BC:
   10 if (CC(1).eq.0D0 .and. CC(2).eq.0D0) then
         write(*,fmt=
     +     '(''Coefficients must not both be zero. Retype: '')')
         if (.not.GETRS(CC,2)) go to 100
         go to 10

      end if
C    Success, so copy into output:
      C1 = CC(1)
      C2 = CC(2)
      GT1BC = .TRUE.
      return

C Respond to user aborting input:
  100 GT1BC = .FALSE.
      end function GT1BC

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function GTINDX(ATYPE,BTYPE,KLO,KHI,MAXEVS)
C GTINDX asks the user to give e-v index range for the current problem.
C If neither endpoint is LCO, e-vs are indexed from k=0.
C If one or both is LCO, arbitrary negative k are allowed.
C It returns TRUE if this is done. If user aborts it returns FALSE and
C the in/out arguments are left unchanged.
C Input arguments:
C   ATYPE,BTYPE
C           The type of the (current) endpoints.
C   MAXEVS  Max no. of e-vs that may be calculated at one solver call.
C In/out arguments:
C   KLO,KHI Set by user to values with the restrictions: if ATYPE or
C           BTYPE is 'LCO' then KLO is an arbitrary integer, otherwise
C           KLO>=0. KHI is >=KLO and <= KLO+MAXEVS-1.

C     .. Scalar Arguments ..
      integer KHI,KLO,MAXEVS
      character*4 ATYPE,BTYPE
C     ..
C     .. Local Arrays ..
      integer KNEW(1:2)
C     ..
   10 continue
      if (ATYPE.ne.'LCO' .and. BTYPE.ne.'LCO') then
        write(*,ADVANCE='NO',fmt=
     +    '(''Give eigenvalue index range klo,khi (0<=klo<=khi): '')')
        if (.not.GETIS(KNEW,2)) go to 100
C      Check for validity:
        if (KNEW(1).lt.0 .or. KNEW(1).gt.KNEW(2)) then
           write(*,fmt=*) 'Invalid input'
           go to 10
        end if
      else
        write(*,ADVANCE='NO',fmt=
     +  '(''Give eigenvalue index range klo,khi '',
     +  ''(klo<=khi, negative values allowed): '')')
        if (.not.GETIS(KNEW,2)) go to 100
C      Check for validity:
        if (KNEW(1).gt.KNEW(2)) then
           write(*,fmt=*) 'Invalid input'
           go to 10
        end if
      end if

      KLO = KNEW(1)
      KHI = KNEW(2)
      if (KHI-KLO+1.gt.MAXEVS) then
         KHI = KLO+MAXEVS-1
         write(*,FMT=*)'  More than',MAXEVS,
     +   ' eigenvalues requested, range reduced to ',KLO,':',KHI
      end if
      GTINDX = .TRUE.
      return

C Respond to user aborting input:
  100 GTINDX = .FALSE.
      end function GTINDX

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function GTTOL(TOLER)
C GTTOL gets a tolerance value
C It returns TRUE if this is done. If user aborts it returns FALSE and
C the in/out arguments are left unchanged.
C In/out arguments:
C   TOLER  Set by user to value with the restriction: TOLER>0.

C     .. Scalar Arguments ..
      double precision TOLER
C     ..
   10 write(*,ADVANCE='NO',
     +fmt='(/''Give tolerance (>0): '')')
      if (.not.GETR(TOLER)) go to 100
C    Check for validity:
      if (TOLER.le.0D0) then
         write(*,fmt=*) 'Invalid input'
         go to 10

      end if

      GTTOL = .TRUE.
      return

C Respond to user aborting input:
  100 GTTOL = .FALSE.
      end function GTTOL

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine CAMENU(CALSET)
C CAMENU handles the 'calculation choices submenu', which looks thus:
C-----------------------------------------------------------------------
C                 ***CALCULATION CHOICES SUBMENU***
C z Back to previous menu
C --Radio buttons controlling solver call & reporting--------
C     --Eigenvalues------------------------------------------
C [X]1   Eigenvalues only
C [ ]2   As (1) & compare with "true" values from database
C     --Eigenvalues + eigenfunctions-------------------------
C [ ]3   Eigenfunction calc using AUTO x-mesh formed by solver
C [ ]4   Eigenfunction calc using 'UNIForm' x-mesh
C      (equally spaced in a transformed variable)
C [ ]5   Eigenfunction calc using USER x-mesh from database
C [ ]6   As (5) & compare with "true" values from database
C     --Spectral density function (SDF)----------------------
C [ ]7   SDF calc using 'UNIForm' lambda-mesh
C        (equally spaced in a transformed variable)
C [ ]8   SDF calc using USER lambda-mesh from database
C [ ]9   As (8) & compare with "true" values from database
C Choose option:
C-----------------------------------------------------------------------
C It sets CALSET = (menu option number).

C     .. Parameters ..
      integer NCHOIC
      parameter (NCHOIC=9)
C     .. Scalar Arguments ..
      integer CALSET
C     .. Local Scalars ..
      integer I
C     .. Local Arrays ..
      character*3 RADIO(1:NCHOIC)
C     ..
C Arguments
C ---------
C CALSET INOUT integer. Calculation-choice flag. Initial default=1, set in
C        calling program. May be set by menu choices 1 to 9 in this
C        routine, or may be left unchanged on exit.

C If we want to have a menu loop, restore this & the GO TO below
C  100 continue
      do 105 I=1,NCHOIC
        RADIO(I) = '[ ]'
105   continue
      RADIO(CALSET) = '[X]'

      write(*,fmt=9999) RADIO
 9999 format (/17x,'***CALCULATION CHOICES SUBMENU***',/
     +' z Back to previous menu',/
     +' --Radio buttons controlling solver call & reporting--------',/
     +5x,'--Eigenvalues------------------------------------------',/
     +1x,a,'1   Eigenvalues only',/
     +1x,a,'2   As (1) & compare with "true" values from database',/
     +5x,'--Eigenvalues + eigenfunctions-------------------------',/
     +1x,a,'3   Eigenfunction calc using AUTO x-mesh formed by solver',/
     +1x,a,'4   Eigenfunction calc using ''UNIForm'' x-mesh',/
     +'      (equally spaced in a transformed variable)',/
     +1x,a,'5   Eigenfunction calc using USER x-mesh from database',/
     +1x,a,'6   As (5) & compare with "true" values from database',/
     +5x,'--Spectral density function (SDF)----------------------',/
     +1x,a,'7   SDF calc using ''UNIForm'' lambda-mesh',/
     +8x,'(equally spaced in a transformed variable)',/
     +1x,a,'8   SDF calc using USER lambda-mesh from database',/
     +1x,a,'9   As (8) & compare with "true" values from database')

      write(*,fmt=9994,advance='NO')
 9994 format(1x,'Choose option: ')
      if (.not.GETIR(CALSET,1,9)) go to 200
C      go to 100
C Respond to Quit:
  200 return

      end subroutine CAMENU

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine DBUPDT(CALSET,IPROB,KLO,KHI,FORCE)
      use DBMOD
C DBUPDT updates one or more of the database object(s) 'EVTRUE',
C 'EFTRUE', 'SDTRUE', 'XMSHUN', 'LMSHUN' according to the calculation
C choice CALSET, the 'ID' info in IPROB,KLO,KHI and the flag FORCE, as
C in this table:
C CALSET value |'ID' info used | Action taken
C     1,3      | none          | none as no stored data
C     2        | IPROB,KLO,KHI | Update 'EVTRUE'
C     4        | none          | Update 'XMSHUN'
C     5,6      | IPROB,KLO,KHI | Update 'EVTRUE','EFTRUE'
C     7        | none          | Update 'LMSHUN'
C     8,9      | IPROB         | Update 'SDTRUE'
C
C If .not.FORCE, nothing is done if the ID matches that of the stored
C data. If FORCE, a user dialogue to re-set the data (change a UNIF
C mesh or choose a different database block) is forced.
C     .. Scalar arguments ..
      integer CALSET,IPROB,KLO,KHI
      logical FORCE
c      print*,'Enter DBUPDT: database identification is'
c      call DBSUMM
      go to(10,20,30,40,50,60,70,80,90),CALSET
C or if out of range:
      write(*,*)'Panic: SLMOD%DBUPDT illegal CALSET value',CALSET
      stop
C e-vs only: no stored data
   10 if (FORCE) write(*,*)
     +  'No stored data is associated with this calculation choice'
      go to 200
C e-vs+"true"
   20 call UPEVTR(IPROB,KLO,KHI,FORCE)
      go to 200
C e-vs+e-fns AUTO x-mesh: no stored data
   30 if (FORCE) write(*,*)
     +  'No stored data is associated with this calculation choice'
      go to 200
C e-vs+e-fns UNIF x-mesh
   40 call UPXMUN(FORCE)
      go to 200
C e-vs+e-fns USER x-mesh
C e-vs+e-fns USER x-mesh+"true"
   50 continue
   60 call UPEVTR(IPROB,KLO,KHI,FORCE)
      call UPEFTR(IPROB,KLO,KHI,FORCE)
      go to 200
C SDF UNIF mesh
   70 call UPLMUN(FORCE)
      go to 200
C SDF USER lambda-mesh
C SDF USER lambda-mesh+"true"
   80 continue
   90 call UPSDTR(IPROB,FORCE)
      go to 200
  200 continue
c      print*,'Exit DBUPDT: database identification is'
c      call DBSUMM
      return
      end subroutine DBUPDT

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine DBDISP(CALSET,IPROB,KLO,KHI)
      use DBMOD
C DBDISP displays one or more of the database object(s) 'EVTRUE',
C 'EFTRUE', 'SDTRUE', 'XMSHUN', 'LMSHUN' according to the calculation
C choice CALSET, and the 'ID' info in IPROB,KLO,KHI (only a subset
C needed in some cases).
C
C     .. Scalar arguments ..
      integer CALSET,IPROB,KLO,KHI
c      print*,'Enter DBDISP: database identification is'
c      call DBSUMM
      go to(10,20,30,40,50,60,70,80,90),CALSET
C or if out of range:
      write(*,*)'Panic: SLMOD%DBDISP illegal CALSET value',CALSET
      stop
C e-vs only: no stored data
   10 continue
      go to 200
C e-vs+"true"
   20 call DIEVTR(IPROB,KLO,KHI)
      go to 200
C e-vs+e-fns AUTO x-mesh: no stored data
   30 continue
      go to 200
C e-vs+e-fns UNIF x-mesh
   40 call DIXMUN()
      go to 200
C e-vs+e-fns USER x-mesh
C e-vs+e-fns USER x-mesh+"true"
   50 continue
   60 call DIEVTR(IPROB,KLO,KHI)
      call SPAUSE
      call DIEFTR(IPROB,KLO,KHI)
      go to 200
C SDF UNIF mesh
   70 call DILMUN()
      go to 200
C SDF USER lambda-mesh
C SDF USER lambda-mesh+"true"
   80 continue
   90 call DISDTR(IPROB)
      go to 200
  200 continue
c      print*,'Exit DBDISP: database identification is'
c      call DBSUMM
      return
      end subroutine DBDISP

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine PRHELP(FNAME)
C PRHELP displays the 'help text' held in file FNAME.hlp (i.e. the name
C held in FNAME, with .hlp appended). FNAME can be a full or relative
C pathname.
C     .. Scalar Arguments ..
      character*(*) FNAME
C     ..
C     .. Local Scalars ..
      integer IOFLAG,LFNAME
      character*75 LINE
      character*79 STARS
C     ..

      STARS = '****************************************'
     +      //'***************************************'
      LFNAME = ITRIM(FNAME)
      open(20,STATUS='OLD',FILE=FNAME(1:LFNAME)//'.hlp',ERR=190,
     +     IOSTAT=IOFLAG)

      write(*,*)
      LINE=STARS
      LINE(25:25+18+LFNAME) = 'Help from file '//FNAME(1:LFNAME)//'.hlp'
  100 write(*,fmt=110) '* '//LINE//' *'
  110 format(a79)
      read (20,fmt='(a)',END=200) LINE
      go to 100

  190 write(*,fmt='(3a,i3)') 'Couldn''t open file ',
     +  FNAME(1:LFNAME)//'.hlp',', error #',IOFLAG
  200 write(*,fmt=110) STARS
      close (20)
      end subroutine PRHELP

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine REPORT(QEVREP,QEFREP,QSDREP,
     *   SOLVER,IPROB,KLO,KHI,TOLER,IFAIL,
     *   EV,
     *   NXMESH,XMESH,EF,PDEF,
     *   NLMESH,LMESH,RHO,
     *   ELAPSE,OKEXIT)

      use DBMOD
C REPORT is called immediately after SOLVIT to display the computed
C e-vs, e-fns and SDF; and if asked, compare them with "true" values.
C It does not check the validity of what is asked to do, e.g. if asked
C to print e-fn values when none have been computed it will probably
C print garbage.
C All needed "true" values must have been read into the database cache
C 'EVTRUE', 'EFTRUE' or 'SDTRUE' from the database BEFORE REPORT is
C called.
C
C Output:
C Main display is to screen.
C Unit 20 is logfile solver//nn//'.m' where solver is 6 character name
C of current solver, nn is current Problem no.
C Unit 21 is auxiliary logfile solver//nn//'.aux' holding eigenfunctions
C in form suitable to include in EFTRU.nn database files.
C Unit 21 also receives logging data from SLEIGN2 and SLEDGE, if this
C has been switched on (by changing the source code of SOLVIT).

C Arguments:
C QEVREP,QEFREP,QSDREP
C            Integer, describe what reporting of evs, e-fns and SDF
C            (respectively) the user requested:
C            =0 Don't report
C            =1 Report
C            =2 Report and compare with "true"  values from database
C            In case of comparing with "true" values, whether this is
C            done depends on whether the "true" data is actually
C            present. So a call to the DBMOD routines EVCHK, EFCHK or
C            SDCHK is made at the appropriate place, to check this.
C SOLVER     Name of Solver, 'SLEDGE','SL02F' etc
C IPROB      Problem number
C KLO,KHI    Eigenvalue index range last used by SOLVIT
C TOLER      Tolerance last used by SOLVIT
C IFAIL      IFAIL(K) holds error flag for calculation of e-v of index K
C EV         EV(K,1) holds computed e-v of index K
C            EV(K,2) holds estimated error/error-bound for this (see
C            info in the various SOLVIT routines)
C            If IFAIL(K).ne.OKEXIT these will be suspect or absent.
C            It is assumed SOLVIT has set sensible values in this case!
C NXMESH     NXMESH=no. of points in x-mesh where e-fn values computed,
C            including endpoints A & B.
C XMESH      x-mesh is in XMESH(1:NXMESH)
C EF,PDEF    (only used when QEFREP is 1 or 2)
C            Indexed 1:MAXMSH,KLO:KHI. EF(I,K),PDEF(I,K) hold value of
C            u(x),pu'(x) at x=XMESH(I), for K-th e-fn, for I=1:NXMESH
C            If IFAIL(K).ne.OKEXIT these will be suspect or absent.
C  !!NOTE!!  It is assumed SOLVIT has set them to all zero if absent.
C NLMESH     NLMESH=no. of points in lambda-mesh where SDF values
C            computed, including endpoints A & B.
C LMESH      x-mesh is in LMESH(1:NLMESH)
C RHO        (only used when QSDREP is 1 or 2)
C            RHO(I) holds value rho(lambda) of SDF at lambda=LMESH(I)
C ELAPSE     Total CPU time for Solver call
C OKEXIT     Value of Solver's error indicator on SUCCESSFUL exit
C
C Imported from DBMOD:
C EVTRU      (only used when QEVREP.eq.2)
C            Array of "true" eigenvalues in given range KLO..KHI,
C            is in EVTRU(1,K). Believed bound on error is in EVTRU(2,K)
C QEVTRU     (only used when QEVREP.eq.2)
C            Flags to say which "true" evs currently provided. If K-th e-v
C            is present it is in EVTRU(1,K), and QEVTRU(K) equals PRESNT;
C            otherwise QEVTRU(K) equals ABSENT.
C EFTR,PDEFTR (only when QEFREP.eq.2)
C            Array of "true" e-functions & derivatives in given range
C            KLO..KHI on given mesh, if provided.
C QEFTRU     (only used when QEFREP.eq.2)
C            Flags to say which "true" e-functions currently provided.
C            If K-th e-fn u_k(x) is present then
C               values of u_k(x(i))   are in EFTR(i,K)
C               values of pu'_k(x(i)) are in PDEFTR(i,K)
C            (i=0..NUMX+1), and QEFTRU(K) equals PRESNT;
C            otherwise QEFTRU(K) equals ABSENT.
C SDTRU      (only when QSDREP is 2)
C            array of "true" rho(lambda) values on lambda-mesh
C NLMESH     no. of points in lambda-mesh where SDF computed.
C LMESH      lambda-mesh is in LMESH(1:NLMESH).

C! QERROR is not handled properly! It works as the program stands
C! but if both evs/efns and SDF were computed in one 'solve', error
C! conditions in evs would be overwritten by success in SDF.
C     ..
C     .. Scalar Arguments ..
      integer IPROB,KLO,KHI,NLMESH,NXMESH,OKEXIT,
     +  QEVREP,QEFREP,QSDREP
      double precision ELAPSE,TOLER
      character*(*) SOLVER
C     ..
C     .. Array Arguments ..
      integer IFAIL(KLO:KHI)
      double precision EV(KLO:KHI,1:2),
     +  XMESH(1:MAXMSH),
     +  EF(1:MAXMSH,KLO:KHI),PDEF(1:MAXMSH,KLO:KHI),
     +  LMESH(1:MAXMSH),RHO(1:MAXMSH)
C     ..
C     .. Local Scalars ..
      integer IK,J,K,NFEVAL,KOFSET
      double precision EFER,PDEFER,ABSERR
      logical QERROR
      character KSTR*12
C     ..
C
C***+****|****+****|****BLOCK 1: EIGENVALUES+****|****+****|****+****|**
        if (QEVREP.gt.0) then
        KOFSET = KLO-1
C       MATLAB initializing statement for e-v variable in .m file:
        write(20,*) 'ev=[];'
        write(*,fmt=300)
300     format('   K',T10,'Eigenvalue',T35,'IFAIL',
     +                                             T52,'Absolute error',
     +         /    T30,'(if error exit)',T50,'(est.)    ("true")')
C       QERROR is set to = any(IFAIL(KLO:KHI).ne.OKEXIT) for PRHELP
C       message below
        QERROR = .FALSE.

C       Loop over eigenvalues:
        do 310 K = KLO,KHI
          if (IFAIL(K).eq.OKEXIT) then
            write(*,advance='NO',fmt='(i4,3x,1p,g22.15,t47,d10.2)')
     +        K,EV(K,1),EV(K,2)
          else
            QERROR = .TRUE.
            write(*,advance='NO',
     +        fmt='(i4,''(?)'',1pg22.15,t35,i4,t47,d10.2)')
     +        K,EV(K,1),IFAIL(K),EV(K,2)
          end if

C         Are e-v errors to be displayed?
C         Note EVTRU(2,K) holds believed bound on abs error in
C         EVTRU(1,K), which we may use in the future.
C         Note output on unit 20 puts out a dummy value for ABSERR if no
C         genuine value present
          if (QEVREP.ge.2 .and. EVCHK(IPROB,KLO,KHI) .and.
     +                QEVTRU(K-KOFSET).eq.PRESNT) then
            ABSERR = (EV(K,1)-EVTRU(1,K-KOFSET))
            write(*,fmt='(1p,d10.2)') ABSERR
          else
            ABSERR = 1.0d20
            write(*,*)
          end if
          write(20,fmt=305) SOLVER,IPROB,TOLER,K,IFAIL(K),EV(K,1),ABSERR
          write(21,fmt=305) SOLVER,IPROB,TOLER,K,IFAIL(K),EV(K,1),ABSERR
  305     format (1x,'%',a,'/0',i3,1p,d10.2,i5,i4,d24.15,d10.2)
          write(20,fmt=306) K,EV(K,1),ABSERR
  306     format(1x,'ev=[ev;',i5,1p,d24.15,d10.2,'];')
          write(21,'(a,i3)') 'U ',K

C***+****|****+****|****BLOCK 2: EIGENFUNCTIONS (inside block 1)+****|**
C         Are e-fn values and errors to be displayed?
C         If error flag raised for K-th e-v (if IFAIL(K).eq.0) user is
C         told values are suspect but they are always displayed.
C         "True" e-fns: are available in EFTR,PDEFTR but in general only
C         for some K in KLO..KHI, namely those for which QEFTRU(K) .eq.
C         PRESNT.
C         Errors in computed e-fns are reported for these K if QEFREP=2

          if (QEFREP.ge.1) then

            call SPAUSE

            if (IFAIL(K).ne.OKEXIT) write(*,*)
     +        ' Eigenfunction values are suspect!'
            write(*,fmt=9999,advance='NO') K
 9999       format (1x,'Eigenfunction Values for k=',i4,/
     +              1x,8x,'x         u(x)              pu''(x)')
C           Error can be displayed if e-fn values present & mesh=USER
            if (QEFREP.ge.2 .and. EFCHK(IPROB,KLO,KHI) .and.
     +                 QEFTRU(K-KOFSET).eq.PRESNT) then
              write(*,fmt=9998)
 9998         format (11x,'u(x) error  pu''(x) error')
            else
              write(*,fmt=*)
            end if

C           Display e-fn values ..
C           To logfile, put start of a MATLAB statement that will
C           form a matrix from the x,u and pu' values:
            write(20,fmt=*) 'a=[ ...'
            do 307 J=1,NXMESH
              write(*,fmt=9997,advance='NO')XMESH(J),EF(J,K),PDEF(J,K)
              write(20,fmt=9997,advance='NO')
     +          XMESH(J),EF(J,K),PDEF(J,K)
 9997         format (1x,0p,f12.7,1p,2e18.10)
              write(21,fmt=99970) EF(J,K),PDEF(J,K)
99970         format (1x,1p,2e18.10)
C             .. and e-fn errors if asked for & available
              if (QEFREP.ge.2 .and. EFCHK(IPROB,KLO,KHI) .and.
     +                  QEFTRU(K-KOFSET).eq.PRESNT) then
                EFER = EF(J,K) - EFTR(J,K-KOFSET)
                PDEFER = PDEF(J,K) - PDEFTR(J,K-KOFSET)
                write(*,fmt=9996) EFER,PDEFER
                write(20,fmt=9996) EFER,PDEFER
 9996           format (2x,1p,2e12.2)
              else
                write(*,fmt=*)
                write(20,fmt=*)
              end if
  307       continue
C           Put out MATLAB statements that extract e-fns and their
C           derivatives to variables u0,pdu0,u1,pdu1, ...
            call INT2ST(K,KSTR,IK)
            write(20,fmt=9995) KSTR(1:IK),KSTR(1:IK)
 9995       format(' ];',/,
     +             ' u',a,'=a(:,2);',/,
     +             ' pdu',a,'=a(:,3);')
            if (K.eq.KLO) write(20,fmt=9994)
 9994       format(' x=a(:,1);')
            call SPAUSE
          end if
  310   continue

      end if

C***+****|****+****|*BLOCK 3: SPECTRAL DENSITY FUNCTION****|****+****|**
      if (QSDREP.ge.1) then
        QERROR = IFAIL(KLO).ne.OKEXIT
        if (QERROR) write(*,*)
     +    '***Exit flag has failure value: SDF values are suspect!'
        write(*,fmt=9993,advance='NO') IFAIL(KLO)
 9993   format (1x,'Spectral Density Values: exit flag=',i3,/
     +          9x,'lambda',18x,'rho(lambda)')
C       Error can be displayed if true SDF present
        if (QSDREP.ge.2 .and. SDCHK(IPROB)) then
          write(*,fmt=9992)
 9992     format (11x,'rho(lambda) error')
        else
          write(*,fmt=*)
        end if
C       To logfile, put start of a MATLAB statement that will
C       form a matrix from the lambda and rho(lambda) values:
        write(20,fmt=*) 'sdf=[ ...'
        do 312 J=1,NLMESH
          write(*,advance='NO',fmt='(i4,1p,2g25.15)')
     +      J,LMESH(J),RHO(J)
          write(20,advance='NO',fmt='(1p,2g25.15)')
     +      LMESH(J),RHO(J)
          if (QSDREP.ge.2 .and. SDCHK(IPROB)) then
            write(*,fmt='(1p,d14.2)') RHO(J)-SDTRU(1,J)
            write(20,fmt='(1p,d10.2)') RHO(J)-SDTRU(1,J)
          else
            write(*,*)
            write(20,*)
          end if
  312   continue
        write(20,fmt=9991)
 9991   format(' ];',/,
     +         ' lambda=sdf(:,1);',/,
     +         ' rho=sdf(:,2);')

        call SPAUSE

      end if

C     Get no. of function evaluations used in SOLVIT call
      NFEVAL = NEVAL(0)

      write(*,fmt=
     + '(1X,"CPU secs:",f10.4,", No. of function evaluations: ",i10)'
     +  ) ELAPSE,NFEVAL

      write(20,fmt=315) SOLVER,IPROB,TOLER,KLO,KHI,ELAPSE,NFEVAL
      write(21,fmt=315) SOLVER,IPROB,TOLER,KLO,KHI,ELAPSE,NFEVAL
315   format (1x,'%',a,'/1',i3,1p,d10.2,2i5,0p,f10.4,i10)

      if (QERROR) call PRHELP('errflags')

      end subroutine REPORT

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine OPNLOG(TSETNM,IPROB,SOLVER)
      integer IPROB
      character TSETNM*8,SOLVER*6
      integer IOERR,IS,IFNAM0,IFNAME
      character STR*10,FNAME0*43,FNAME1*43

C     OPNLOG opens logfiles for new Problem. Note 'append' status in
C     case same Problem is used several times in a run.
C     Main logfile name, in FNAME, is 'solvernn.m' where
C       solver = name held in SOLVER variable
C       nn = decimal representation of IPROB
C     Log file is intended as an executable MATLAB file, hence the '.m'
C     Auxiliary logfile is 'solvernn.aux'

C     This is a f77-ish trick for getting nn with leading zeros
      write(STR,'(i3)') 100+IPROB
      call STTRIM(SOLVER,IS)
      FNAME0 = ADDDIR(TSETNM, SOLVER(1:IS)//STR(2:3))
      call STTRIM(FNAME0,IFNAM0)

      FNAME1 = FNAME0(1:IFNAM0)//'.m'
      call STTRIM(FNAME1,IFNAME)
      close(20)
      open(20,file=FNAME1, status='OLD', action='WRITE',
     + position='APPEND', iostat=IOERR)
      if (IOERR.eq.0) then
        write(*,*) 'Appending to log-file ',FNAME1(1:IFNAME)
      else
        open(20,file=FNAME1, status='NEW', action='WRITE',
     +    iostat=IOERR)
        if (IOERR.eq.0) then
          write(*,*) 'Opening new log-file ',FNAME1(1:IFNAME)
        else
          write(*,*) 'Unable to open log-file unit 20, file ',
     +     FNAME1(1:IFNAME),', IOSTAT=',IOERR
          stop
        end if
      end if
C ..and auxiliary logfile intended to report e-fns suitably for
C including in database files
      FNAME1 = FNAME0(1:IFNAM0)//'.aux'
      call STTRIM(FNAME1,IFNAME)
      close(21)
      open(21,file=FNAME1, status='OLD', action='WRITE',
     + position='APPEND', iostat=IOERR)
      if (IOERR.eq.0) then
        write(*,*) 'Appending to log-file ',FNAME1(1:IFNAME)
      else
        open(21,file=FNAME1, status='NEW', action='WRITE',
     +    iostat=IOERR)
        if (IOERR.eq.0) then
          write(*,*) 'Opening new log-file ',FNAME1(1:IFNAME)
        else
          write(*,*) 'Unable to open log-file unit 21, file ',
     +     FNAME1(1:IFNAME),', IOSTAT=',IOERR
          stop
        end if
      end if
      write(*,*)
      end subroutine OPNLOG

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      character*80 function ADDDIR(DIR,PATH)
      character*(*) DIR,PATH
C DIR is a directory name, PATH is a relative path name (possibly
C just a file name). ADDDIR returns the result of prepending DIR
C to PATH following the path-naming conventions of the computer system
C described by OPSYS. The allowed values of OPSYS are currently
C       'MSDOS', 'UNIX', 'VMS' or MAC'.
C
C Trailing spaces on DIR are checked for and removed.
C Total length of result mustn't exceed 80 chars.
C
C On first call, the user is asked to specify OPSYS.
C! If user aborts this, ADDDIR returns all spaces as result: klugey way
C! to abort whole program in this case. IMPROVE!!
C
C To avoid the OPSYS enquiry, the implementation on a specific platform
C can simply initialize OPSYS to the appropriate value instead of UNDEF.

      integer UNDEF,MSDOS,UNIX,MAC,VMS
      parameter (UNDEF=0,MSDOS=1,UNIX=2,MAC=3,VMS=4)
      integer LENDIR
      character*80 A
C Backslash: a literal one is treated as escape char by some compilers
      character*1,parameter:: BSLASH=ACHAR(92)
      integer,save:: OPSYS=UNDEF

      if (OPSYS .eq. UNDEF) then
        write(*,'(a)',advance='NO')
     +   'Please describe path-name conventions of your system.',
     +   '1 MSDOS e.g. DIR1'//BSLASH//'DIR2'//BSLASH//'FILE',
     +   '2 Unix  e.g. dir1/dir2/file',
     +   '3 Mac   e.g. :dir1:dir2:file',
     +   '4 VMS   e.g. [DIR1.DIR2]FILE',
     +   'Choose option: '
        if (.not. GETIR(OPSYS,1,4)) OPSYS = UNDEF
      end if

      LENDIR = ITRIM(DIR)
      if (OPSYS.eq.MSDOS) then
        ADDDIR = DIR(1:LENDIR)//BSLASH//PATH
      else if (OPSYS.eq.UNIX) then
        ADDDIR = DIR(1:LENDIR)//'/'//PATH
      else if (OPSYS.eq.MAC) then
        ADDDIR = DIR(1:LENDIR)//':'//PATH
      else if (OPSYS.eq.VMS) then
        A = PATH
        if (A(1:1).ne.'[') A='[]'//A
        ADDDIR = '[.'//DIR(1:LENDIR)//A(2:LEN(A))
      else if (OPSYS.eq.UNDEF) then
        ADDDIR = ' '
      else
        write(*,*) 'ADDDIR: Name of operating system "',OPSYS,
     +  '" unknown - unable to construct filenames'
        stop
      end if
      end function ADDDIR

      subroutine MBROWS(TSETNM,NPROB)
      character TSETNM*8,BROWSF*12,PARNM*72,TITLE*72
      integer IOERR,IPROB,IT,NEPRM,NPARM,NPROB
      call STTRIM(TSETNM,IT)
      BROWSF = TSETNM(1:IT)//'.lst'
      open(21,file=BROWSF,status='OLD',iostat=IOERR)
      if (IOERR.ne.0) then
        open(21,file=BROWSF,status='NEW',iostat=IOERR)
        if (IOERR.ne.0) then
          write(*,*) 'Unable to open file ',BROWSF,', IOSTAT=',IOERR
        else
          call PRHELP('slbrows')
          call SPAUSE
          write(21,*) 'Contents of Test Set ',TSETNM
          do 50 IPROB=1,NPROB
            call SETUP0(IPROB,TITLE,NPARM,NEPRM,PARNM)
            write(21,fmt='(/1x,''No.'',i2,'':'',a)') IPROB,TITLE
            if (NPARM.le.0) then
              write(21,fmt='(1x,''No parameters'')')
            else
              write(21,fmt='(1x,''Parms '',a)') PARNM
            end if
   50     continue
          close(21)
          write(*,*) 'Created Browse-File ',BROWSF
          write(*,*) 'Now look at it with your favourite editor'
        end if
      end if
      end subroutine MBROWS

      end module SLMOD

C***+****|****+****|****+*THE MAIN PROGRAM**+****|****+****|****+****|**
      program SLDRIVER

      use SLMOD !also imports SLCONSTS, SLPSET, SLUTIL, SAFEIO, SLTSTPAK
      use DBMOD !also imports SLCONSTS, SLPSET, SLUTIL, SAFEIO
      use SOLVRS, only: DEPREC,FORBID,NSOLVR,SLVNAM,SOLVIT

      implicit none
C
C     .. Parameters ..
      character*(*) EOFMSG
      parameter (EOFMSG='z or Z followed by ENTER')
C   Used to translate 'calculation-mode' values to SOLVIT & REPORT flags:
      integer,parameter:: QCALC (1:9)=(/1,1,2,3,3,3,4,4,4/)
     +                   ,QEVREP(1:9)=(/1,2,1,1,1,2,0,0,0/)
     +                   ,QEFREP(1:9)=(/0,0,1,1,1,2,0,0,0/)
     +                   ,QSDREP(1:9)=(/0,0,0,0,0,0,1,1,2/)
      character*24,parameter::ACMENU(1:9)=
     +                  (/'                EVs only'
     +                   ,'              EVs+"true"'
     +                   ,'       EV+EFn, AUTO mesh'
     +                   ,'       EV+EFn, UNIF mesh'
     +                   ,'       EV+EFn, USER mesh'
     +                   ,'EV+EFn, USER mesh+"true"'
     +                   ,'          SDF, UNIF mesh'
     +                   ,'          SDF, USER mesh'
     +                   ,'   SDF, USER mesh+"true"'/)
C     ..
C     .. Local Scalars ..
      double precision ELAPSE,TOLER
      integer CHOICE,FLAG,I,IDUM,IPARM,IPROB,ISOLVR,IT,K,KHI,KLO,
     +        NEPRM,NLMESH,NXMESH,NPARM,NPROB,OKEXIT,CALSET
     +       ,CHOIC9
      logical OK,SYM
      character ATYPE*4,BTYPE*4,ATORIG*4,BTORIG*4,BCSTRA*18,BCSTRB*18,
     +  PARNM*72,TITLE*72,SOLVER*6
C     ..
C     .. Local Arrays ..
      integer IFAIL(1:MAXEVS)
      double precision PARM(0:10),
     +   EV(1:2*MAXEVS),
     +   XMESH(1:MAXMSH),
     +   EF(1:MAXMSH,1:MAXEVS),PDEF(1:MAXMSH,1:MAXEVS),
     +   LMESH(1:MAXMSH),
     +   RHO(1:MAXMSH)

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C Global Initializations ...
C Initialize SLTSTPAK package:
      call TSTINI(TSETNM,NPROB)
C Initialize DBMOD package:
      call DBINIT
C Initialize Playback facility
      call OPNPBK
C Intro & Help message:
      write(*,fmt=19999) TSETNM,NPROB,NSOLVR,(SLVNAM(I),I=1,NSOLVR)
19999 format(/,79('*'),/,
     + '*',t79,'*',/,
     + '*',t24,'SLDRIVER Version 4.1 June 1998',t79,'*',/,
     + '*',t16,'John Pryce''s Driver for Sturm-Liouville solvers',
     +   t79,'*',/,
     + '*',t79,'*',/,
     + '*',t15,'Running with Problem Set "',a8,'" of',i3,' problems',
     +   t79,'*',/,
     + '*',t79,'*',/,
     + '*',t24,'and the following',i3,' Solvers:',t79,'*',/,
     + ('*',t79,'*',t22,5a8))

      write(*,fmt=19998) EOFMSG
19998 format('*',t79,'*',/,
     + 79('*'),//,
     + ' NOTE: Typing ',a,' aborts any data-input operation')

      call PRHELP('initpuff')
      call SPAUSE

      call PRHELP('slhelp0')
      call SPAUSE

C Check for existence of browse file, create it if not:
      call MBROWS(TSETNM,NPROB)

C Initial setting of eigenvalue & eigenfunction & SDF flags:
C     Calculation mode is 'eigenvalues only'
      CALSET = 1
C     No points in x-mesh for e-fns or lambda-mesh for SDF
      NXMESH = 0
      NLMESH = 0

C Initial input of Problem:
C   Quit/no-quit loop for initial input starts here:
   10 continue

C     Initialize TRUEDB ("true" e-v & e-fn database) package.
C     Input is the (relative) Pathname where database is to be found.
C     ADDDIR will ask user to specify Op Sys pathname conventions here
        call STTRIM(TSETNM,IT)
C!    WARNING! unsafe access of DBMOD variable:
        DBPATH = ADDDIR('truevals',' ')
C!    ' ' means abort,see doc of ADDDIR:
        if (DBPATH.eq.' ') go to 90
        DBPATH = ADDDIR(TSETNM, DBPATH)

      if (.not. GTSLVR(ISOLVR,SOLVER,OKEXIT)) go to 90

      if (.not.GTPROB(NPROB,IPROB,TITLE,NPARM,NEPRM,PARNM,PARM,
     +    AORIG,BORIG,ATORIG,BTORIG,A1ORIG,A2ORIG,B1ORIG,B2ORIG,SYM))
     +    go to 90

C     Initialize 'endpoints used' to the defaults for this Problem
      A = AORIG
      B = BORIG
C     Also the 'types of the endpoints used'
      ATYPE = ATORIG
      BTYPE = BTORIG
C     and the BC coefficients
      A1 = A1ORIG
      A2 = A2ORIG
      B1 = B1ORIG
      B2 = B2ORIG

C     Set initial e-v index-range and tolerance
C      if (.not.GTINDX(ATYPE,BTYPE,KLO,KHI,MAXEVS)) go to 90
C      if (.not.GTTOL(TOLER)) go to 90
C     Set default initial e-v index-range and tolerance
      KLO = 0
      KHI = 0
      TOLER = 1d-4

      go to 99
C Handle user-aborts:
   90 if (YESNO('Do you want to quit? (y/n) ').eq.'n') go to 10
      call CLOPBK
      stop

   99 continue

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C                      *** Main Menu loop ***
  100 continue
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      write(*,fmt=9999) TSETNM, ACMENU(CALSET), SOLVER, IPROB, TITLE,
     +            PARNM,(PARM(IPARM),IPARM=1,NPARM)

 9999 format (//,'Testset: ',a,t33,'*** MAIN MENU ***',
     +       t66,'Current values',/
     +       ' z Quit (terminate run)',/
     +       ' 0 Calculation choices',33x,a24,/
     +       '   eigenvalues; eigenvalues + eigenfunctions;',
     +        ' spectral density function',/
     +       ' 1 Choose solver',57x,a6,/
     +       ' 2 Choose problem',59x,i3,/
     +       '     (',a72,')',:,/
     +       ' 3 Give parameters (if any) for current problem',/
     +       3x,a72,:,/
     +       1p,3x,'Currently',7x,4d15.8,(:,19x,4d15.8))
      if (ATYPE.eq.'R' .or. ATYPE.eq.'WR') then
        BCSTRA = " A1*u(a)=A2*pu'(a)"
      elseif (ATYPE(1:2).eq.'LC') then
        BCSTRA = '[u,A1*f+A2*g](a)=0'
      else
        BCSTRA = '                  '
      end if
      if (BTYPE.eq.'R' .or. BTYPE.eq.'WR') then
        BCSTRB = " B1*u(b)=B2*pu'(b)"
      elseif (BTYPE(1:2).eq.'LC') then
        BCSTRB = '[u,B1*f+B2*g](b)=0'
      else
        BCSTRB = '                  '
      end if
      write(*,fmt=9998) A,B,A1,B1,BCSTRA,BCSTRB,A2,B2,KLO,KHI,TOLER

 9998 format(' 4 Give endpoints for current problem',4x,1p,' A =',
     +       d15.8,' B =',d15.8,/
     +       ' 5 Give one/both BCs for current problem','  A1=',d15.8,
     +       ' B1=',d15.8,/
     +       '  (',a18,',',a18,') A2=',d15.8,' B2=',d15.8,/
     +       ' 6 Give eigenvalue index range',35x,i5,'  to',i5,/
     +       ' 7 Give tolerance',52x,d10.2,/
     +       ' 8 Solve problem as currently given',/
     +       ' 9 Display/change mesh & "true" values from database')

      write(*,fmt=9996,advance='NO') AORIG,BORIG,ATORIG,BTORIG,EOFMSG
 9996 format ('(Extra information:',4x,'Original endpoints',1p,
     +       ' A =',d15.8,' B =',d15.8,/30x,'Their types',12x,a4,15x,
     +       a4,'  )',//'To abort any data-input operation type ',a,
     +       /'Choose option: ')

      if (GETIR(CHOICE,0,9)) go to (200,210,220,230,240,250,260,
     +                               270,280,290) CHOICE + 1
C ..else drop through into Quit option..

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C                    *** Main Menu choices ***
C
C In each choice, if the flag OK is returned as .false. this signals
C user aborted the option, nothing was changed & we return to main menu
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**

C     CHOICE=Z, Quit
      if (YESNO('Do you want to quit? (y/n) ').eq.'n') go to 1001
C     Close the Playback file
      call CLOPBK
      stop

C     CHOICE=0, Calculation choices
  200 call CAMENU(CALSET)
      go to 1001

C     CHOICE=1, Input Solver
  210 if (.not. GTSLVR(ISOLVR,SOLVER,OKEXIT)) go to 1001
      write(*,*) '  Solver is now ',SOLVER
      go to 1000

C     CHOICE=2, Input Problem number
  220 if (GTPROB(NPROB,IPROB,TITLE,NPARM,NEPRM,PARNM,PARM,AORIG,BORIG,
     +     ATORIG,BTORIG,A1ORIG,A2ORIG,B1ORIG,B2ORIG,SYM)) then
C       Initialize 'endpoints used' to the defaults for this Problem
        A = AORIG
        B = BORIG
C       Also the 'types of the endpoints used'
        ATYPE = ATORIG
        BTYPE = BTORIG
C       and the BC coefficients
        A1 = A1ORIG
        A2 = A2ORIG
        B1 = B1ORIG
        B2 = B2ORIG
      end if

      go to 1001

C   CHOICE=3, Input Parameters
  230 if (GTPARM(NPARM,PARNM,PARM,AORIG,BORIG,ATORIG,BTORIG,
     +        A1ORIG,A2ORIG,B1ORIG,B2ORIG,SYM)) then
C     This resets the 'ORIG' endpts & BCs to the defaults for this value of
C     PARMs. If A is now <=AORIG and the latter is not Regular, set it
C     =AORIG and re-set its properties to the defaults, else set its type
C     to Regular and leave it and its BCs as is.
C     Similarly at B.
        if (A.le.AORIG .and. ATORIG.ne.'R') then
          A = AORIG
          ATYPE = ATORIG
          A1 = A1ORIG
          A2 = A2ORIG
        else
          ATYPE = 'R'
        end if
        if (B.ge.BORIG .and. BTORIG.ne.'R') then
          B = BORIG
          BTYPE = BTORIG
          B1 = B1ORIG
          B2 = B2ORIG
        else
          BTYPE = 'R'
        end if
      end if

      go to 1000

C   CHOICE=4, new values of one or both endpts
  240 call GTENDS(A,B,ATYPE,BTYPE,AORIG,BORIG,ATORIG,BTORIG)
C    A,B are NOT recorded in SLTSTPAK private memory.
      go to 1000

C   CHOICE=5, new values of one or both BCs
C     Left BC:
  250 write(*,advance='NO',fmt=
     +  '(/1x,''New A1,A2 (z or Z to leave unchanged): '')')
      OK = GT1BC(A1,A2)
C    If changed, copy the new values into SLTSTPAK internal memory:
      if (OK) call SETBCS(0,A1,A2)
C     Right BC:
      write(*,advance='NO',fmt=
     +  '(1x,''New B1,B2 (z or Z to leave unchanged): '')')
      OK = GT1BC(B1,B2)
C    If changed, copy the new values into SLTSTPAK internal memory:
      if (OK) call SETBCS(1,B1,B2)
      go to 1000

C   CHOICE=6, new eigenvalue index range
  260 OK = GTINDX(ATYPE,BTYPE,KLO,KHI,MAXEVS)
      go to 1000

C   CHOICE=7, new tolerance
  270 OK = GTTOL(TOLER)
      go to 1000

C   CHOICE=8, solve problem as currently given
  280 continue
C     Open logfile(s) for current Solver & Problem either 'new' or
C     'old,append'
      call OPNLOG(TSETNM,IPROB,SOLVER)

C     Update the cached "true" values and mesh from database
      call DBUPDT(CALSET,IPROB,KLO,KHI,FORCE=NO)

C     If an x-mesh is required, copy from 'XMSHUN' or 'EFTRUE'
      if (CALSET.eq.4) then
        if (.not.GTXMUN(NXMESH,XMESH)) then
          write(*,fmt=9289)
 9289     format('Sorry: UNIForm x-mesh, needed for this ',
     +    'calculation choice, isn''t set.')
          go to 288
        end if
      end if
      if (CALSET.eq.5 .or. CALSET.eq.6) then
        if (.not.GTXMTR(NXMESH,XMESH)) then
          write(*,fmt=9288)
 9288     format('Sorry: USER x-mesh from file, needed for this ',
     +    'calculation choice, isn''t set',/
     +    6x,'..probably because no data file exists for this problem.')
          go to 288
        end if
      end if

C     If a lambda-mesh is required, copy from 'LMSHUN' or 'EFTRUE'
      if (CALSET.eq.7) then
        if (.not.GTLMUN(NLMESH,LMESH)) then
          write(*,fmt=9287)
 9287     format('Sorry: UNIForm lambda-mesh, needed for this ',
     +    'calculation choice, isn''t set.')
          go to 288
        end if
      end if
      if (CALSET.eq.8 .or. CALSET.eq.9) then
        if (.not.GTLMTR(NLMESH,LMESH)) then
          write(*,fmt=9286)
 9286     format('Sorry: USER lambda-mesh from file, needed for ',
     +    'this calc. choice, isn''t set',/
     +    6x,'..probably because no data file exists for this problem.')
          go to 288
        end if
      end if

C     Use arrays QCALC etc. to convert different calculation-mode values
C       CALSET = 1 2|3 4 5 6|7 8 9
C     to values needed by SOLVIT:
C       QCALC  = 1 1|2 3 3 3|4 4 4
C     and by REPORT:
C       QEVREP = 1 2|1 1 1 2|0 0 0
C       QEFREP = 0 0|1 1 1 2|0 0 0
C       QSDREP = 0 0|0 0 0 0|1 1 2

C     Clear function evaluation counter
      IDUM = NEVAL(0)

C     Ensure that ev-index range is valid in case user switched from
C     a problem with an LCO end to one without:
      if (KLO.lt.0 .and. ATYPE.ne.'LCO' .and. BTYPE.ne.'LCO') then
        write(*,
     +    '(/"LCO-style eigenvalue index range no longer valid...")')
        OK = GTINDX(ATYPE,BTYPE,KLO,KHI,MAXEVS)
      end if

C     Set all evs in range to easily recognized dummy value in case
C     severe failure leaves them without even an estimate
      do 282 K=KLO,KHI
        EV(K-KLO+1) = XINFTY
 282  continue

      write(*,fmt=9284)
 9284 format(/,15('*'),
     +' CALLING SOLVER, SOME DIAGNOSTICS MAY BE PRINTED ',15('*'))

      call SOLVIT(ISOLVR,QCALC(CALSET),A,B,A1,A2,B1,B2,ATYPE,BTYPE,
     +     KLO,KHI,TOLER,
     +     NXMESH,XMESH,NLMESH,LMESH,
     +     IFAIL,EV,EF,PDEF,RHO,
     +     ELAPSE,FLAG)

C FLAG=OKCALC means solver at least had a go at requested computation.
C FLAG=DEPREC means it wasn't recommended but user made it try anyway.
C FLAG=FORBID means it couldn't even try (see comments in SOLVIT):
      if (FLAG.ne.FORBID) then
C       Now we know REPORT will be called, output 'About the Problem'
C       info (which main program knows but REPORT doesn't):
      write(*,fmt=9283) TSETNM, IPROB, SOLVER, ACMENU(CALSET),
     +            PARM(1:NPARM)
 9283 format(31('*'),' SOLUTION REPORT ',31('*'),/
     +       ,'TESTSET: ',a8,', PROBLEM:',i3,', SOLVER: ',a6,
     +          ', COMPUTE',a24,/
     +       'Parameter values: ',4d15.8,(:,19x,4d15.8))
      write(*,fmt=9282) A,B,TOLER,A1,B1,A2,B2
 9282 format(4x,'Endpoints: A =',1p,d15.8,' B =',d15.8,8x,
     +       'Tolerance:',d9.2/
     +       10x,'BCs: A1=',d15.8,' B1=',d15.8,/
     +       15x,     'A2=',d15.8,' B2=',d15.8,/
     +       79('*'))

        if (FLAG.eq.DEPREC) write(*,'(/,15x,a,/)')
     +    '  ***Results of "deprecated" calculation'
        call REPORT(QEVREP(CALSET),QEFREP(CALSET),QSDREP(CALSET),
     *     SOLVER,IPROB,KLO,KHI,TOLER,IFAIL,
     *     EV,
     *     NXMESH,XMESH,EF,PDEF,
     *     NLMESH,LMESH,RHO,
     *     ELAPSE,OKEXIT)
      end if
C     Close the two log-files:
      close(20)
      close(21)
      go to 1000

C Unable to call SOLVIT because needed mesh was absent:
  288 write(*,'(6x,a)')'Please set needed data or choose another option'
      go to 1000

C   CHOICE=9, display/change mesh & "true" values
  290 write(*,advance='NO',fmt=9299)
 9299 format (/17x,'***SUBMENU: DISPLAY/CHANGE MESH or "TRUE" DATA***',/
     +       ' z Return to previous menu',/
     +       ' 1 Display current mesh or "true" data from file',/
     +       ' 2 Change/re-read current mesh or "true" data from file',/
     +       'Choose option: ')

      if(GETIR(CHOIC9,1,2)) go to (291,292),CHOIC9
      go to 1001

C     Display current cached mesh/"true" data:
  291 write(*,*)
      call DBDISP(CALSET,IPROB,KLO,KHI)
      go to 290
C     Update cached mesh/"true" data from database:
  292 write(*,*)
      call DBUPDT(CALSET,IPROB,KLO,KHI,FORCE=YES)
      call SPAUSE
      go to 290

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C                   *** End of Main Menu loop ***
C Use 1000 to give a pause before redisplaying menu, else 1001
 1000 call SPAUSE
 1001 go to 100

      end program SLDRIVER

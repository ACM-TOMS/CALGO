C***+****|****+****|* COPYRIGHT J D PRYCE 1998 **|****+****|****+****|**
      module DBMOD
      use SLCONSTS
      use SLPSET
      use SAFEIO
      use SLUTIL
      implicit none
C Unit for reading data files:
      integer, parameter:: INFIL=18
C***+****|****+****|****+ GENERAL DESCRIPTION ***|****+****|****+****|**
C DBMOD contains the routines
C   DBINIT
C   EVCHK
C   EFCHK
C   SDCHK
C   UPEVTR
C   UPEFTR
C   UPSDTR
C   UPXMUN
C   UPLMUN
C   GTXMTR
C   GTXMUN
C   GTLMTR
C   GTLMUN
C   EVSCAN
C   EVDAT
C   EFSCAN
C   EFDAT
C   SDSCAN
C   SDDAT
C & debugging routine
C   DBSUMM
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
C DBMOD handles 3 data objects that are read from an external file of
C "true" data values and 'cached' within the program:
C Name of |
C object  |Contents
C --------+--------------------------------------------------------------
C 'EVTRUE'|A set of "true" e-vs of a given problem: abstractly, a vector
C         |of k-values and a corresponding vector of lambda(k), and of
C         |believed errors in the lambda(k), where the
C         |k's are a subset of the current range KLO,KHI.
C         |
C 'EFTRUE'|An x-mesh (x(i), i=1, ..., n) plus a set of "true" e-fns of a
C         |given problem: abstractly, a matrix with rows labelled by k
C         |values, and columns by i which indexes the x(i). The k,i
C         |position holds u(x(i)) and pu'(x(i)) for the k-th e-fn. Again
C         |the k's are a subset of the current range KLO,KHI. They need
C         |have no relation to the k's of any "true" e-vs for the same
C         |problem.
C         |
C 'SDTRUE'|A lambda-mesh (lambda(j), j=1, ..., m) plus a corresponding
C         |vector of "true" Spectral Density Function (SDF) values
C         |rho(lambda(j) of a given problem, and a vector of believed
C         |errors in these values.
C         |
C --------+--------------------------------------------------------------
C Also 2 data objects not connected to an external file:
C --------+--------------------------------------------------------------
C 'XMSHUN'|A "uniform" x-mesh
C         |
C 'LMSHUN'|A "uniform" lambda-mesh
C --------+--------------------------------------------------------------
C
C Once an object has got a value, this value is not discarded till it
C "has to be". This is for user convenience: meshes and "true" data act
C as 'settings' that can be re-used many times. Just when a value is
C discarded, is driven by what objects are used by different values of
C Calculation-Choices setting CALSET, shown by the Y entries in this
C Table below.
C For each object we can also determine if 'no value is defined'.
C This is shown by the logical expression in the "'Undefd' flag" column
C in the table.
C          |           USAGE             |      ID data     |'Undefd'
C CALSET = | 1  2 | 3  4  5  6 | 7  8  9 |                  |  flag
C ---------+------+------------+---------+------------------+----------+
C 'EVTRUE' |    Y |          Y |         | %IPROB %KLO %KHI | %IPROB=-1
C 'EFTRUE' |      |       Y  Y |         | %IPROB %KLO %KHI | %IPROB=-1
C 'SDTRUE' |      |            |    Y  Y | %IPROB           | %IPROB=-1
C ---------+------+------------+---------+------------------+----------+
C 'XMSHUN' |      |    Y       |         | none             | %NMESH=-1
C 'LMSHUN' |      |            | Y       | none             | %NMESH=-1
C ---------+------+------------+---------+------------------+----------+
C where %IPROB in the ID data of 'EVTRUE' means a 'EVTRUE'%IPROB field
C in the (conceptual) record 'EVTRUE'.
C When the user invokes 'Solve' from the main menu the program checks
C whether the object(s) required by the current CALSET have values that
C are usable.
C E.g if CALSET=6, the table shows 'EVTRUE' and 'EFTRUE' will be used.
C For each, the "true" data cached within the program belongs to a
C specific problem & k-range, so the program checks that the %IPROB
C field *equals* the currently set IPROB and the %KLO:%KHI range
C *contains* the currently set KLO:KHI (so that any "true" e-v, e-fn
C data held on file for this k-range must be in the cache). If not, the
C database is accessed.
C
C It is not necessary for any "true" values for this IPROB and k-range
C to be present, but a valid x-mesh must be read, or the operation is
C aborted. Other CALSET values are similar.
C
C The items with no ID data can be used *unconditionally* provided a
C value is present (in case of 'XMSHUN' the program checks if the mesh
C  lies wholly within the current interval (A,B) and warns if not, but it
C can use the mesh anyway).
C
C When the user invokes 'Display/change' from the main menu, the
C program checks whether the object(s) required by the current CALSET
C (1) have defined values, (2) are usable, and reports summary data. The
C user then has the options
C - update the object(s), as in case of 'Solve'. If this fails for any
C   object, the current value is left unchanged. (For CALSET=6, the
C   result may be that 'EVTRUE' is updated and 'EFTRUE' not, or v.v.)
C - show values in more detail, if defined. E.g. for CALSET=5, the
C   x-mesh would be reported. For CALSET=6, each of "true" e-vs, x-mesh,
C   and "true" e-fn values can be requested.
C
C   x-mesh design points
C   --------------------
C   We expect users to wish to compare different e-fns on the *same*
C   x-mesh provided they come from the same Problem (same IPROB value)
C   in the following circumstances:
C   a. Different solver used
C   b. Changed e-v index
C   c. Changed BCs but same endpoints
C   d. Changed values of Problem-parameters
C   e. Changing the endpoints a,b from their "original" values.
C
C   - For case d recall that changing parameters can change an
C     endpoint's classification, e.g. from regular to singular.
C   - For case e, recall an endpoint that is singular in the Problem
C     as originally posed can only be moved *inward* from its original
C     value, but an originally regular point can be moved in either
C     direction.
C
C   Also, some solvers crash if an e-fn is evaluated *at* a singular
C   endpoint. But SLEDGE is explicitly designed not to do so, and in
C   fact the endpoints *must* be given it as part of any mesh, or it
C   will give an error exit.
C
C   To cope with cases d,e above, SOLVIT *copies* the x-mesh passed to
C   it and alters it before passing it to the solver as follows.
C   Remove all points <=a and >=b (where a,b are the *current*
C   endpoints). Then add a (and similarly b) in the mesh if it is
C   regular (in particular if it has been moved inward from the
C   original a value), or if the solver is SLEDGE. The resulting mesh
C   is output from SOLVIT for use by REPORT.
C
C   Inserting a,b like this improves the look of e-fn graphs!
C
C   AUTO meshes are treated as if the mesh is formed by the SOLVIT call,
C   used in the REPORT call and forgotten. It can't be re-used; nor
C   displayed, except as part of the e-fn printout in REPORT.
C
C***+****|****+****| Details of module variables |****+****|****+****|**
C DBPATH: Pathname relative to directory of executable program, where
C     data files are held.

C Concrete storage representing 'EVTRUE':
C   QEVPRB,QEVKLO,QEVKHI:
C     'ID-card' of "true" ev data currently cached in memory.
C   EVTRU: Row 1 holds "true" e-vs, row 2 holds believed absolute errors
C   Data for index k is in EVTRU(:,k-QEVKLO+1)
C   QEVTRU(K) says whether K-th "true" e-v is ABSENT or PRESNT
C Ancillary data for 'EVTRUE', used for reporting
C   NEVBLK: No. of data blocks in file EVTRU.nn
C   IEVBLK: index of block chosen by user, in range 1:NEVBLK

C Concrete storage representing 'EFTRUE':
C   QEFPRB,QEFKLO,QEFKHI:
C   'ID-card' of "true" efn data currently cached in memory.
C   NXMTR: no. of points in XMSHTR
C   XMSHTR,EFTR,PDEFTR:
C   x-mesh, "true" values of u(x), pu'(x) on x-mesh for each k
C   QEFTRU(K) says whether K-th "true" e-fn is ABSENT or PRESNT
C Ancillary data for 'EFTRUE', used for reporting
C   NEFBLK: No. of data blocks in file EFTRU.nn
C   IEFBLK: index of block chosen by user, in range 1:NEFBLK

C Concrete storage representing 'SDTRUE':
C   QSDPRB: 'ID card' of lambda-mesh currently cached in memory
C   NLMTR: no. of points in lambda-mesh
C   LMSHTR,SDTRU
C   lambda-mesh, "true" values of SDF rho(lambda) on mesh
C Ancillary data for 'SDTRUE', used for reporting
C   NSDBLK: No. of data blocks in file SDTRU.nn
C   ISDBLK: index of block chosen by user, in range 1:NSDBLK

C Concrete storage representing 'XMSHUN':
C   NXMUN: no. of points in UNIForm x-mesh
C   XMSHUN UNIForm x-mesh

C Concrete storage representing 'LMSHUN':
C   NLMUN: no. of points in UNIForm lambda-mesh
C   LMSHUN UNIForm lambda-mesh
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**

      character*40 DBPATH
C 'EVTRUE':
      integer QEVPRB,QEVKLO,QEVKHI
      double precision EVTRU(1:2,1:MAXEVS)
      integer QEVTRU(1:MAXEVS)
      integer NEVBLK,IEVBLK
C 'EFTRUE':
      integer:: QEFPRB,QEFKLO,QEFKHI
      integer NXMTR
      double precision XMSHTR(1:MAXMSH),EFTR(1:MAXMSH,1:MAXEVS)
     +  ,PDEFTR(1:MAXMSH,1:MAXEVS)
      integer QEFTRU(1:MAXEVS)
      integer NEFBLK,IEFBLK
C 'SDTRUE':
      integer QSDPRB
      integer NLMTR
      double precision LMSHTR(1:MAXMSH),SDTRU(1:2,1:MAXMSH)
      integer NSDBLK,ISDBLK
C 'XMSHUN':
      integer NXMUN
      double precision XMSHUN(1:MAXMSH)
C 'LMSHUN':
      integer NLMUN
      double precision LMSHUN(1:MAXMSH)
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      contains
C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function EVCHK(IPROB,KLO,KHI)
      integer IPROB,KLO,KHI
      EVCHK = QEVPRB.eq.IPROB .and. QEVKLO.eq.KLO .and. QEVKHI.eq.KHI
      end function EVCHK

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function EFCHK(IPROB,KLO,KHI)
      integer IPROB,KLO,KHI
      EFCHK = QEFPRB.eq.IPROB .and. QEFKLO.eq.KLO .and. QEFKHI.eq.KHI
      end function EFCHK

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function SDCHK(IPROB)
      integer IPROB
      SDCHK = QSDPRB.eq.IPROB
      end function SDCHK

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine DBINIT
C Set initial state of 'EVTRUE', 'EFTRUE', 'SDTRUE', 'XMSHUN', 'LMSHUN'
C to have no data.
      QEVPRB =-1
      QEFPRB =-1
      QSDPRB =-1
      NXMUN = -1
      NLMUN = -1
C and default k-range: does this make sense?
      QEVKLO = 0
      QEVKHI = 9
      QEFKLO = 0
      QEFKHI = 9
      end subroutine DBINIT

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine UPEVTR(IPROB,KLO,KHI,FORCE)
C UPEVTR updates DBMOD's 'EVTRUE'=(EVTRU,QEVTRU,QEVPRB,QEVKLO,QEVKHI).
C All arguments are input.
C If .not.FORCE and already present with 'ID'=(IPROB,KLO,KHI), nothing
C   is done.
C If .not.FORCE and present with same IPROB & enclosing [KLO,KHI] it
C   doesn't read database file but shifts existing data appropriately.
C Else, it gets from database file a block of "true" data with
C   this ID.
C   If no file exists, or user aborts, it returns without changing
C   anything.
C It sets ancillary 'EVTRUE' info if it reads file or tries to:
C   NEVBLK=no. of blocks on data file, 0 if no file
C   IEVBLK=index. of user-selected block, 0 if no file or user-abort

C     .. Scalar arguments ..
      integer IPROB,KLO,KHI
      logical FORCE
C ..  Local scalars ..
      integer I,OFFSET

      if (QEVPRB.eq.IPROB .and. .not.FORCE) then
        if (QEVKLO.eq.KLO .and. QEVKHI.eq.KHI) return
        if (QEVKLO.le.KLO .and. QEVKHI.ge.KHI) then
          OFFSET = KLO-QEVKLO
          do 100 I=1,KHI-KLO+1
            EVTRU(1,I) = EVTRU(1,I+OFFSET)
            EVTRU(2,I) = EVTRU(2,I+OFFSET)
            QEVTRU(I) = QEVTRU(I+OFFSET)
  100     continue
          QEVKLO = KLO
          QEVKHI = KHI
          return
        end if
      end if
C Else access database:
      call EVSCAN(IPROB,NEVBLK)
      IEVBLK = 0
      if (NEVBLK.eq.0) return
      write(*,'(''Which block do you want (1 to'',i2,
     +  ''), z to abort? '')',advance='NO') NEVBLK
      if (.not.GETIR(IEVBLK,1,NEVBLK)) return

C Else OK with NEVBLK, IEVBLK both >0 so:
C     Reset 'EVTRUE's ID data
      QEVPRB = IPROB
      QEVKLO = KLO
      QEVKHI = KHI
C     Read the actual data:
c     print*,'UPEVTR calling EVDAT'
      call EVDAT(IPROB,KLO,KHI,EVTRU,QEVTRU)

      end subroutine UPEVTR

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine UPEFTR(IPROB,KLO,KHI,FORCE)
C UPEFTR updates DBMOD's
C   'EFTRUE'=(NXMTR,XMSHTR,EFTR,PDEFTR,QEFTRU,QEFPRB,QEFKLO,QEFKHI).
C All arguments are input.
C If .not.FORCE and already present with 'ID'=(IPROB,KLO,KHI), nothing
C   is done.
C If .not.FORCE and present with same IPROB & *wider* k-range it
C   doesn't read database file but shifts existing data appropriately
C   to match the given k-range.
C Else, it gets from database file a block of "true" data with
C   this ID.
C   If no file exists, or user aborts, it returns without changing
C   'EFTRUE'.
C It sets ancillary 'EFTRUE' info if it reads file or tries to:
C   NEFBLK=no. of blocks on data file, 0 if no file
C   IEFBLK=index. of user-selected block, 0 if no file or user-abort

C     .. Scalar arguments ..
      integer IPROB,KLO,KHI
      logical FORCE
C     .. Local scalars ..
      integer I,J,OFFSET

      if (QEFPRB.eq.IPROB .and. .not.FORCE) then
        if (QEFKLO.eq.KLO .and. QEFKHI.eq.KHI) return
        if (QEFKLO.le.KLO .and. QEFKHI.ge.KHI) then
          OFFSET = KLO-QEFKLO
          do 30 I=1,KHI-KLO+1
            do 20 J=1,NXMTR
              EFTR(J,I) = EFTR(J,I+OFFSET)
              PDEFTR(J,I) = PDEFTR(J,I+OFFSET)
              QEFTRU(I) = QEFTRU(I+OFFSET)
   20       continue
   30     continue
          QEFKLO = KLO
          QEFKHI = KHI
          return
        end if
      end if

C Else access database:
      call EFSCAN(IPROB,NEFBLK)
      IEFBLK = 0
      if (NEFBLK.eq.0) return
      write(*,'(''Which block do you want (1 to'',i2,
     +  ''), z to abort? '')',advance='NO') NEFBLK
      if (.not.GETIR(IEFBLK,1,NEFBLK)) return

C Else OK with NEFBLK, IEFBLK both >0 so:
C     Reset 'EFTRUE's ID data
      QEFPRB = IPROB
      QEFKLO = KLO
      QEFKHI = KHI
C     Read the actual data:
c     print*,'UPEFTR calling EFDAT'
      call EFDAT(IPROB,KLO,KHI,NXMTR,XMSHTR,EFTR,PDEFTR,QEFTRU)

      end subroutine UPEFTR

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine UPSDTR(IPROB,FORCE)
C All arguments are input.
C UPSDTR updates DBMOD's 'SDTRUE'=(NLMTR,LMSHTR,SDTRU,QSDPRB).
C If .not.FORCE and already present with 'ID'=(IPROB), nothing is done.
C Else, it gets from database file a block of "true" data with
C   this ID.
C   If no file exists, or user aborts, it pretends it does and puts
C   empty info in 'SDTRUE'.
C It sets ancillary 'SDTRUE' info if it reads file or tries to:
C   NSDBLK=no. of blocks on data file, 0 if no file
C   ISDBLK=index. of user-selected block, 0 if no file or user-abort

C     .. Scalar arguments ..
      integer IPROB
      logical FORCE

      if (QSDPRB.eq.IPROB .and. .not.FORCE) return
C Else access database:
      call SDSCAN(IPROB,NSDBLK)
      ISDBLK = 0
      if (NSDBLK.eq.0) return
      write(*,'(''Which block do you want (1 to'',i2,
     +  ''), z to abort? '')',advance='NO') NSDBLK
      if (.not.GETIR(ISDBLK,1,NSDBLK)) return

C Else OK with NSDBLK, ISDBLK both >0 so:
C     Reset 'SDTRUE's ID data
      QSDPRB = IPROB
C     Read the actual data:
c     print*,'UPSDTR calling SDDAT'
      call SDDAT(IPROB,NLMTR,LMSHTR,SDTRU)
      end subroutine UPSDTR

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine UPXMUN(FORCE)
C Argument FORCE is input.
C UPXMUN updates DBMOD's 'XMSHUN'=(NXMUN,XMSHUN).
C It gets NXMSHI from user, default 10 if abort, forms uniform mesh of
C NXMUN=NXMSHI+2 points in XMSHUN.
C     .. Scalar arguments ..
      logical FORCE
C     .. Local scalars ..
      integer NXMSHI
      double precision A,B
      if (NXMUN.le.0 .or. FORCE) then
C       extract current endpoints of SLP from module SLPSET
        call GTABCU(A,B)
        write(*,9991,advance='NO') MAXMSH-2
 9991   format('How many interior points in x-mesh (not counting A,B)',/
     +   ' in range 10 to',i3,': ')
        if (.not. GETIR(NXMSHI,10,MAXMSH-2)) then
          NXMSHI = 10
          write(*,'(3x,a)') 'No. of points set to default of 10'
        end if
        NXMUN = NXMSHI+2
        write(*,'(3x,a,i3,a)')
     +    'Forming "UNIForm" mesh of ',NXMUN,' points'
C       Revised x-meshes spec means we always include A,B in UNIF mesh:
        call MKMESH(NXMUN,A,B,.TRUE.,.TRUE.,XMSHUN)
      end if
      end subroutine UPXMUN

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine UPLMUN(FORCE)
C Argument FORCE is input.
C UPLMUN updates DBMOD's 'LMSHUN'=(NLMUN,LMSHUN).
C It gets NLMUN from user, default 10 if abort, forms uniform mesh of
C NLMUN points, over [LLO,LHI], in LMSHUN. If LLO or LHI is infinite,
C it is excluded from the mesh but the no. of points stays the same.
C     .. Scalar arguments ..
      logical FORCE
C     .. Local scalars ..
      double precision LLO,LHI
      if (NLMUN.le.0 .or. FORCE) then
        write(*,9999,advance='NO') -XINFTY
 9999   format('Give lower endpoint of lambda-mesh',/
     +   ' (use ',1p,d10.2,'for minus infinity: ')
        if (.not. GETR(LLO)) then
          LLO=0d0
          write(*,*) '  Lower endpoint set to default of ',LLO
        end if
        write(*,9998,advance='NO') XINFTY
 9998   format('Give upper endpoint of lambda-mesh',/
     +   ' (use ',1p,d10.2,'for plus infinity: ')
        if (.not. GETRR(LHI,LLO,XINFTY)) then
          LHI=XINFTY
          write(*,*) '  Upper endpoint set to default of ',LHI
        end if
        write(*,9991,advance='NO') MAXMSH
 9991   format('How many points in lambda-mesh',/
     +   ' in range 1 to',i3,': ')
        if (.not. GETIR(NLMUN,1,MAXMSH)) then
          NLMUN = 10
          write(*,*) 'No. of points set to default of ',NLMUN
        end if
        write(*,'(3x,a,i3,a)')
     +    'Forming "UNIForm" mesh of ',NLMUN,' points'
        call MKMESH(NLMUN,LLO,LHI,LLO.gt.-XINFTY,LHI.lt.XINFTY,LMSHUN)
      end if
      end subroutine UPLMUN

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine DIEVTR(IPROB,KLO,KHI)
C DIEVTR displays DBMOD's 'EVTRUE'=(EVTRU,QEVTRU,QEVPRB,QEVKLO,QEVKHI).
C Also ancillary 'EVTRUE' info
C   NEVBLK=no. of blocks on data file, 0 if no file
C   IEVBLK=index of user-selected block, 0 if no block selected

C     .. Scalar arguments ..
      integer IPROB,KLO,KHI
C ..  Local scalars ..
      integer IK,K
      if (QEVPRB.eq.-1) then
        write(*,*)'No "true" eigenvalue data is present'
        return
      end if
      write(*,*)'Stored "true" eigenvalue data is block ',IEVBLK,
     +  ' of ',NEVBLK,' on file'
      if (QEVPRB.eq.IPROB) then
        if (QEVKLO.eq.KLO .and. QEVKHI.eq.KHI) then
          write(*,*)
     +      'Data matches current problem and k range'
        elseif (QEVKLO.le.KLO .and. QEVKHI.ge.KHI) then
          write(*,*)
     +      'Data matches current problem and ',
     +      'encloses current k range'
        else
          write(*,*)
     +      'Data matches current problem but ',
     +      'for different k range'
        end if
      else
        write(*,*) 'Data is not for current problem'
      end if
      write(*,*) 'Problem ',QEVPRB,' from k= ',KLO,' to k= ',KHI
      write(*,fmt="(3x,'K',5x,'Eigenvalue',13x,'Believed error')")
      do 100 K=KLO,KHI
        IK = K-KLO+1
        if (QEVTRU(IK).eq.PRESNT) write(*,
     +    fmt='(i4,3x,1p,g22.15,3x,d10.2)')K,EVTRU(1,IK),EVTRU(2,IK)
  100   continue
      end subroutine DIEVTR

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine DIEFTR(IPROB,KLO,KHI)
C DIEFTR displays DBMOD's 'EFTRUE'=
C                             (EFTR,PDEFTR,QEFTRU,QEFPRB,QEFKLO,QEFKHI).
C Also ancillary 'EFTRUE' info
C   NEFBLK=no. of blocks on data file, 0 if no file
C   IEFBLK=index of user-selected block, 0 if no block selected

C     .. Scalar arguments ..
      integer IPROB,KLO,KHI
C ..  Local scalars ..
      integer IK,J,K
      if (QEFPRB.eq.-1) then
        write(*,*)'No x-mesh / "true" eigenfunction data is present'
        return
      end if
      write(*,*)'Stored "true" eigenfunction data is block ',IEFBLK,
     +  ' of ',NEFBLK,' on file'
      if (QEFPRB.eq.IPROB) then
        if (QEFKLO.eq.KLO .and. QEFKHI.eq.KHI) then
          write(*,*)
     +      'Data matches current problem and k range'
        elseif (QEFKLO.le.KLO .and. QEFKHI.ge.KHI) then
          write(*,*)
     +      'Data matches current problem and ',
     +      'encloses current k range'
        else
          write(*,*)
     +      'Data matches current problem but ',
     +      'for different k range'
        end if
      else
        write(*,*) 'Data is not for current problem'
      end if
      write(*,*) 'Problem ',QEFPRB,' from k= ',KLO,' to k= ',KHI
      do 100 K=KLO,KHI
        IK = K-KLO+1
        if (QEFTRU(IK).eq.PRESNT) then
          write(*,fmt=9999) K
 9999     format (1x,'Stored Eigenfunction Values for k=',i4,/
     +            1x,8x,'x         u(x)              pu''(x)')
            do 90 J=1,NXMTR
              write(*,fmt=9997)XMSHTR(J),EFTR(J,IK),PDEFTR(J,IK)
 9997         format (1x,0p,f12.7,1p,2e18.10)
  90        continue
        end if
  100 continue
      end subroutine DIEFTR

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine DISDTR(IPROB)
C DISDTR displays DBMOD's 'SDTRUE'=(NLMTR,LMSHTR,SDTRU,QSDPRB).
C Also ancillary 'SDTRUE' info
C   NSDBLK=no. of blocks on data file, 0 if no file
C   ISDBLK=index. of user-selected block, 0 if no file or user-abort

C     .. Scalar arguments ..
      integer IPROB
C     .. Local scalars ..
      integer I

      if (QSDPRB.eq.-1) then
        write(*,*)
     +    'No lambda-mesh / "true" spectral density data is present'
        return
      end if
      write(*,*)'Stored "true" spectral density data is block ',
     +  ISDBLK,' of ',NSDBLK,' on file'
      if (QSDPRB.eq.IPROB) then
        write(*,*)'Data matches current problem'
      else
        write(*,*) 'Data is not for current problem'
      end if
      write(*,*) 'Problem ',QSDPRB
      write(*,fmt="(6x,'lambda',15x,'rho(lambda)',9x,'Believed error')")
      do 100 I=1,NLMTR
        write(*,fmt='(i2,3x,1p,2g22.15,3x,d10.2)')
     +    I,LMSHTR(I),SDTRU(1,I),SDTRU(2,I)
  100   continue
      end subroutine DISDTR

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine DIXMUN()
C DIXMUN displays DBMOD's 'XMSHUN'=(NXMUN,XMSHUN).
C     .. Local scalars ..
      integer I
      if (NXMUN.eq.-1) then
        write(*,*)'No "uniform" x-mesh data is present'
        return
      end if
      write(*,fmt=9999) NXMUN,XMSHUN(1),XMSHUN(NXMUN)
 9999 format('Stored ''UNIForm'' x-mesh has',i3,' points from x =',
     +  g10.5,' to x =',g10.5,/
     +  '   i          x(i)')
      do 100 I=1,NXMUN
        write(*,fmt='(i3,3x,1p,g22.15)') I,XMSHUN(I)
  100   continue
      end subroutine DIXMUN

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine DILMUN()
C DILMUN displays DBMOD's 'LMSHUN'=(NLMUN,LMSHUN).
C     .. Local scalars ..
      integer I
      if (NLMUN.eq.-1) then
        write(*,*)'No "uniform" lambda-mesh data is present'
        return
      end if
      write(*,fmt=9999) NLMUN,LMSHUN(1),LMSHUN(NLMUN)
 9999 format('Stored ''UNIForm'' lambda-mesh has',i3,
     +  ' points from lambda =',g10.5,' to lambda =',g10.5,/
     +  '   i          lambda(i)')
      do 100 I=1,NLMUN
        write(*,fmt='(i3,3x,1p,g22.15)') I,LMSHUN(I)
  100   continue
      end subroutine DILMUN

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function GTXMTR(NXMESH,XMESH)
C GTXMTR tries to copy 'EFTRUE's x-mesh to NXMESH,XMESH. If this mesh
C exists (QEFPRB>0), it copies it and returns TRUE, else it leaves
C NXMESH,XMESH unchanged and returns FALSE.

C     .. Scalar arguments ..
      integer NXMESH
C     .. Array arguments ..
      double precision XMESH(1:MAXMSH)
C     .. Local scalars ..
      integer I

      if (QEFPRB.gt.0) then
        GTXMTR = .TRUE.
        NXMESH = NXMTR
        do 102 I=1,NXMESH
          XMESH(I) = XMSHTR(I)
  102   continue
        return
      else
        GTXMTR = .FALSE.
      end if
      end function GTXMTR

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function GTXMUN(NXMESH,XMESH)
C GTXMUN tries to copy 'XMSHUN's x-mesh to NXMESH,XMESH. If this mesh
C exists (NXMUN>0), it copies it and returns TRUE, else it leaves
C NXMESH,XMESH unchanged and returns FALSE.

C     .. Scalar arguments ..
      integer NXMESH
C     .. Array arguments ..
      double precision XMESH(1:MAXMSH)
C     .. Local scalars ..
      integer I

      if (NXMUN.gt.0) then
        GTXMUN = .TRUE.
        NXMESH = NXMUN
        do 102 I=1,NXMESH
          XMESH(I) = XMSHUN(I)
  102   continue
        return
      else
        GTXMUN = .FALSE.
      end if
      end function GTXMUN

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function GTLMTR(NLMESH,LMESH)
C GTLMTR tries to copy 'SDTRUE's lambda-mesh to NLMESH,LMESH. If this mesh
C exists (QSDPRB>0), it copies it and returns TRUE, else it leaves
C NLMESH,LMESH unchanged and returns FALSE.

C     .. Scalar arguments ..
      integer NLMESH
C     .. Array arguments ..
      double precision LMESH(1:MAXMSH)
C     .. Local scalars ..
      integer I

      if (QSDPRB.gt.0) then
        GTLMTR = .TRUE.
        NLMESH = NLMTR
        do 102 I=1,NLMESH
          LMESH(I) = LMSHTR(I)
  102   continue
        return
      else
        GTLMTR = .FALSE.
      end if
      end function GTLMTR

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      logical function GTLMUN(NLMESH,LMESH)
C GTLMUN tries to copy 'LMSHUN's lambda-mesh to NLMESH,LMESH. If this
C mesh exists (NLMUN>0), it copies it and returns TRUE, else it leaves
C NLMESH,LMESH unchanged and returns FALSE.

C     .. Scalar arguments ..
      integer NLMESH
C     .. Array arguments ..
      double precision LMESH(1:MAXMSH)
C     .. Local scalars ..
      integer I

      if (NLMUN.gt.0) then
        GTLMUN = .TRUE.
        NLMESH = NLMUN
        do 102 I=1,NLMESH
          LMESH(I) = LMSHUN(I)
  102   continue
        return
      else
        GTLMUN = .FALSE.
      end if
      end function GTLMUN

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine EVSCAN(IPROB,NBLOCK)
      integer IPROB,NBLOCK
C EVSCAN scans the file of "true" eigenvalue data for Problem no.
C IPROB in the current test set, & displays on-screen the description
C line for each data-block in the file, also any Comment lines.
C It returns NBLOCK=no. of blocks in file, or 0 if file doesn't
C exist or has no data.
C
C The database for Problem no. nn (01<=nn<=60) is the textfile
C                          EVTRU.nn
C in directory "testset\truevals" where "testset" denotes the name of
C current test set. It need not exist: if it does not, it is treated
C as if it is empty.
C
C The file comprises a sequence of blocks, and may also have Comment
C lines with C in column 1 (these must not occur in the middle of a
C sequence of free-format data).
C Each block starts with a record having letter B in column 1, the rest
C of the line being treated as comment & displayed on screen.
C E.g.         B Evs with alpha=1, BCs u'(0)=0, u(pi)=0
C The block continues with a sequence of records with V in col. 1
C the rest of the line containing
C        Index k   Eigenvalue eig_k   Believed bound on absolute error
C in free-format, with the k values >=0 and sorted in increasing order.
C E.g.         V  1 4.1234567890 2.0d-10
C indicates that we believe ev no. 1 equals 4.1234567890 with error at
C most 2 units in the last place.
C
C Creating and adding to these files is the user's responsibility.
C The log-file output from the program, using suitably high
C accuracy, can often be trusted to provide the needed data.
C
C The file is identified by the Problem No. The individual data
C blocks carry no identification of the parameter values, range A,B and
C boundary conditions used, so if you use various values for these,
C identify them by comment lines! The program can't understand the
C comments, but displays them when it reads the file.
C
C Verifying the file is also the user's responsibility. In particular,
C make sure the eigenvalue records for each block have
C k's in ascending order. If records for a block are for k = 4 5 6 0 1 2
C in that order, and KLO..KHI is 2..5, then the routine will see the data
C for k= 4 and 5, and miss that for k = 2.
C
C Algorithm summary:
C Set block-count NBLOCK & record-count ILINE to 0
C Open input file, if it doesn't exist count as end-of-file
C loop
C   Read record-type (=1st character of record)
C   exit loop on end-of-file
C   if 'C' (comment), read & display rest of line
C   elseif not 'X', skip record
C   else increase NBLOCK & read description line & display it
C   end if
C end loop
C close input file & return

      integer INFIL
      parameter(INFIL=18)

      integer IP, IREC, IOERR
      character FNAME*43, REC_ID*1, LINE*72

      NBLOCK = 0
C    Set FNAME to 'evtru.nn' where nn = decimal representation of IPROB
C    This is a f77-ish trick for getting nn with leading zeros
C    Pathname is then added.
      write(FNAME,'(i3)') 100+IPROB
      call STTRIM(DBPATH,IP)
      FNAME = DBPATH(1:IP)//'evtru.'//FNAME(2:3)
      call STTRIM(FNAME,IP)
      open(INFIL,file=FNAME,status='OLD',action='READ',
     +  position='REWIND',iostat=IOERR)
      if (IOERR.ne.0) then
        write(*,'("! Database file ",a," not found",t79,"!")')
     +    FNAME(1:IP)
      else
        write(*,'("! Reading database file ",a,t79,"!")') FNAME
        IREC = 0
C       Main Loop over blocks, displaying B and V data for each:
150     continue
          read(INFIL,'(a1)',advance='NO',iostat=IOERR) REC_ID
          if (IOERR.lt.0 .or. REC_ID.eq.'E') goto 199
          IREC = IREC+1
c          print*, 'Record',IREC,' type "',REC_ID,'"'

          if (REC_ID.eq.'C') then
            read(INFIL,'(a)') LINE
            write(*,'("!   ",a,t79,"!")') LINE
          elseif (REC_ID.ne.'B') then
C           skip rest of line:
            read(INFIL,*)
          else
            read(INFIL,'(a)') LINE
            NBLOCK = NBLOCK+1
            write(*,'("! Block",i2,":",t79,"!",/
     +                "! ",a,t79,"!")') NBLOCK,LINE
          end if
C       end of loop
        goto 150
199     continue
      end if
      close(INFIL)
      end subroutine EVSCAN

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine EVDAT(IPROB,KLO,KHI,EVTRU,QEVTRU)
      integer IBLOCK,IPROB,KLO,KHI,QEVTRU(KLO:KHI)
      double precision EVTRU(1:2,KLO:KHI)
C EVDAT reads the file of "true" eigenvalue data for Problem no.
C IPROB in the current test set, up to the IEVBLK-th block where IEVBLK
C is a module variable.
C EVDAT returns in EVTRU a list comprising those ev's in the block
C whose index k lies in the range KLO..KHI.
C
C Specifically, if the database file for this problem doesn't exist
C nothing is altered. Otherwise, on exit from EVDAT:
C - For each ev in database whose k is in KLO..KHI:
C     EVTRU(1,K) holds eig_k;
C     EVTRU(2,K) holds a believed bound >=0 for the error in eig_k;
C     QEVTRU(K) holds the value PRESNT (=1)
C - For each other k in KLO..KHI
C     QEVTRU(K) holds the value ABSENT (=0)
C
C Algorithm summary:
C Open input file
C loop
C   Read to IEVBLK-th 'B'-record as in EVSCAN
C end loop
C loop
C   Read record-type ('B' or 'V' expected but could be garbage)
C   if end-of-file then count it as record-type 'E'
C                  else increase line-count
C   end if
C   if type is 'B' or 'E'
C       exit loop
C   elsif type 'V'
C     read V-data & insert in output arrays & display
C   else
C       ERROR message: unrecognized record-type
C   end if
C end loop
C Close input file

      integer IP, IREC, IOERR, K
      double precision EV,EVERR
C!      character FNAME*43, LINE*72, REC_ID*1
      character FNAME*43, REC_ID*1
C Make this a fatal error becuse it should never happen, and if it does,
C an array bound exception will soon occur in main program!
      if (KHI-KLO+1 .gt. MAXEVS) then
         write(*,'(/,a)') '*** EVDAT: Fatal error, KHI-KLO+1>MAXEVS'
         stop
      end if

C    Set FNAME to 'evtru.nn' where nn = decimal representation of IPROB
C    This is a f77-ish trick for getting nn with leading zeros
C    Pathname is then added.
      write(FNAME,'(i3)') 100+IPROB
      call STTRIM(DBPATH,IP)
      FNAME = DBPATH(1:IP)//'evtru.'//FNAME(2:3)
      call STTRIM(FNAME,IP)
      open(INFIL,file=FNAME,status='OLD',action='READ',
     +  position='REWIND',iostat=IOERR)

C     Initialize ev-data status info
      do 20 K=KLO,KHI
        QEVTRU(K) = ABSENT
 20   continue
C
C     Now count up to IEVBLK-th block:
      IBLOCK = 0
      IREC = 0
250   continue
        read(INFIL,'(a1)',advance='NO',iostat=IOERR) REC_ID
        IREC = IREC+1
        if (REC_ID.eq.'B') IBLOCK = IBLOCK+1
c      print*,'IREC IBLOCK REC_ID: ',IREC,IBLOCK,REC_ID
        if (IBLOCK.ge.IEVBLK) goto 299
        read(INFIL,*)
      goto 250
C     Now we are at 2nd character of desired block, skip rest of line:
299   continue
      read(INFIL,'(a)')
C     Read V-data in the block, skipping comment lines but
C     terminating at any other record-type:
      write(*,'(t1,"! Block",i2,": K values",t79,"!")') IEVBLK
350   continue
        read(INFIL,'(a1)',advance='NO',iostat=IOERR) REC_ID
        if (IOERR.ne.0) REC_ID = 'E'
        IREC = IREC+1
c        print*, 'Record',IREC,' type "',REC_ID,'"'
        if (REC_ID.eq.'C') then
          read(INFIL,*)
          goto 350
        end if
        if (REC_ID.ne.'V') goto 900
        read(INFIL,*,iostat=IOERR) K,EV,EVERR
        if (IOERR.ne.0) then
          write(*,*)
          write(*,*)'***EVDAT: error reading V-data at record',IREC
        elseif (K.ge.KLO .and. K.le.KHI) then
          EVTRU(1,K) = EV
          EVTRU(2,K) = EVERR
          QEVTRU(K) = PRESNT
          write(*,fmt='(i5)',advance='NO') K
        end if
      goto 350

900   write(*,*)
      close(INFIL)
      write(*,*)'... Done'
      end subroutine EVDAT

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine EFSCAN(IPROB,NBLOCK)
      integer IPROB,NBLOCK
C EFSCAN scans the file of "true" eigenfunction data for Problem no.
C IPROB in the current test set, & displays on-screen the Mesh
C belonging to each data-block in the file, also any Comment lines.
C It returns NBLOCK=no. of blocks in file, or 0 if file doesn't
C exist or has no data.
C
C The database for Problem no. nn (01<=nn<=60) is the textfile
C                          EFTRU.nn
C in directory "testset\truevals" where "testset" denotes the name of
C current test set. It need not exist: if it does not, it is treated
C as if it is empty.
C
C The file comprises a sequence of blocks, and may also have Comment
C lines with C in column 1 (these must not occur in the middle of a
C sequence of free-format data).
C Each block starts with a record having letter X in column 1
C and continuing with the free-format data (which can spill over onto
C further lines)
C       m   x(1)    ...    x(m)
C i.e. the value of NMESH and the points x(i) of the mesh.
C
C Convention:
C   x(1) may be given as -1.0d35 and x(m) as 1.0d35 (SLTSTPAK's
C   values for -infinity, +infinity), these will be set by the code
C   to the current value of a resp. b.
C
C Restrictions:
C   If you are setting up meshes to experiment with, you can do what you
C   like! But if setting up "true" eigenfunction values for general
C   use, observe these rules to avoid situations which some solvers
C   can't handle:
C 1.    a <= x(1) < x(2) < ... < x(m) <= b
C   where a, b are the endpoints of the problem (the current ones, which
C   may have been changed from the default ones!)
C 2.Equality at either end is forbidden UNLESS that end is regular or
C   weakly regular (of type 'R' or 'WR' in SLTSTPAK). Of course any
C   end that has been moved inward from the default end is of type 'R'.
C NOTE this excludes, for instance, x=-1, x=1 being included in "true"
C   data for the Legendre equation (Standard#20) even though all
C   efns (for the default BCs) are polynomials, and finite at these
C   points.

C The rest of the block is a sequence of records having letter U in
C column 1 and continuing with free-format data
C      k   u(1) pu'(1)    ...    u(m) pu'(m)
C i.e. the eigenvalue index k and the values of the k-th eigenfunction
C u_k(x) and derivative p(x)u_k'(x) at the points x=x(i) of the mesh.
C An example block is (with the X and U's assumed to be in col 1)
C    C This is a sample block
C    X 5 1.000 1.500 2.000 2.500 3.000
C    U 0 0.12345678 1.3567890
C        0.35673892 0.2121212
C        0.50012098 0
C        0.35673892 -0.2121212
C        0.12345678 -1.3567890
C
C Creating and adding to these files is the user's responsibility.
C The log-file output from the program, using suitably high
C accuracy, can sometimes be trusted to provide the needed data:
C but often it can't! Codes' eigenfunction values are still erratic.
C
C The file is identified by the Problem No. The individual data
C blocks carry no identification of the parameter values, range A,B and
C boundary conditions used, so if you use various values for these,
C identify them by comment lines! The program can't understand the
C comments, but displays them when it reads the file.
C
C Verifying the file is also the user's responsibility. In particular,
C make sure the eigenfunction records (U-records) for each block have
C k's in ascending order. If U-records for a block are for k = 4 5 6 0 1 2
C in that order, and KLO..KHI is 2..5, then the routine will see the data
C for k= 4 and 5, and miss that for k = 2.
C
C Algorithm summary:
C Set block-count NBLOCK & record-count ILINE to 0
C Open input file, if it doesn't exist count as end-of-file
C loop
C   Read record-type (=1st character of record)
C   exit loop on end-of-file
C   if 'C' (comment), read & display rest of line
C   elseif not 'X', skip record
C   else increase NBLOCK & read mesh & display it
C   end if
C end loop
C close input file & return

      integer IP, IREC, IOERR, J, NMESH
      double precision XMSHTR(1:MAXMSH)
      character FNAME*43, REC_ID*1, LINE*72

      NBLOCK = 0
C    Set FNAME to 'eftru.nn' where nn = decimal representation of IPROB
C    This is a f77-ish trick for getting nn with leading zeros
C    Pathname is then added.
      write(FNAME,'(i3)') 100+IPROB
      call STTRIM(DBPATH,IP)
      FNAME = DBPATH(1:IP)//'eftru.'//FNAME(2:3)
      call STTRIM(FNAME,IP)
      open(INFIL,file=FNAME,status='OLD',action='READ',
     +  position='REWIND',iostat=IOERR)
      if (IOERR.ne.0) then
        write(*,'("! Database file ",a," not found",t79,"!")')
     +    FNAME(1:IP)
      else
        write(*,'("! Reading database file ",a,t79,"!")') FNAME
        IREC = 0
C       Main Loop over blocks, displaying X and U data for each:
150     continue
          read(INFIL,'(a1)',advance='NO',iostat=IOERR) REC_ID
          if (IOERR.lt.0 .or. REC_ID.eq.'E') goto 199
          IREC = IREC+1
c          write(*,*)'Record',IREC,' type "',REC_ID,'"'

          if (REC_ID.eq.'C') then
            read(INFIL,'(a)') LINE
            write(*,'("!   ",a,t79,"!")') LINE
          elseif (REC_ID.ne.'X') then
C           skip rest of line:
            read(INFIL,*)
          else
            read(INFIL,*,iostat=IOERR) NMESH
     *          ,(XMSHTR(J),J=1,MIN(NMESH,MAXMSH))
            if (IOERR.ne.0) then
              write(*,*)'*** EFSCAN: error reading X-data at record'
     *          ,IREC
            else
              NBLOCK = NBLOCK+1
              write(*,'("! Block",i2,":",t79,"!")') NBLOCK
              if (NMESH.gt.MAXMSH) then
                write(*,*)'Mesh specified with',NMESH,' points,',
     *                    ' reduced to maximum allowed:',MAXMSH
                NMESH = MAXMSH
              end if
              write(*,'("!   has mesh of",i3," points:",t79,"!")') NMESH
              write(*,'("!",t79,"!",t4,1p6e12.4)')(XMSHTR(J),J=1,NMESH)
            end if
          end if
C       end of loop
        goto 150
199     continue
      end if
      close(INFIL)
      end subroutine EFSCAN

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine EFDAT(IPROB,KLO,KHI,NXMTR,XMSHTR,EFTR,PDEFTR,QEFTRU)
      integer IBLOCK,IPROB,KLO,KHI,NXMTR,QEFTRU(KLO:KHI)
      double precision EFTR(1:MAXMSH,KLO:KHI),
     *  PDEFTR(1:MAXMSH,KLO:KHI),XMSHTR(1:MAXMSH)
C EFDAT reads the file of "true" eigenfunction data for Problem no.
C IPROB in the current test set, up to the IEFBLK-th block where IEFBLK
C is a module variable, and returns in XMSHTR (starting at
C XMSHTR(1)) the x-mesh for that block, in EFTR a list of vectors
C comprising the values on the x-mesh of those eigenfunctions in the
C block for which data exists and whose index k lies in the range
C KLO..KHI, in PDEFTR the corresponding derivative values. QEFTRU(K)
C says whether data for index K is present.
C
C Specifically, on exit from EFDAT
C    NXMTR holds no. of points in x-mesh of user-chosen block
C    XMSHTR(1:NXMTR) holds their values
C while for each such k, as above:
C    EFTR(k,1:NXMTR) holds values of k-th eigenfunction on mesh
C    PDEFTR(k,1:NXMTR) holds values of k-th eigenfunction on mesh
C    QEFTRU(k) holds the value PRESNT (=1)
C For each other k in KLO..KHI
C    QEFTRU(k) holds the value ABSENT (=0)
C
C Algorithm summary:
C Open input file
C loop
C   Read to IEFBLK-th 'X'-record as in EFSCAN
C end loop
C Read X-mesh & store
C loop
C   Read record-type ('X' or 'U' expected but could be garbage)
C   if end-of-file then count it as record-type 'E'
C                  else increase line-count
C   end if
C   if type is 'X' or 'E'
C       exit loop
C   elsif type 'U'
C     read U-data & insert in output arrays & display
C   else
C       ERROR message: unrecognized record-type
C   end if
C end loop
C Close input file

      integer IP, IREC, IOERR, J, K
      double precision U(1:MAXMSH),PDU(1:MAXMSH)
      character FNAME*43, REC_ID*1

C    Set FNAME to 'eftru.nn' where nn = decimal representation of IPROB
C    This is a f77-ish trick for getting nn with leading zeros
C    Pathname is then added.
      write(FNAME,'(i3)') 100+IPROB
      call STTRIM(DBPATH,IP)
      FNAME = DBPATH(1:IP)//'eftru.'//FNAME(2:3)
      call STTRIM(FNAME,IP)
      open(INFIL,file=FNAME,status='OLD',action='READ',
     +  position='REWIND',iostat=IOERR)

C Make this a fatal error becuse it should never happen
      if (IOERR.ne.0 .or. IEFBLK.le.0) then
        write(*,*)
     +  '*** EFDAT: Fatal error, cannot open file or Block no. <= 0'
      stop
      end if
C Make this a fatal error becuse it should never happen, and if it does,
C an array bound exception will soon occur in main program!
      if (KHI-KLO+1 .gt. MAXEVS) then
        write(*,*)'*** EFDAT: Fatal error, KHI-KLO+1>MAXEVS'
        stop
      end if

C     Initialize e-fn data status info:
      do 20 K=KLO,KHI
        QEFTRU(K) = ABSENT
 20   continue
C
C     Now count up to IEFBLK-th block:
      IBLOCK = 0
      IREC = 0
250   continue
        read(INFIL,'(a1)',advance='NO',iostat=IOERR) REC_ID
        IREC = IREC+1
        if (REC_ID.eq.'X') IBLOCK = IBLOCK+1
c      print*,'IREC IBLOCK REC_ID: ',IREC,IBLOCK,REC_ID
        if (IBLOCK.ge.IEFBLK) goto 299
        read(INFIL,*)
      goto 250
C     Now we are at 2nd character of desired block:
299   continue
      read(INFIL,*,iostat=IOERR) NXMTR
     *      ,(XMSHTR(J),J=1,MIN(NXMTR,MAXMSH))
c      print*,'NXMTR: ',NXMTR
      if (IOERR.ne.0) then
        write(*,*)'EFDAT: error reading X-data at record',IREC
      else
        if (NXMTR.gt.MAXMSH) then
          NXMTR = MAXMSH
        end if
      end if
C     Read U-data in the block, skipping comment lines but
C     terminating at any other record-type:
      write(*,'(t1,"! Block",i2,": K values",t79,"!")') IEFBLK
350   continue
        read(INFIL,'(a1)',advance='NO',iostat=IOERR) REC_ID
        if (IOERR.ne.0) REC_ID = 'E'
        IREC = IREC+1
c        write(*,*)'Record',IREC,' type "',REC_ID,'"'
        if (REC_ID.eq.'C') then
          read(INFIL,*)
          goto 350
        end if
        if (REC_ID.ne.'U') goto 900
        read(INFIL,*,iostat=IOERR) K,(U(J),PDU(J),J=1,NXMTR)
        if (IOERR.ne.0) then
          write(*,*)'EFDAT: error reading U-data at record',IREC
        elseif (K.ge.KLO .and. K.le.KHI) then
          do 400 J=1,NXMTR
            EFTR(J,K) = U(J)
            PDEFTR(J,K) = PDU(J)
400       continue
          QEFTRU(K) = PRESNT
          write(*,fmt='(i5)',advance='NO') K
        end if
      goto 350

900   write(*,*)
c      print*,'exit EFDAT: NXMTR=',NXMTR
      close(INFIL)
      write(*,*)'... Done'
      end subroutine EFDAT

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine SDSCAN(IPROB,NBLOCK)
      integer IPROB,NBLOCK
C SDSCAN scans the file of "true" data for Spectral Density Function
C (SDF) rho(lambda) of Problem no. IPROB in the current test set, &
C displays on-screen the description line for each data-block in the
C file, also any Comment lines. It returns NBLOCK=no. of blocks in file,
C or 0 if file doesn't exist or has no data.
C
C The database for Problem no. nn (01<=nn<=60) is the textfile
C                          SDTRU.nn
C in directory "testset\truevals" where "testset" denotes the name of
C current test set. It need not exist: if it does not, it is treated
C as if it is empty.
C
C The file comprises a sequence of blocks, and may also have Comment
C lines with C in column 1 (these must not occur in the middle of a
C sequence of free-format data).
C Each block starts with a record having letter B in column 1, the rest
C of the line being treated as comment & displayed on screen.
C E.g.         B Problem with lambda-dependent BC coeffs = 1 0 0 1
C The block continues with a sequence of records with R in col. 1
C the rest of the line containing
C        lambda rho(lambda) Believed bound on absolute error
C in free-format, with the lambda values sorted in increasing order.
C E.g.         R  0.5 0.7324 2d-4
C indicates that we believe rho(0.5) equals 0.7324 with error at
C most 2 units in the last place.
C
C Creating and adding to these files is the user's responsibility.
C
C The file is identified by the Problem No. The individual data
C blocks carry no identification of the parameter values, range A,B and
C boundary conditions used, so if you use various values for these,
C identify them by comment lines! The program can't understand the
C comments, but displays them when it reads the file.
C
C Verifying the file is also the user's responsibility. In particular,
C make sure the SDF records for each block have lambda's in ascending
C order.
C
C Algorithm summary:
C Set block-count NBLOCK & record-count ILINE to 0
C Open input file, if it doesn't exist count as end-of-file
C loop
C   Read record-type (=1st character of record)
C   exit loop on end-of-file
C   if 'C' (comment), read & display rest of line
C   elseif not 'B', skip record
C   else increase NBLOCK & read description line & display it
C   end if
C end loop
C close input file & return

      integer IP, IREC, IOERR
      character FNAME*43, REC_ID*1, LINE*72

      NBLOCK = 0
C    Set FNAME to 'sdtru.nn' where nn = decimal representation of IPROB
C    This is a f77-ish trick for getting nn with leading zeros
C    Pathname is then added.
      write(FNAME,'(i3)') 100+IPROB
      call STTRIM(DBPATH,IP)
      FNAME = DBPATH(1:IP)//'sdtru.'//FNAME(2:3)
      call STTRIM(FNAME,IP)
      open(INFIL,file=FNAME,status='OLD',action='READ',
     +  position='REWIND',iostat=IOERR)
      if (IOERR.ne.0) then
        write(*,'("! Database file ",a," not found",t79,"!")')
     +    FNAME(1:IP)
      else
        write(*,'("! Reading database file ",a,t79,"!")') FNAME
        IREC = 0
C       Main Loop over blocks, displaying B and R data for each:
150     continue
          read(INFIL,'(a1)',advance='NO',iostat=IOERR) REC_ID
          if (IOERR.lt.0 .or. REC_ID.eq.'E') goto 199
          IREC = IREC+1
c          print*, 'Record',IREC,' type "',REC_ID,'"'

          if (REC_ID.eq.'C') then
            read(INFIL,'(a)') LINE
            write(*,'("!   ",a,t79,"!")') LINE
          elseif (REC_ID.ne.'B') then
C           skip rest of line:
            read(INFIL,*)
          else
            read(INFIL,'(a)') LINE
            NBLOCK = NBLOCK+1
            write(*,'("! Block",i2,":",t79,"!",/
     +                "! ",a,t79,"!")') NBLOCK,LINE
          end if
C       end of loop
        goto 150
199     continue
      end if
      close(INFIL)
      end subroutine SDSCAN

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine SDDAT(IPROB,NLMTR,LMSHTR,SDTRU)
      integer IPROB,NLMTR
      double precision LMSHTR(1:MAXMSH),SDTRU(1:2,1:MAXMSH)
C SDDAT reads the file of "true" SDF data for Problem no. IPROB in the
C current test set, up to the ISDBLK-th block where ISDBLK is a module
C variable.
C SDDAT returns in LMSHTR a mesh of lambda values in ascending order, in
C SDTRU(1,:) the corresponding SDF values rho(lambda) and in SDTRU(2,:),
C believed bounds on the absolute errors in the SDF values.
C
C If the database file for this problem doesn't exist nothing is
C altered.
C
C Algorithm summary:
C Open input file
C loop
C   Read to ISDBLK-th 'B'-record as in SDSCAN
C end loop
C loop
C   Read record-type ('B' or 'R' expected but could be garbage)
C   if end-of-file then count it as record-type 'E'
C                  else increase line-count
C   end if
C   if type is 'B' or 'E'
C       exit loop
C   elsif type 'R'
C     read R-data & insert in output arrays & display
C   else
C       ERROR message: unrecognized record-type
C end if
C end loop
C Close input file

      integer IBLOCK, IP, IREC, IOERR
      double precision LAMBDA, SD,SDERR
      character FNAME*43, LINE*72, REC_ID*1
      logical RECOK

C    Set FNAME to 'sdtru.nn' where nn = decimal representation of IPROB
C    This is a f77-ish trick for getting nn with leading zeros
C    Pathname is then added.
      write(FNAME,'(i3)') 100+IPROB
      call STTRIM(DBPATH,IP)
      FNAME = DBPATH(1:IP)//'sdtru.'//FNAME(2:3)
      call STTRIM(FNAME,IP)
      open(INFIL,file=FNAME,status='OLD',action='READ',
     +  position='REWIND',iostat=IOERR)

C
C     Now count up to ISDBLK-th block:
      IBLOCK = 0
      IREC = 0
250   continue
        read(INFIL,'(a1)',advance='NO',iostat=IOERR) REC_ID
        IREC = IREC+1
        if (REC_ID.eq.'B') IBLOCK = IBLOCK+1
c   print*,'IREC IBLOCK REC_ID: ',IREC,IBLOCK,REC_ID
        if (IBLOCK.ge.ISDBLK) goto 299
        read(INFIL,*)
      goto 250
C     Now we are at 2nd character of desired block:
299   continue
      read(INFIL,'(a)') LINE
      write(*,'(t1,"! Block",i2,": lambda values",t79,"!")') ISDBLK
C     Read R-data in the block, skipping comment lines but
C     terminating at any other record-type:
      NLMTR = 0
350   continue
        read(INFIL,'(a1)',advance='NO',iostat=IOERR) REC_ID
        if (IOERR.ne.0) REC_ID = 'E'
        IREC = IREC+1
c        print*, 'Record',IREC,' type "',REC_ID,'"'
        if (REC_ID.eq.'C') then
          read(INFIL,*)
          goto 350
        end if
        if (REC_ID.ne.'R') goto 900
        RECOK = .TRUE.
        read(INFIL,*,iostat=IOERR) LAMBDA,SD,SDERR
        if (IOERR.ne.0) then
          RECOK = .FALSE.
          write(*,*)'*** SDDAT: error reading R-data at record',IREC
        elseif (NLMTR .ge. MAXMSH) then
          RECOK = .FALSE.
          write(*,*)'SDDAT: more than ',MAXMSH,
     +      ' values in SDF Data Block, excess ignored'
        elseif (NLMTR.gt.0) then
          if (LAMBDA.le.LMSHTR(NLMTR)) then
            RECOK = .FALSE.
            write(*,*)'SDDAT: lambda value ',LAMBDA,
     +        ' in SDF Data Block not in asc. order will be ignored'
          end if
        end if
        if (RECOK) then
          NLMTR = NLMTR+1
          LMSHTR(NLMTR) = LAMBDA
          SDTRU(1,NLMTR) = SD
          SDTRU(2,NLMTR) = SDERR
          write(*,fmt='(1pg10.4)',advance='NO') LAMBDA
        end if
      goto 350

900   write(*,*)
      close(INFIL)
      write(*,*)'... Done'
      end subroutine SDDAT

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**
      subroutine DBSUMM
      use SAFEIO
C Crude routine to display ID information for cached data

      print*,'EVTRUE:'
      print*,'   QEVPRB,QEVKLO,QEVKHI=',QEVPRB,QEVKLO,QEVKHI
C      double precision EVTRU(1:2,1:MAXEVS)
      print*,'   QEVTRU(1:QEVKHI-QEVKLO+1)=',QEVTRU(1:QEVKHI-QEVKLO+1)
      print*,'   NEVBLK,IEVBLK=',NEVBLK,IEVBLK

      print*,'EFTRUE:'
      print*,'   QEFPRB,QEFKLO,QEFKHI=',QEFPRB,QEFKLO,QEFKHI
      print*,'   NXMTR=',NXMTR
      print*,'   XMSHTR(1:NXMTR)=',XMSHTR(1:NXMTR)
C     +  ,EFTR(1:MAXMSH,1:MAXEVS),PDEFTR(1:MAXMSH,1:MAXEVS)
      print*,'   QEFTRU(1:QEFKHI-QEFKLO+1)=',QEFTRU(1:QEFKHI-QEFKLO+1)
      print*,'   NEFBLK,IEFBLK=',NEFBLK,IEFBLK

      print*,'SDTRUE:'
      print*,'   QSDPRB=',QSDPRB
      print*,'   NLMTR=',NLMTR
      print*,'   LMSHTR(1:NLMTR)=',LMSHTR(1:NLMTR)
      print*,'   NSDBLK,ISDBLK=',NSDBLK,ISDBLK

      print*,'XMSHUN:'
      print*,'   NXMUN=',NXMUN
      print*,'   XMSHUN(1:NXMUN)=',XMSHUN(1:NXMUN)

      print*,'LMSHUN:'
      print*,'   NLMUN=',NLMUN
      print*,'   LMSHUN(1:NLMUN)=',LMSHUN(1:NLMUN)

      call SPAUSE
      end subroutine DBSUMM

      end module DBMOD

C***+****|****+****|* COPYRIGHT J D PRYCE 1998 **|****+****|****+****|**
C A simple driver to show the use of SLTSTPAK. It uses SL02F as solver.
C The problem-set is arbitrary and can be chosen at link-time.
      program EVSIMP
      use SLTSTPAK
      use MARCOMOD
      implicit none
      integer IWK,NTEMP
      parameter (IWK=10000,NTEMP=400)
      double precision A,A1,A2,B,B1,B2,CURREN,ELAPSE,TOL,XINFTY
      integer I,IFAIL,IPROB,K,N,NEPRM,NPARM,NPROB,KNOBS(3)
      logical SYM
      character ATYPE*4,BTYPE*4,PARNM*72,TITLE*72,TSETNM*8
      double precision ELAM(0:1),PARM(0:10),WK(0:IWK,1:4),
     +  WKSMAL(0:NTEMP,1:7)
      data KNOBS/3*0/
      common/SL02CM/ AINFO,BINFO
      character AINFO*1,BINFO*1
      external GETBC

C SOLVER-INDEPENDENT STARTUP PHASE
C   Call the general initialization routine:
      call TSTINI(TSETNM,NPROB,XINFTY)
C   Select the Problem
      write (*,FMT='('' Problem number (1 to'',i3,''): '')') NPROB
      read (*,FMT=*) IPROB
C   Call the first section of set-up code:
      call SETUP0(IPROB,TITLE,NPARM,NEPRM,PARNM)
C   Write the problem title and ask for parameters if any:
      write (*,FMT=*) TITLE
      if (NPARM.gt.0) then
         write (*,FMT=*) 'Give values of these parameters: ', PARNM
         read (*,FMT=*) (PARM(I),I=1,NPARM)
      end if
C   Call the second section of set-up code:
      call SETUP1(PARM,A,B,ATYPE,BTYPE,A1,A2,B1,B2,SYM)
C   Set eigenvalue index & tolerance:
      write (*,FMT='('' Eigenvalue index: '')')
      read (*,FMT=*) K
      write (*,FMT='('' Tolerance: '')')
      read (*,FMT=*) TOL

C  SOLVER-SPECIFIC CODE
C   Give initial guess of ELAM(0) & its error:
      ELAM(0) = 0d0
      ELAM(1) = 1d0
C   Convert information about endpoints into form needed by SL02F:
      call TRANSL(A.ne.-XINFTY,ATYPE,AINFO)
      call TRANSL(B.ne.XINFTY,BTYPE,BINFO)
C   Give max no. of meshpoints to be used by SL02F
      N = 3000
C   Select soft failure option:
      IFAIL = -1
C   Set the CPU clock running:
      call CPU(0,CURREN,ELAPSE)
C   Solve the problem:
      call SL02F(ELAM,A,B,K,AINFO,BINFO,SYM,TOL,COEFFN,GETBC,N,WK,IWK,
     +           WKSMAL,NTEMP,KNOBS,IFAIL)
C   Read the CPU clock:
      call CPU(1,CURREN,ELAPSE)
C   Report results:
      write (*,FMT=*) 'Final estimate of eigenvalue:',ELAM(0)
      write (*,FMT='('' Error estimate:'',1pd9.2)') ELAM(1)
      write (*,FMT=*) 'Number of function evaluations:',NEVAL(0)
      write (*,FMT=*) 'Number of meshpoints used during solution:',N
      write (*,FMT=*) 'Exit value of IFAIL:',IFAIL
      write (*,FMT='('' CPU time (sec): '',f10.5)') ELAPSE
      if (IFAIL.ne.0) write (*,FMT=*) ' WARNING: IFAIL IS NONZERO '
      end

      subroutine TRANSL(XFINIT,XTYPE,XINFO)
C Translates SLTSTPAK endpoint-type info into SL02F form.
      logical XFINIT
      character XINFO*1,XTYPE*4
      if (XTYPE.eq.'R   ' .or. XTYPE.eq.'RS  ') then
         XINFO = 'R'
      else if (XFINIT) then
         XINFO = 'S'
      else
         XINFO = 'I'
      end if
      end

      subroutine GETBC(Y,PDY,EIG,X,IEND,ISING)
C GETBC interfaces between SL02F and GETBCS, putting the arguments in a
C different order. It also extracts the ISING data via SL02CM.
      use SLTSTPAK,only:GETBCS
      implicit none
      double precision EIG,PDY,X,Y
      integer IEND
      logical ISING
      common/SL02CM/ AINFO,BINFO
      character AINFO*1,BINFO*1

      call GETBCS(IEND,X,EIG,PDY,Y)
      if (IEND.eq.0) then
        ISING = AINFO.ne.'R'
      else
        ISING = BINFO.ne.'R'
      end if
      end

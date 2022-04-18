C***+****|****+****|* COPYRIGHT J D PRYCE 1998 **|****+****|****+****|**
      module SLTSTVAR
      double precision,parameter::
     +         PI=3.14159265358979,
     +         XINFTY=HUGE(1D0)
c      'tstcom1.f'
      double precision A1,A2,A,B1,B2,B,P,PDU,Q,U,W,X
      integer IPROB,NFEVAL,NPROB
      logical SYM
c      'tstcom2.f'
      character ATYPE*4,BTYPE*4,PARNM*72,TITLE*72
c      'tstcom3.f'
      integer IPARM,NEPRM,NPARM
      double precision PARM(0:10),EIG
      equivalence(EIG,PARM(0))
      end module SLTSTVAR

C***+****|****+****|****+****|****+****|****+****|****+****|****+****|**

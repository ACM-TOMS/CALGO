!  updated 04/26/2004
!    This is a program which tests the installation 
!    The Airy function pacakge by B.R.Fabijonas.  
!      reference:  B.R. Fabijonas, Algorithm XXX: Airy functions,
!       ACM Trans. Math. Softw., Vol. xxx, pp. xxx-xxx.
!    Please report all bug reports to bfabi@mail.smu.edu
!
!  This programs tests all featuers of the package:  
!    Airy functions of real argument computed by integration of ODE
!    Airy functions or real argument computed by summation of 
!      asymptotic expansions
!    Airy functions of complex argument computed by integration of ODE
!    Airy functions or complex argument computed by summation of 
!      asymptotic expansions
!    Real and complex zeros of the Airy functions computed by table lookup
!    Real and complex zeros of the Airy functions computed by 
!      summation of asymptotic expansions
!    Vectorized version of zeros of Ai(x) and its derivative with 
!      associated values
!    Vectorized version of zeros of Bi(x) and its derivative with 
!      associated values
!    Vectorized version of zeros of Bi(z) and its derivative with
!      associated values
!    Auxiliary functions for negative real arguments computed by 
!      integration of ODE
!    Auxiliary functions for negative real arguments computed by 
!      summation of asymptotic expansions
!
! A module to interface the error computation and switch
!   between relative precision and absolute error.
!   the endings:
!    s = real single
!    d = real double
!    c = complex single
!    z = complex double
!
       module error
       implicit none
       private
       public errg
       integer, parameter :: prs = kind(0.0e0)
       integer, parameter :: prd = kind(0.0d0)
       real(prd), parameter :: zerod = 0.0_prd, &
            pid = 3.1415926535897932384626433832795028841971_prd
       real(prs), parameter :: zeros = 0.0_prs, pis = pid
       interface errg
         module procedure errds
         module procedure errs
         module procedure errd
         module procedure errzc
         module procedure errz
         module procedure errc
       end interface
       interface psi
	 module procedure psis
	 module procedure psid
       end interface
       contains
!
! real double and single
         real(prd) function errds(arg1,arg2)
         real(prd) arg1
         real(prs) arg2
         errds = abs(psi(abs(arg1)) - psi(abs(arg2)))
         end function errds
!
! real single
         real(prs) function errs(arg1,arg2)
         real(prs) arg1, arg2
         errs = abs(psi(arg1) - psi(arg2))
         end function errs
!
! real double
         real(prd) function errd(arg1,arg2)
         real(prd) arg1, arg2
         errd = abs(psi(arg1) - psi(arg2))
         end function errd
!
! complex double and single
         real(prd) function errzc(arg1,arg2)
         complex(prd) arg1, arg3
         complex(prs) arg2
         real(prd) x1,x2
         arg3 = arg2
         if ( abs(arg1) < 1.0_prd .or. abs(arg3) < 1.0_prd ) then
           errzc = abs(arg1-arg2)
         else
           x1 = atan2(aimag(arg1),real(arg1,prd))
           x2 = atan2(aimag(arg3),real(arg3,prd))
           x2 = mod(abs(x1-x2),2.0_prd*pid)
           errzc = abs(cmplx(log(abs(arg1/arg3)),min(x2,2.0_prd*pid-x2),prd))
         end if
         end function errzc
!
! complex single
         real(prs) function errc(arg1,arg2)
         complex(prs) arg1, arg2
         real(prs) x1,x2
         if ( abs(arg1) < 1.0_prs .or. abs(arg2) < 1.0_prs ) then
           errc = abs(arg1-arg2)
         else
           x1 = atan2(aimag(arg1),real(arg1))
           x2 = atan2(aimag(arg2),real(arg2))
           x2 = mod(abs(x1-x2),2.0_prs*pis)
           errc = abs(cmplx(log(abs(arg1/arg2)),min(x2,2.0_prs*pis-x2)))
         end if
         end function errc
!
! complex double
         real(prd) function errz(arg1,arg2)
         complex(prd) arg1, arg2
         real(prd) x1,x2
         if ( abs(arg1) < 1.0_prd .or. abs(arg2) < 1.0_prd ) then
           errz = abs(arg1-arg2)
         else
           x1 = atan2(aimag(arg1),real(arg1,prd))
           x2 = atan2(aimag(arg2),real(arg2,prd))
           x2 = mod(abs(x1-x2),2.0_prd*pid)
           errz = abs(cmplx(log(abs(arg1/arg2)),min(x2,2.0_prd*pid-x2),prd))
         end if
         end function errz
!
! real single
         real(prs) function psis(arg)
         real(prs) arg, at
	 at = abs(arg)
         if (at > 1.0_prs) then
           psis = sign(1.0_prs + log(at),arg)
	 else
	   psis = arg
	 end if
         end function psis
!
! real double
         real(prd) function psid(arg)
         real(prd) arg, at
         at = abs(arg)
         if (at > 1.0_prd) then
           psid = sign(1.0_prd + log(at),arg)
	 else
	   psid = arg
	 end if
         end function psid
       end module error
!
!
! program to test installation
!
      program test
      use airy_functions_real_single
      use airy_functions_real_double
      use airy_functions_complex_single
      use airy_functions_complex_double
      use error
      implicit none
      integer, parameter :: prs = kind(0.0e0)
      integer, parameter :: prd = kind(0.0d0)
      real(prd), parameter :: oned = 1.0_prd, thvd = 1.50_prd
      real(prs), parameter :: ones = 1.0_prs, thvs = 1.50_prs
      real(prs) aias, daias, aibs, daibs, rms, pis, g1, tols, ffacs, hs, rs
      real(prd) x, aix, daix, aiass, daiass, rmd, pid, f1, told, ffacd, hd, rd
      real (prd), dimension(5) ::  aizv, daizv, aiassv, daiassv
      complex(prs) zs, ciunits, aiasz, daiasz, aibsz, daibsz
      complex(prd) z, aiz, daiz, aiassz, daiassz, ciunitd
      complex(prd), dimension(5) ::   bizv, dbizv, biassv, dbiassv
      real(prd), dimension(20) :: ray
      complex(prd), dimension(20) :: airay, dairay
      character(3) ai, dai, bi, dbi
      integer ierr, itol, itol2
      integer ierrcnt, i, iptcnt, ilocerrcnt, ilocptcnt, j, nt, na, np, nam, nap
      open(19,file="airyerr.log",status='unknown')

      pid = 3.14159265358979323846264338327950288419716939937510582097494459230782_prd
      pis = real(pid,prs)
      ciunits = cmplx(0.0_prs,1.0_prs,prs)
      ciunitd = cmplx(0.0_prd,1.0_prd,prd)
!
! factor for testing computed values agains stored values or for computing
!   values at the same point in different precisions
      itol = 100
!
! factor for extremely close, but different, argument values
      itol2 = 100
!
! machine epsilons
      tols = exp((1-digits(tols))*log(real(radix(tols),prs)))
      told = exp((1-digits(told))*log(real(radix(told),prd)))
!
! set the counters and declare a few variables
      ierrcnt = 0
      iptcnt = 0
      ai  = ' Ai'
      dai = 'dAi'
      bi  = ' Bi'
      dbi = 'dBi'
      print*, ' ---------------------------------------------------'
      print*, ' '
      print*, 'This is a program which tests the installation '
      print*, '  the Airy function pacakge by B.R. Fabijonas.'
      print*, '  Please send all bug reports to bfabi@smu.edu'
      print*, ' '
!
! write some preliminaries to the output file
      write(19,*) ' This is a log file for the installation of the Airy function'
      write(19,*) ' package by B.R. Fabijonas.  Please send all bug reports'
      write(19,*) ' to bfabi@smu.edu.  The reference for this package is '
      write(19,*) ' B.R. Fabijonas 2004. Algorithm XXX:  Airy functions. '
      write(19,*) ' ACM Trans. Math. Softw., Vol.xxx, pp.xxx-xxx.'
      write(19,*) ' A detailed description of the computation method and parameter'
      write(19,*) ' determination can be found in B.R. Fabijonas, D.W. Lozier,'
      write(19,*) ' and F.W.J. Olver, Computation of complex Airy functions and '
      write(19,*) ' their zeros using asymptotics and the differential equation. '
      write(19,*) ' ACM Trans. Math. Softw., Vol.xxx, pp. xxx-xxx.'
      write(19,*) ' '
      write(19,'("***Parameter values used in the code which are compiler specific.")')
      write(19,'("  We use the abbreviations no. for number of terms and ")')
      write(19,'("   a.e. for asymptotic expansion.")')
      write(19,'("  Single precision computations: standard Airy functions")')
      call airy_info(rms, hs, nt, na, np)
      write(19,'(A45,e24.16)') '  machine epsilon', epsilon(rms)
      write(19,'(A45,e24.16)') '  radius of the delineating circle rho', rms
      write(19,'(A45,e24.16)') '  integration step size h', hs
      write(19,'(A60,i4)') '  no. in the Taylor series N',nt
      write(19,'(A60,i4)') '  no. in the a.e.s S', na
      write(19,'(A60,i4/)') '  number of partitions of the integration ray P', np
      write(19,'("  Single precision computations: auxiliary functions")')
      call airy_aux_info ( rs, hs, nt, nam, nap, np)
      write(19,'(A45,e24.16)') '  radius of the delineating circle rho_m ',rs
      write(19,'(A45,e24.16)') '  integrationstep size h_m ',hs
      write(19,'(A60,i4)') '  no. of terms in the Taylor series N_m ',nt
      write(19,'(A60,i4)') '  no. of terms in the a.e.s of the modulus functions S_m ',nam
      write(19,'(A60,i4)') '  no. of terms in the a.e.s of the phase functions S_theta ',nap
      write(19,'(A60,i4/)') '  no. of partitions of the integration ray P_m ',np
      write(19,'("  Double precision computations: standard Airy functions")')
      call airy_info(rmd, hd, nt, na, np)
      write(19,'(A45,e24.16)') '  machine epsilon ', epsilon(rmd)
      write(19,'(A45,e24.16)') '  the parameter rho ', rmd
      write(19,'(A45,e24.16)') '  step size h ', hd
      write(19,'(A60,i4)') '  no. in the Taylor series N ',nt
      write(19,'(A60,i4)') '  no. used in the a.e.s S ', na
      write(19,'(A60,i4/)') '  number of partitions of the integration ray P ', np
      write(19,'("  Double precision computations: auxiliary functions")')
      call airy_aux_info ( rd, hd, nt, nam, nap, np)
      write(19,'(A45,e24.16)') '  the parameter rho_m ',rd
      write(19,'(A45,e24.16)') '  step size h_m ',hd
      write(19,'(A60,i4)') '  no. in the t.s. N_m ',nt
      write(19,'(A60,i4)') '  no. in the a.e. of the modulus functions S_m ',nam
      write(19,'(A60,i4)') '  no. in the a.e. of the phase functions S_theta ',nap
      write(19,'(A60,i4//)') '  number of partitions of the integration ray P_m ',np
!
!*****************
! testing real routines
!****************
      ilocerrcnt = 0
      ilocptcnt = 0
      write(6,'(" Comparing the function values returned by the")')
      write(6,'("   real routines against stored values")',advance="NO")
      write(19,'( "***Comparing the function values returned by the")')
      write(19,'("  real routines against stored values")')
      x = 5.82340_prd  ! integration routine
      ffacd = oned !10.0_prd**max(oned,thvd*log(abs(x)))
      call airy_ai(x, aix, daix, ierr)
      call airy_ai(real(x,prs),aias,daias,ierr)
      write(19,1000) ai, x
      write(19,1001) aias, aix, 0.1539272803479132e-04_prd
      f1 = errg(aix,0.1539272803479132e-04_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      f1 = errg(daix,-0.3777943926223622e-04_prd)
      write(19,1000) dai, x
      write(19,1001) daias, daix, -0.3777943926223622e-04_prd
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      call airy_bi(x, aix, daix, ierr)
      call airy_bi(real(x,prs),aias,daias,ierr)
      write(19,1000) bi, x
      write(19,1001) aias, aix, 0.4288114816683048e+04_prd
      f1 = errg(aix,0.4288114816683048e+04_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1000) dbi, x
      write(19,1001) daias, daix, 0.1015462058214280e+05_prd
      f1 = errg(daix,0.1015462058214280e+05_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      ilocptcnt = ilocptcnt + 4
      x = -17.345890_prd  ! asymptotic expansion routine
      ffacd = oned !10.0_prd**max(oned,thvd*log(abs(x)))
      call airy_ai(x, aix, daix, ierr)
      call airy_ai(real(x,prs),aias,daias,ierr)
      write(19,1000) ai, x
      write(19,1001) aias, aix, -0.2677772589719993e+00_prd
      f1 = errg(aix,-0.2677772589719993e+00_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1000) dai, x
      write(19,1001) daias, daix, -0.2900295189410030e+00_prd
      f1 = errg(daix,-0.2900295189410030e+00_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      call airy_bi(x, aix, daix, ierr)
      call airy_bi(real(x,prs),aias,daias,ierr)
      write(19,1000) bi, x
      write(19,1001) aias, aix, 0.6870907209572491e-01_prd
      f1 = errg(aix,0.6870907209572491e-01_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1000) dbi, x
      write(19,1001) daias, daix, -0.1114292633371775e+01_prd
      f1 = errg(daix,-0.1114292633371775e+01_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      ilocptcnt = ilocptcnt + 4
      write(19,1500) ilocptcnt
      write(19,1501) ilocerrcnt
      iptcnt = iptcnt + ilocptcnt
      ierrcnt = ierrcnt + ilocerrcnt
      write(6,'("...done with ",i3," errors")') ilocerrcnt
!
!*****************
! testing complex routines
!****************
      ilocerrcnt = 0
      ilocptcnt = 0
      print*, ' '
      write(6,'(" Comparing the function values returned by the")')
      write(6,'("   complex routines against stored values")',advance="NO")
      write(19,'( "***Comparing the function values returned by the complex routines")')
      write(19,'("  against stored values")')
      z = cmplx(1.5662_prd,3.681_prd,prd) ! integration routine
      call airy_ai(z, aiz, daiz, ierr)
      call airy_ai(cmplx(1.5662,3.681),aiasz,daiasz,ierr)
      ffacd = oned !10.0_prd**max(oned,thvd*log(abs(z)))
      write(19,2000) ai, real(z,prd),aimag(z)
      write(19,2001) real(aiasz), aimag(aiasz), real(aiz,prd), aimag(aiz), &
              0.3807338825307106e+00_prd,0.3605171603516821e+00_prd
      f1 = errg(aiz,&
        cmplx(0.3807338825307106e+00_prd,0.3605171603516821e+00_prd,prd))
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,2000) dai, real(z,prd),aimag(z)
      write(19,2001) real(daiasz), aimag(daiasz), real(daiz,prd), aimag(daiz), &
             -0.2685171251031442e+00_prd,-0.1010759493443000e+01_prd
      f1 = errg(daiz,&
        cmplx(-0.2685171251031442e+00_prd,-0.1010759493443000e+01_prd,prd))
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      ilocptcnt = ilocptcnt + 2
      z = cmplx(-9.8342_prd,4.2811_prd,prd) ! asymptotic expansion routine
      call airy_ai(z, aiz, daiz, ierr)
      call airy_ai(cmplx(-9.8342,4.2811),aiasz,daiasz,ierr)
      write(19,2000) ai, real(z,prd),aimag(z)
      write(19,2001) real(aiasz), aimag(aiasz), real(aiz,prd), aimag(aiz), &
             0.1069401631365655e+06_prd,-0.4772398758353830e+05_prd
      f1 = errg(aiz,&
        cmplx(0.1069401631365655e+06_prd,-0.4772398758353830e+05_prd,prd))
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,2000) dai, real(z,prd),aimag(z)
      write(19,2001) real(daiasz), aimag(daiasz), real(daiz,prd), aimag(daiz), &
             -0.2216534220094920e+06_prd,-0.3110801974925891e+06_prd
      f1 = errg(daiz,&
        cmplx(-0.2216534220094920e+06_prd,-0.3110801974925891e+06_prd,prd))
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      ilocptcnt = ilocptcnt + 2
      write(19,1500) ilocptcnt
      write(19,1501) ilocerrcnt
      iptcnt = iptcnt + ilocptcnt
      ierrcnt = ierrcnt + ilocerrcnt
      write(6,'("...done with ",i3," errors")') ilocerrcnt
!
!*****************
! testing interface of the routines
!****************
!  a perturbation analysis on the asymptotic expansion shows that 
!    |Ai(x(1+eps)) - Ai(x)| <= |Ai(x)|*(.25+x**(3/2))*eps.  the 
!    same holds for the derivative.
!    we will use this tolerance
      print*, ' '
      write(6,'(" Comparing the function values on both side of the delineating circle")')
      write(6,'("   and verifying continuity within a tolerance")',advance="NO")
      write(19,'("***Testing function values on both side of the delineating circle")')
      write(19,'("  and verifying continuity within a tolerance")')
      write(19,'("  (no output printed unless an error is found)")')
      ilocerrcnt = 0
      ilocptcnt = 0
! for real variables
      call airy_info(rms)
      call airy_info(rmd)
      call airy_ai(rms+epsilon(rms)*2,aias,daias,ierr)
      call airy_ai(rms-epsilon(rms)*2,aibs,daibs,ierr)
      ffacs = abs(aias)*(.250_prs + rms*sqrt(rms))*4*epsilon(rms)
      g1 = errg(aias,aibs)
      if ( g1 > itol*ffacs ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,3000) ai, rms
         write(19,3001) aias, aibs
      end if 
      ffacs = abs(daias)*(.250_prs + rms*sqrt(rms))*4*epsilon(rms)
      g1 = errg(daias,daibs)
      if ( g1 > itol*ffacs ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,3000) dai, rms
         write(19,3001) daias, daibs
      end if 
      call airy_ai(rmd+epsilon(rmd)*2,aix,daix,ierr)
      call airy_ai(rmd-epsilon(rmd)*2,aiass,daiass,ierr)
      ffacd = abs(aix)*(.250_prd + rmd*sqrt(rmd))*4*epsilon(rmd)
      f1 = errg(aix,aiass)
      if ( f1 > itol*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,3000) ai, rmd
         write(19,3001) aix, aiass
      end if 
      f1 = errg(daix,daiass)
      ffacd = abs(daix)*(.250_prd + rmd*sqrt(rmd))*4*epsilon(rmd)
      if ( f1 > itol*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,3000) dai, rmd
         write(19,3001) daix, daiass
      end if 
      ilocptcnt = ilocptcnt + 4
! for complex variables
      do i = -99,99
        zs = (rms + epsilon(rms)*2)*exp(ciunits*i*0.01_prs*pis)
        call airy_ai(zs,aiasz,daiasz,ierr)
        zs = (rms - epsilon(rms)*2)*exp(ciunits*i*0.01_prs*pis)
        call airy_ai(zs,aibsz,daibsz,ierr)
        ffacs = abs(aiasz)*(0.250_prs+rms*sqrt(rms))*4*epsilon(rms)
        g1 = errg(aiasz,aibsz)
        if ( g1 > itol*ffacs ) then 
          ilocerrcnt = ilocerrcnt+1
          write(19,4000) ai, rms*exp(ciunits*i*0.01_prs*pis)
          write(19,4001) aiasz,aibsz
write(19,*) g1, itol*ffacs
        end if 
        g1 = errg(daiasz,daibsz)
        ffacs = abs(daiasz)*(0.250_prs+rms*sqrt(rms))*4*epsilon(rms)
        if ( g1 > itol*ffacs ) then 
          ilocerrcnt = ilocerrcnt+1
          write(19,*) g1, ffacs
          write(19,4000) dai, rms*exp(ciunits*i*0.01_prs*pis)
          write(19,4001) daiasz,daibsz
        end if 
        z = (rmd + epsilon(rmd)*2)*exp(ciunitd*i*0.01_prd*pid)
        call airy_ai(z,aiz,daiz,ierr)
        z = (rmd - epsilon(rmd)*2)*exp(ciunitd*i*0.01_prd*pid)
        call airy_ai(z,aiassz,daiassz,ierr)
        f1 = errg(aiz,aiassz) 
        ffacd = abs(aiz)*(0.250_prd+rmd*sqrt(rmd))*4*epsilon(rmd)
        if ( f1 > itol*ffacd ) then 
          ilocerrcnt = ilocerrcnt+1
          write(19,4000) ai, rmd*exp(ciunitd*i*0.01_prd*pid)
          write(19,4002) aiasz,aiassz, aiassz
        end if 
        f1 = errg(daiz,daiassz)
        ffacd = abs(daiz)*(0.250_prd+rmd*sqrt(rmd))*4*epsilon(rmd)
        if ( f1 > itol*ffacd ) then 
          ilocerrcnt = ilocerrcnt+1
          write(19,4000) dai, rmd*exp(ciunitd*i*0.01_prd*pid)
          write(19,4002) daiz,daiassz
        end if 
        ilocptcnt = ilocptcnt + 4
      end do
      write(19,1500) ilocptcnt
      write(19,1501) ilocerrcnt
      iptcnt = iptcnt + ilocptcnt
      ierrcnt = ierrcnt + ilocerrcnt
      write(6,'("...done with ",i3," errors")') ilocerrcnt
!
!*****************
! testing different precisions
!****************
      write(6,'(" ")')
      write(6,'(" Comparing the computation in different precisions at random points")')
      write(6,'("   on the real line and in the complex plane and verifying the")')
      write(6,'("   accuracy of the lower precision function values.")',advance='no')
      write(19,'("***Comparing the computation in different precisions at random points ")')
      write(19,'("  on the real line and in the complex plane and verifying the accuracy of")')
      write(19,'("  the lower precision function values (no output printed unless an error is found).")')
      ilocerrcnt = 0
      ilocptcnt = 0
! for real variables
      do i = 1,200
        call random_number(rmd)
        rms = (rmd-0.50_prd)*20.0_prd
        rmd = rms
        call airy_ai(rmd, aix, daix, ierr)
        call airy_ai(rms, aias, daias, ierr)
        f1 = errg(aix,aias)
        ffacs = ones !10.0_prs**max(ones,thvs*log(abs(rms)))
        if ( f1 > itol2*tols*ffacs ) then 
          ilocerrcnt = ilocerrcnt+1
          write(19,5000) ai, rmd
          write(19,5001) aias, aix
        end if 
        f1 = errg(daix,daias)
        if ( f1 > itol2*tols*ffacs ) then
          ilocerrcnt = ilocerrcnt+1
          write(19,5000) dai, rmd
          write(19,5001) daias, daix
        end if 
        ilocptcnt = ilocptcnt + 2
      end do
! for complex variables
      do i = 1,200
        call random_number(rmd)
        call random_number(x)
        zs = rmd*10*exp(ciunitd*(x-0.50_prd)*2.0_prd*pid)
        z = zs
        call airy_ai(zs,aiasz,daiasz,ierr)
        call airy_ai(z,aiz,daiz,ierr)
        ffacs = ones !10.0_prs**max(ones,thvs*log(abs(zs)))
        f1 = errg(aiz,aiasz)
        if ( f1 > itol*tols*ffacs ) then 
          ilocerrcnt = ilocerrcnt+1
          write(19,6000) ai, z
          write(19,6001) aiasz,aiz
        end if 
        f1 = errg(daiz,daiasz)
        if ( f1 > itol*tols*ffacs ) then 
          ilocerrcnt = ilocerrcnt+1
          write(19,6000) dai, z
          write(19,6001) daiasz,daiz
        end if 
        ilocptcnt = ilocptcnt + 2
      end do
      write(19,1500) ilocptcnt
      write(19,1501) ilocerrcnt
      iptcnt = iptcnt + ilocptcnt
      ierrcnt = ierrcnt + ilocerrcnt
      write(6,'("...done with ",i3," errors")') ilocerrcnt
!
!*****************
! testing vectorized airy_ai routine
!****************
      print*, ' '
      write(6,'(" Testing the vectorized routine of the Airy functions")',advance='no')
      write(19,'("***Testing the vectorized routine of the Airy functions (no output printed")')
      write(19,'("  unless an error is found).")')
      ilocerrcnt = 0
      ilocptcnt = 0
      do j = -20,20
        do i = 1,20
          ray(i) = 0.5_prd*i
          ilocptcnt = ilocptcnt + 1
        end do
        call airy_ai(ray,0.05_prd*j*pid,airay,dairay,ierr,.true.)
        do i = 1,20
          call airy_ai(ray(i)*exp(ciunitd*j*0.05_prd*pid), aiz,daiz, ierr)
          g1 = errg (aiz,airay(i))
          ffacd = oned !10.0_prd**max(oned,thvd*log(ray(i)))
          if ( g1 > itol2*told*ffacd ) then 
            ilocerrcnt = ilocerrcnt+1
            write(19,2002) ai, ray(i), j*0.05_prd*pid
            write(19,2003) aiz,  airay(i)
          end if
          g1 = errg (daiz,dairay(i))
          if ( g1 > itol2*told*ffacd ) then 
            ilocerrcnt = ilocerrcnt+1
            write(19,2002) dai, ray(i), j*0.05_prd*pid
            write(19,2003) daiz, dairay(i)
          end if 
        end do
      end do
      write(19,1500) ilocptcnt
      write(19,1501) ilocerrcnt
      iptcnt = iptcnt + ilocptcnt
      ierrcnt = ierrcnt + ilocerrcnt
      write(6,'("...done with ",i3," errors")') ilocerrcnt
!
!*****************
! testing zero routines
!****************
      print*, ' '
      print*, 'Testing the routines for the zeros of the Airy functions'
      write(19,'("***Testing the routines for the zeros of the Airy functions")')
      write(19,'("  real and complex zeros--scalar version")')
      ilocerrcnt = 0
      ilocptcnt = 0
      write(6,'("   real zeros by table lookup")',advance='no')
      write(19,'("  **Zeros by table lookup")')
      call airy_ai_zero (18, aias, ierr, ai_assoc=daias, dai_zero=aibs, &
                          dai_assoc=daibs)
      call airy_ai_zero (18, aix, ierr, ai_assoc=aiass, dai_zero=daix, &
                          dai_assoc=daiass)
      ffacd = 10.0_prd**max(oned,thvd*log(abs(18.0_prd)))
      write(19,1002) 18, ai
      write(19,1001) aias, aix, -19.126380474246954_prd
      f1 = errg(aix, -19.126380474246954_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1003) ai, 18
      write(19,1001) daias, aiass,-1.179880729870146_prd
      f1 = errg(aiass,-1.179880729870146_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1002) 18, dai
      write(19,1001) aibs, daix, -18.764798437665956_prd
      f1 = errg(daix,-18.764798437665956_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1003) dai, 18
      write(19,1001) daibs, daiass,-0.2710702785769711_prd
      f1 = errg(daiass,-0.2710702785769711_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      ilocptcnt = ilocptcnt + 4
      call airy_bi_zero (18, aias, ierr, bi_assoc=daias, dbi_zero=aibs, &
                          dbi_assoc=daibs)
      call airy_bi_zero (18, aix, ierr, bi_assoc=aiass, dbi_zero=daix, &
                          dbi_assoc=daiass)
      write(19,1002) 18, bi
      write(19,1001) aias, aix, -18.765508284480081_prd
      f1 = errg(aix,-18.765508284480081_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1003) bi, 18
      write(19,1001) daias, aiass, -1.174276253118059_prd
      f1 = errg(aiass,-1.174276253118059_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1002) 18, dbi
      write(19,1001) aibs, daix, -19.125697156412638_prd
      f1 = errg(daix,-19.125697156412638_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1003) dbi, 18
      write(19,1001) daibs, daiass,0.2697826139865426_prd
      f1 = errg(daiass,0.2697826139865426_prd)
      if ( f1 > itol*told*ffacd )  ilocerrcnt = ilocerrcnt+1
      ilocptcnt = ilocptcnt + 4
      call airy_bi_zero (18, aiasz, ierr, bi_assoc=daiasz, dbi_zero=aibsz, &
                          dbi_assoc=daibsz)
      call airy_bi_zero (18, aiz, ierr, bi_assoc=aiassz, dbi_zero=daiz, &
                          dbi_assoc=daiassz)
      write(19,1004) 18, bi
      write(19,2001) aiasz, aiz, 9.494603679187176_prd,16.603624606898762_prd
      f1 = errg(aiz,cmplx(9.494603679187176_prd,16.603624606898762_prd,prd))
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1005) bi, 18
      write(19,2001) daiasz, aiassz, 1.445920793639133_prd,-0.8328073286646578_prd
      f1 = errg(aiassz,cmplx(1.445920793639133_prd,-0.8328073286646578_prd,prd))
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1004) 18, dbi
      write(19,2001) aibsz, daiz, 9.313152405593771_prd,16.290870307133684_prd
      f1 = errg(daiz,cmplx(9.313152405593771_prd,16.290870307133684_prd,prd))
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1005) dbi, 18
      write(19,2001) daibsz, daiassz, -0.3321948851433578_prd,-0.1913210621443901_prd
      f1 = errg(daiassz,&
           cmplx(-0.3321948851433578_prd,-0.1913210621443901_prd,prd))
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      ilocptcnt = ilocptcnt + 4
      write(19,1500) ilocptcnt
      write(19,1501) ilocerrcnt
      iptcnt = iptcnt + ilocptcnt
      ierrcnt = ierrcnt + ilocerrcnt
      write(6,'("...................done with ",i3," errors")') ilocerrcnt
!
!*****************
! testing zero routines
!****************
      write(6,'("   real zeros by asymptotic expansion")',advance='no')
      write(19,'("  **Zeros by asymptotic expansion")')
      call airy_ai_zero (86, aias, ierr, ai_assoc=daias, dai_zero=aibs, &
                          dai_assoc=daibs)
      call airy_ai_zero (86, aix, ierr, ai_assoc=aiass, dai_zero=daix, &
                          dai_assoc=daiass)
      ffacd = oned !10.0_prd**max(oned,thvd*log(abs(86.0_prd)))
      write(19,1002) 86, ai
      write(19,1001) aias, aix, -54.657586491868699_prd
      f1 = errg(aix,-54.657586491868699_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1003) ai, 86
      write(19,1001) daias, aiass,-1.534044239238300_prd
      f1 = errg(aiass,-1.534044239238300_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1002) 86, dai
      write(19,1001) aibs, daix, -54.444826792609810_prd
      f1 = errg(daix,-54.444826792609810_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1003) dai, 86
      write(19,1001) daibs, daiass,-0.2076995780278862_prd
      f1 = errg(daiass,-0.2076995780278862_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      ilocptcnt = ilocptcnt + 4
      call airy_bi_zero (86, aias, ierr, bi_assoc=daias, dbi_zero=aibs, &
                          dbi_assoc=daibs)
      call airy_bi_zero (86, aix, ierr, bi_assoc=aiass, dbi_zero=daix, &
                          dbi_assoc=daiass)
      write(19,1002) 86, bi
      write(19,1001) aias, aix, -54.444911130586874_prd
      f1 = errg(aix,-54.444911130586874_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1003) bi, 86
      write(19,1001) daias, aiass,-1.532549805062612_prd
      f1 = errg(aiass,-1.532549805062612_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1002) 86, dbi
      write(19,1001) aibs, daix,-54.657502808936158_prd
      f1 = errg(daix,-54.657502808936158_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1003) dbi, 86
      write(19,1001) daibs, daiass, 0.2074972409267743_prd
      f1 = errg(daiass, 0.2074972409267743_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      ilocptcnt = ilocptcnt + 4
      call airy_bi_zero (86, aiasz, ierr, bi_assoc=daiasz, dbi_zero=aibsz, &
                          dbi_assoc=daibsz)
      call airy_bi_zero (86, aiz, ierr, bi_assoc=aiassz, dbi_zero=daiz, &
                          dbi_assoc=daiassz)
      write(19,1004) 86, bi
      write(19,2001) aiasz, aiz, 27.288200667644368_prd,47.358306153862507_prd
      f1 = errg(aiz,cmplx(27.288200667644368_prd,47.358306153862507_prd,prd))
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1005) bi, 86
      write(19,2001) daiasz, aiassz, 1.879045614304762_prd,-1.084330361798879_prd
      f1 = errg(aiassz,cmplx(1.879045614304762_prd,-1.084330361798879_prd,prd))
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1004) 86, dbi
      write(19,2001) aibsz, daiz, 27.181741517075782_prd,47.174096724925626_prd
      f1 = errg(daiz,cmplx(27.181741517075782_prd,47.174096724925626_prd,prd))
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1005) dbi, 86
      write(19,2001) daibsz, daiassz, -0.2544106266601735_prd,-0.1468108932911459_prd
      f1 = errg(daiassz,cmplx(-0.2544106266601735_prd,-0.1468108932911459_prd,prd))
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      ilocptcnt = ilocptcnt + 4
      write(19,1500) ilocptcnt
      write(19,1501) ilocerrcnt
      iptcnt = iptcnt + ilocptcnt
      ierrcnt = ierrcnt + ilocerrcnt
      write(6,'("...........done with ",i3," errors")') ilocerrcnt
!
!*****************
! testing vectorized zero routines
!****************
      write(6,'("   real and complex zeros--vectorized version")',advance="no")
      write(19,'("  **Zeros--vectorized version (no output printed unless an error is found)")')
      ilocerrcnt = 0
      ilocptcnt = 0
      call airy_ai_zero (23, aizv, ierr, ai_assoc=aiassv)
      ffacd = oned !10.0_prd**max(oned,thvd*log(abs(27.0_prd)))
      f1 = errg(aizv(1),-22.567612917496504_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 23, ai
         write(19,1008) aizv(1), -22.567612917496504_prd
      end if 
      f1 = errg(aizv(2),-23.224165001121680_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 24, ai
         write(19,1008) aizv(2), -23.224165001121680_prd
      end if 
      f1 = errg(aizv(3),-23.871564455535918_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 25, ai
         write(19,1008) aizv(3), -23.871564455535918_prd
      end if 
      f1 = errg(aizv(4),-24.510301236589672_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 26, ai
         write(19,1008) aizv(4), -24.510301236589672_prd
      end if 
      f1 = errg(aizv(5),-25.140821166148957_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 27, ai
         write(19,1008) aizv(5), -25.140821166148957_prd
      end if 
      f1 = errg(aiassv(1),1.229700701509681_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) ai, 23
         write(19,1008) aiassv(1),1.229700701509681_prd
      end if 
      f1 = errg(aiassv(2),-1.238547875329632_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) ai, 24
         write(19,1008) aiassv(2),-1.238547875329632_prd
      end if 
      f1 = errg(aiassv(3),1.247089945259408_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) ai, 25
         write(19,1008) aiassv(3),1.247089945259408_prd
      end if 
      f1 = errg(aiassv(4),-1.255349140475735_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) ai, 26
         write(19,1008) aiassv(4),-1.255349140475735_prd
      end if 
      f1 = errg(aiassv(5),1.263345282750799_prd)
      if ( f1 > itol*told*ffacd )then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) ai, 27
         write(19,1008) aiassv(5),1.263345282750799_prd
      end if 
      ilocptcnt = ilocptcnt + 10
      call airy_ai_zero (23, aizv, ierr, dai_zero=daizv, dai_assoc=daiassv)
      f1 = errg(daizv(1),-22.235232285348914_prd)
      if ( f1 > itol*told*ffacd )then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 23, dai
         write(19,1008) daizv(1),-22.235232285348914_prd
      end if 
      f1 = errg(daizv(2),-22.896588738874620_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 24, dai
         write(19,1008) daizv(2),-22.896588738874620_prd
      end if 
      f1 = errg(daizv(3),-23.548526295928802_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 25, dai
         write(19,1008) daizv(3),-23.548526295928802_prd
      end if 
      f1 = errg(daizv(4),-24.191559709526349_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 26, dai
         write(19,1008) daizv(4),-24.191559709526349_prd
      end if 
      f1 = errg(daizv(5),-24.826156425921152_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 27, dai
         write(19,1008) daizv(5),-24.826156425921152_prd
      end if 
      f1 = errg(daiassv(1),0.259812670151466_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) dai, 23
         write(19,1008) daiassv(1),0.259812670151466_prd
      end if 
      f1 = errg(daiassv(2),-0.257916075332572_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) dai, 24
         write(19,1008) daiassv(2),-0.257916075332572_prd
      end if 
      f1 = errg(daiassv(3),0.256112333779654_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) dai, 25
         write(19,1008) daiassv(3),0.256112333779654_prd
      end if 
      f1 = errg(daiassv(4),-0.254393342646825_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) dai, 26
         write(19,1008) daiassv(4),-0.254393342646825_prd
      end if 
      f1 = errg(daiassv(5),0.252751992576574_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) dai, 27
         write(19,1008) daiassv(5),0.252751992576574_prd
      end if 
      ilocptcnt = ilocptcnt + 10
!
!*****************
! testing vectorized zero routines
!****************
      call airy_bi_zero (23, aizv, ierr, bi_assoc=aiassv)
      f1 = errg(aizv(1),-22.235737881803384_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 23, bi
         write(19,1008) aizv(1),-22.235737881803384_prd
      end if 
      f1 = errg(aizv(2),-22.897065554219793_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 24, bi
         write(19,1008) aizv(2),-22.897065554219793_prd
      end if 
      f1 = errg(aizv(3),-23.548977079642448_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 25, bi
         write(19,1008) aizv(3),-23.548977079642448_prd
      end if 
      f1 = errg(aizv(4),-24.191986850648995_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 26, bi
         write(19,1008) aizv(4),-24.191986850648995_prd
      end if 
      f1 = errg(aizv(5),-24.826562012152888_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 27, bi
         write(19,1008) aizv(5),-24.826562012152888_prd
      end if 
      f1 = errg(aiassv(1),1.225154995846960_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) bi, 23
         write(19,1008) aiassv(1),1.225154995846960_prd
      end if 
      f1 = errg(aiassv(2),-1.234163920487302_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) bi, 24
         write(19,1008) aiassv(2),-1.234163920487302_prd
      end if 
      f1 = errg(aiassv(3),1.242855598092317_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) bi, 25
         write(19,1008) aiassv(3),1.242855598092317_prd
      end if 
      f1 = errg(aiassv(4),-1.251253611221631_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) bi, 26
         write(19,1008) aiassv(4),-1.251253611221631_prd
      end if 
      f1 = errg(aiassv(5),1.259378938688538_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) bi, 27
         write(19,1008) aiassv(5),1.259378938688538_prd
      end if 
      ilocptcnt = ilocptcnt + 10
      call airy_bi_zero (23, aizv, ierr, dbi_zero=daizv, dbi_assoc=daiassv)
      f1 = errg(daizv(1),-22.567122080497199_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 23, dbi
         write(19,1008) daizv(1),-22.567122080497199_prd
      end if 
      f1 = errg(daizv(2),-23.223701521208962_prd)
      if ( f1 > itol*told*ffacd )then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 24, dbi
         write(19,1008) daizv(2),-23.223701521208962_prd
      end if 
      f1 = errg(daizv(3),-23.871125771677974_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 25, dbi
         write(19,1008) daizv(3),-23.871125771677974_prd
      end if 
      f1 = errg(daizv(4),-24.509885117016239_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 26, dbi
         write(19,1008) daizv(4),-24.509885117016239_prd
      end if 
      f1 = errg(daizv(5),-25.140425655367874_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 27, dbi
         write(19,1008) daizv(5),-25.140425655367874_prd
      end if 
      f1 = errg(daiassv(1),-0.258852215916922_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) dbi, 23
         write(19,1008) daiassv(1),-0.258852215916922_prd
      end if 
      f1 = errg(daiassv(2),0.257003129647316_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) dbi, 24
         write(19,1008) daiassv(2),0.257003129647316_prd
      end if 
      f1 = errg(daiassv(3),-0.255242710064815_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) dbi, 25
         write(19,1008) daiassv(3),-0.255242710064815_prd
      end if 
      f1 = errg(daiassv(4),0.253563372437697_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) dbi, 26
         write(19,1008) daiassv(4),0.253563372437697_prd
      end if 
      f1 = errg(daiassv(5),-0.251958444330784_prd)
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) dbi, 27
         write(19,1008) daiassv(5),-0.251958444330784_prd
      end if 
      ilocptcnt = ilocptcnt + 10
!
!*****************
! testing vectorized zero routines
!****************
      call airy_bi_zero (23, bizv, ierr, bi_assoc=biassv)
      f1 = errg(bizv(1),&
        cmplx(11.220656371214879_prd,19.580653883038824_prd,prd))
      if ( f1 > itol*told*ffacd )  then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1004) 23, bi
         write(19,2003) bizv(1),11.220656371214879_prd,19.580653883038824_prd
      end if 
      f1 = errg(bizv(2),&
        cmplx(11.549830143889196_prd,20.148722567309367_prd,prd))
      if ( f1 > itol*told*ffacd )  then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1004) 24, bi
         write(19,2003) bizv(2), 11.549830143889196_prd,20.148722567309367_prd
      end if 
      f1 = errg(bizv(3),&
        cmplx(11.874378641347221_prd,20.708893464731311_prd,prd))
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1004) 25, bi
         write(19,2003) bizv(3),11.874378641347221_prd,20.708893464731311_prd
      end if 
      f1 = errg(bizv(4),&
        cmplx(12.194551330833875_prd,21.261588251953778_prd,prd))
      if ( f1 > itol*told*ffacd )  then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1004) 26, bi
         write(19,2003) bizv(4),12.194551330833875_prd,21.261588251953778_prd
      end if 
      f1 = errg(bizv(5),&
        cmplx(12.510575046777394_prd,21.807190718498749_prd,prd))
      if ( f1 > itol*told*ffacd )  then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1004) 27, bi
         write(19,2003) bizv(5),12.510575046777394_prd,21.807190718498749_prd
      end if 
      f1 = errg(biassv(1),&
        cmplx(-1.506774749698633_prd,0.868314074899650_prd,prd))
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1005) bi, 23
         write(19,2003) biassv(1),-1.506774749698633_prd,0.868314074899650_prd
      end if 
      f1 = errg(biassv(2),&
        cmplx(1.517585357164310_prd,-0.874612709515018_prd,prd))
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1005) bi, 24
         write(19,2003) biassv(2),1.517585357164310_prd,-0.874612709515018_prd
      end if 
      f1 = errg(biassv(3),&
        cmplx(-1.528024149210944_prd,0.880692431621191_prd,prd))
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1005) bi, 25
         write(19,2003) biassv(3),-1.528024149210944_prd,0.880692431621191_prd
      end if 
      f1 = errg(biassv(4),&
        cmplx(1.538118145012280_prd,-0.886569309158123_prd,prd))
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1005) bi, 26
         write(19,2003) biassv(4),1.538118145012280_prd,-0.886569309158123_prd
      end if 
      f1 = errg(biassv(5),&
        cmplx(-1.547891444975066_prd,0.892257657608977_prd,prd))
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1005) bi, 27
         write(19,2003) biassv(5),-1.547891444975066_prd,0.892257657608977_prd
      end if 
      ilocptcnt = ilocptcnt + 10
      call airy_bi_zero (23, bizv, ierr, dbi_zero=dbizv, dbi_assoc=dbiassv)
      f1 = errg(dbizv(1),&
        cmplx(11.053994360667899_prd,19.293078213959429_prd,prd))
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 23, dbi
         write(19,2003) dbizv(1),11.053994360667899_prd,19.293078213959429_prd
      end if 
      f1 = errg(dbizv(2),&
        cmplx(11.385596969302076_prd,19.865292017137630_prd,prd))
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 24, dbi
         write(19,2003) dbizv(2),11.385596969302076_prd,19.865292017137630_prd
      end if 
      f1 = errg(dbizv(3),&
        cmplx(11.712438641126221_prd,20.429378925421950_prd,prd))
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 25, dbi
         write(19,2003) dbizv(3),11.712438641126221_prd,20.429378925421950_prd
      end if 
      f1 = errg(dbizv(4),&
        cmplx(12.034781571133943_prd,20.985781895965161_prd,prd))
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 26, dbi
         write(19,2003) dbizv(4),12.034781571133943_prd,20.985781895965161_prd
      end if 
      f1 = errg(dbizv(5),&
        cmplx(12.352863681490561_prd,21.534903282918666_prd,prd))
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1002) 27, dbi
         write(19,2003) dbizv(5),12.352863681490561_prd,21.534903282918666_prd
      end if 
      f1 = errg(dbiassv(1),&
        cmplx(0.318355274563740_prd,0.183451937317410_prd,prd))
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) dbi, 23
         write(19,2003) dbiassv(1),0.318355274563740_prd,0.183451937317410_prd
      end if 
      f1 = errg(dbiassv(2),&
        cmplx(-0.316024910474316_prd,-0.182124025592498_prd,prd))
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) dbi, 24
         write(19,2003) dbiassv(2),-0.316024910474316_prd,-0.182124025592498_prd
      end if 
      f1 = errg(dbiassv(3),&
        cmplx(0.313808934629479_prd,0.180860595992543_prd,prd))
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) dbi, 25
         write(19,2003) dbiassv(3),0.313808934629479_prd,0.180860595992543_prd
      end if 
      f1 = errg(dbiassv(4),&
        cmplx(-0.311697340122842_prd,-0.179656065952459_prd,prd))
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) dbi, 26
         write(19,2003) dbiassv(4),-0.311697340122842_prd,-0.179656065952459_prd
      end if 
      f1 = errg(dbiassv(5),&
        cmplx(0.309681349894832_prd,0.178505532129885_prd,prd))
      if ( f1 > itol*told*ffacd ) then 
         ilocerrcnt = ilocerrcnt+1
         write(19,1003) dbi, 27
         write(19,2003) dbiassv(5),0.309681349894832_prd,0.178505532129885_prd
      end if 
      ilocptcnt = ilocptcnt + 10
      write(19,1500) ilocptcnt
      write(19,1501) ilocerrcnt
      iptcnt = iptcnt + ilocptcnt
      ierrcnt = ierrcnt + ilocerrcnt
      write(6,'("...done with ",i3," errors")') ilocerrcnt
!
!*****************
! testing modulus and phase routines
!****************
      print*, ' '
      write(6,'(" Comparing the computed modulus, phase, and associated values")')
      write(6,'("   against stored values")',advance='no')
      write(19,'("***Comparing the computed modulus, phase, and associated values against stored values")')
      ilocerrcnt = 0
      ilocptcnt = 0
      x = -5.234_prd
      write(19,'("  **by ODE integration")')
      call airy_aux ( real(x,prs), aias, daias, ierr, aibs, daibs)
      call airy_aux ( x, aix, daix, ierr, aiass, daiass)
      ffacd = oned !10.0_prd**max(oned,thvd*log(abs(x)))
      write(19,1006) 'M(x)', x
      write(19,1001) aias, aix,0.3728081698402362_prd
      f1 = errg(aix,0.3728081698402362_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1007) 'theta(x)', x
      write(19,1001) daias, daix,8.759640547244738_prd
      f1 = errg(daix,8.759640547244738_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1006) 'M(x)', x
      write(19,1001) aibs, aiass,0.8540001773409781_prd
      f1 = errg(aiass,0.8540001773409781_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1007) 'phi(x)', x
      write(19,1001) daibs, daiass,7.209566736849307_prd
      f1 = errg(daiass,7.209566736849307_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      ilocptcnt = ilocptcnt + 4
!
!*****************
! testing modulus and phase routines
!****************
      x = -25.392_prd
      write(19,'("  **by asymptotic expansion")')
      call airy_aux ( real(x,prs), aias, daias, ierr, aibs, daibs)
      call airy_aux ( x, aix, daix, ierr, aiass, daiass)
      ffacd = oned !10.0_prd**max(oned,thvd*log(abs(x)))
      write(19,1006) 'M(x)', x
      write(19,1001) aias, aix,0.2513325654640693_prd
      f1 = errg(aix,0.2513325654640693_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1007) 'theta(x)', x
      write(19,1001) daias, daix,86.085580681739742_prd
      f1 = errg(daix,86.085580681739742_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1006) 'N(x)', x
      write(19,1001) aibs, aiass,1.2664912447881540_prd
      f1 = errg(aiass,1.2664912447881540_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      write(19,1007) 'phi(x)', x
      write(19,1001) daibs, daiass,84.516738087389200_prd
      f1 = errg(daiass,84.516738087389200_prd)
      if ( f1 > itol*told*ffacd ) ilocerrcnt = ilocerrcnt+1
      ilocptcnt = ilocptcnt + 4
!*****************
! testing modulus and phase in different precisions
!****************
      write(19,'("  **in different precisions (no output printed unless an error is found)")')
      do i = 1,200
        call random_number(rmd)
        rms = -rmd*10.0_prd
        rmd = rms
        call airy_aux ( rmd, aix, daix, ierr, aiass, daiass)
        call airy_aux ( rms, aias, daias, ierr, aibs, daibs)
        ffacd = ones !10.0_prs**max(ones,thvs*log(abs(rms)))
        f1 = errg(aix,aias)
        if ( f1 > itol*tols*ffacs ) then 
          ilocerrcnt = ilocerrcnt+1
         write(19,5002) ai, rmd
         write(19,5001) aias,aix
        end if 
        f1 = errg(daix,daias)
        if ( f1 > itol*tols*ffacs ) then 
          ilocerrcnt = ilocerrcnt+1
         write(19,5002) ai, rmd
         write(19,5001) daias,daix
        end if 
        f1 = errg(aiass,aibs)
        if ( f1 > itol*tols*ffacs ) then 
          ilocerrcnt = ilocerrcnt+1
         write(19,5003) ai, rmd
         write(19,5001) aibs,aiass
        end if 
        f1 = errg(daiass,daibs)
        if ( f1 > itol*tols*ffacs ) then 
          ilocerrcnt = ilocerrcnt+11
         write(19,5003) ai, rmd
         write(19,5001) daibs,daiass
        end if 
        ilocptcnt = ilocptcnt + 4
      end do
      write(19,1500) ilocptcnt
      write(19,1501) ilocerrcnt
      iptcnt = iptcnt + ilocptcnt
      ierrcnt = ierrcnt + ilocerrcnt
      write(6,'("...done with ",i3," errors")') ilocerrcnt    
      write(19,'("Number of tests performed:",i6)') iptcnt
      write(19,'("Number of errors found:",i6)') ierrcnt
      print*, ' '
      print*, 'Number of tests performed:',iptcnt
      print*, 'Number of errors found:',ierrcnt
      print*, ' '
      print*, ' '
      print*, 'End of run'
      print*, ' '
      print*, ' '
      close(21)      
1500 format('  The number of points tested is ',i4)
1501 format('  The number of errors found is  ',i4//)
1000 format(2x,'Computing',2x,A3,2x,'at x=',f9.5)
1001 format(17x,'We computed ',e15.8,' in single precision'/,17x,'We computed ', &
                  e23.16,' in double precision',/,5x,'and the stored value is ',e23.16/)
1002 format('  Testing the ',i2,'th zero of ',A3)
1008 format(17x,'We computed ',e23.16,' in double precision',/,5x,'and the stored value is ',e23.16/)
1003 format('  Testing the associated function value of ',A3,' at the ',i2,' zero')
1004 format('  Testing the ',i2,'th complex zero of ',A3)
1005 format('  Testing the associated function value of ',A3,' at the ',i2,' complex zero')
1006 format('  Testing the modulus ',A4,' at x= ',f9.5)
1007 format('  Testing the phase ',A8,' at x= ',f9.5)
2000 format(2x,'Computing',2x,A3,2x,'at z=',f9.5,' +I*(',f9.5,')')
2001 format(2x,'In single precision, we computed ',e15.8,8x,' +I*(',e15.8,')',/,2x, &
                  'In double precision, we computed ',e23.16,' +I*(',e23.16,')',/,11x, &
                  'and the stored value is ',e23.16,' +I*(',e23.16,')',/)
2002 format('Error encountered testing'/'   the vectorized version of ',2x,A3, &
                   2x,'at (r,theta)= (',f9.5,',',f9.5,')')
2003 format(2x,'In double precision, we computed ',e23.16,' +I*(',e23.16,')',/,11x, &
                  'and the stored value is ',e23.16,' +I*(',e23.16,')',/)
3000 format('Error encountered testing each side of ',2x,A3,2x,'at x=',f9.5)
3001 format(12x,'Below we got ',e23.16,/8x,'and above we got ',e23.16/)
4000 format('Error encountered testing each side of ',2x,A3,2x,'at z=',f9.5, &
             ' +I*(',f9.5,')')
4001 format(12x,'Below we got ',e15.8,' +I*(',e15.8,')',/&
          8x,'and above we got ',e15.8,' +I*(',e15.8,')'/)
4002 format(12x,'Below we got ',e23.16,' +I*(',e23.16,')',/&
          8x,'and above we got ',e23.16,' +I*(',e23.16,')'/)
5000 format('Error encountered testing',2x,A3,2x,'at x=',f9.5)
5001 format(9x,'Single returned ',e15.8,/5x,'and double returned ',e23.16/)
5002 format('Error encountered testing the modulus of ',2x,A3,2x,'at x=',f9.5)
5003 format('Error encountered testing the phase of   ',2x,A3,2x,'at x=',f9.5)
6000 format('Error encountered testing',2x,A3,2x,'at z=',f9.5, &
             ' +I*(',f9.5,')')
6001 format(8x,'Single returned ',e15.8,9x,' +I*(',e15.8,')',/&
          4x,'and double returned ',e23.16,' +I*(',e23.16,')'/)
      end program test


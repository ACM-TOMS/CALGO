*************************************************************************
**                 TEST PROGRAMME FOR THE FGMRES CODE
*************************************************************************

      program validation
*
      integer lda, ldstrt, lwork
      parameter (lda = 1000, ldstrt = 60)
      parameter (lwork = ldstrt**2 + ldstrt*(2*lda+6) + 6*lda + 1)
*
      integer i, j, n, m, m2
      integer revcom, colx, coly, colz, nbscal
      integer revcom2, colx2, coly2, colz2, nbscal2
      integer irc(7), icntl(7), info(3)
      integer irc2(5), icntl2(8), info2(3)
*
      integer matvec, precondLeft, precondRight, dotProd
      parameter (matvec=1, precondLeft=2, precondRight=3, dotProd=4)
*
      integer nout
*
      complex  a(lda,lda), work(lwork)
      real  cntl(5), rinfo, rn, rx, rc
      real  cntl2(5), rinfo2(2)
*
      complex ZERO, ONE
      parameter (ZERO = (0.0e0, 0.0e0), ONE = (1.0e0, 0.0e0))
*
***************************************************************
** Generate the test matrix a and set the right-hand side
** in positions (n+1) to 2n of the array work.
** The right-hand side is chosen such that the exact solution
** is the vector of all ones.
***************************************************************
*
      write(*,*) '***********************************************'
      write(*,*) 'This code is an example of use of FGMRES'
      write(*,*) 'in single precision complex arithmetic'
      write(*,*) 'Results are written in output files'
      write(*,*) 'fort.20 : log file of FGMRES iterations '
      write(*,*) 'fort.30 : log file of inner GMRES iterations '
      write(*,*) 'and sol_Testfgmres : output of the computation.'
      write(*,*) '***********************************************'
      write(*,*)
      write(*,*) 'Matrix size < ', lda
      read(*,*) n
      if (n.gt.lda) then
        write(*,*) 'You are asking for a too large matrix'
        goto 100
      endif
*
      do j = 1,n
        do i = 1,n
          a(i,j) = ZERO
        enddo
        work(j) = ONE
      enddo
*
      do i = 1,n
        a(i,i) = cmplx(1.0+mod(i,10), 0.0)
      enddo
      do i = 1,n-1
        a(i,i+1) = (-3.0, 1.0)
        a(i+1,i) = (-2.0, 1.0)
      enddo
*
      call CGEMV('N',n,n,ONE,A,lda,work(1),1,ZERO,work(n+1),1)
*
      do j = 1,n
        work(j) = ONE/2.0
      enddo
*
*********************************
** Choose the restart parameter
*********************************
*
      write(*,*) 'Restart  <', ldstrt
      read(*,*) m
*
*******************************************************
** Initialize the control parameters to default value
*******************************************************
*
      call init_cfgmres(icntl,cntl)
      call init_cgmres(icntl2,cntl2)
*
*************************
*c Tune some parameters for FGMRES
*************************
*
* Tolerance
      cntl(1) = 1.e-5
      icntl(3) = 20
*  orthogonalization
      print *,' FGMRES orthogonalization 0:MGS, 1:IMGS, 2:CGS, 3:ICGS'
      read(*,*) icntl(4)
      print *,' FGMRES initial guess 0: zero, 1: user supplied '
      read(*,*) icntl(5)
* Maximum number of iterations
      icntl(6) = 100 
      print *,' Init residual at FGMRES restart 0:implicit, 1:explicit'
      read(*,*) icntl(7)
*
*************************
*c Tune some parameters for GMRES
*************************
*
* Tolerance
      cntl2(1) = 1.e-2
* warning output stream
      icntl2(2) = 0
* Save the convergence history in file fort.30
      icntl2(3) = 30
* Left preconditioning
      icntl2(4) = 1
c     print *,' Inner GMRES precond 0-none, 1: left, 2: right '
* ICGS orthogonalization
      icntl2(5) = 3
*
*****************************************
** Reverse communication implementation
*****************************************
*
10     call drive_cfgmres(n,n,m,lwork,work,
     &         irc,icntl,cntl,info,rinfo)
       revcom = irc(1)
       colx   = irc(2)
       coly   = irc(3)
       colz   = irc(4)
       nbscal = irc(5)
*
       if (revcom.eq.matvec) then
* perform the matrix vector product for the FGMRES iteration
*        work(colz) <-- A * work(colx)
         call cgemv('N',n,n,ONE,a,lda,work(colx),1,
     &            ZERO,work(colz),1)
         goto 10
*
       else if (revcom.eq.precondRight) then
* perform the right preconditioning for the FGMRES iteration
*
* Check if there is enough space left in the workspace to perfrom
* few steps of GMRES as right preconditioner
         rn = float(n)
         rx         = rn + 5.0
         rc         = 5.0*rn + 1 - float(irc(7))
*
* Update the linear part of the second order equation to be solved
         if ((icntl2(5).eq.2).or.(icntl2(5).eq.3)) then
           rx = rx + 1
         endif
* Update the constant part of the second order equation to be solved
*             
         if (icntl2(8).eq.0) then
           rc = rc + rn
         endif
         m2 = ifix((-rx+sqrt(rx**2-4.0*rc))/2.0)
* Perform at most two restarts in the inner GMRES
         icntl2(7) = 2*m2
*
         if (m2.gt.0) then
* copy colx in the workspace (right hand side location) of the inner
* gmres iteration 
           call ccopy(n,work(colx),1,work(irc(6)+n),1)
 20        call drive_cgmres(n,n,m2,irc(7),
     &             work(irc(6)),irc2,icntl2,cntl2,info2,rinfo2) 
*
           revcom2 = irc2(1)
           colx2   = irc2(2) + irc(6) -1
           coly2   = irc2(3) + irc(6) -1
           colz2   = irc2(4) + irc(6) -1
           nbscal2 = irc2(5)
           if (revcom2.eq.matvec) then
* Perform the matrix vector product for the inner GMRES iteration
             call cgemv('N',n,n,ONE,a,lda,work(colx2),1,
     &                 ZERO,work(colz2),1)
             goto 20
           else if (revcom2.eq.precondRight) then
* perform the preconditioning for the inner GMRES iteration
             do i =0,n-1
               work(colz2+i) = work(colx2+i)/a(i+1,i+1)  
             enddo
             goto 20
           else if (revcom2.eq.precondleft) then
* perform the preconditioning for the inner GMRES iteration
             do i =0,n-1
               work(colz2+i) = work(colx2+i)/a(i+1,i+1)  
             enddo
             goto 20
           else if (revcom2.eq.dotProd) then
*      perform the scalar product for the inner GMRES iteration
*      work(colz) <-- work(colx) work(coly)
             call cgemv('C',n,nbscal2,ONE,work(colx2),n,
     &               work(coly2),1,ZERO,work(colz2),1)
             goto 20
           endif
           call ccopy(n,work(irc(6)),1,work(colz),1)
           goto 10
         else
* (m2.le.0)
           print *, ' Not enough space for inner gmres '
           call ccopy(n,work(colx),1,work(colz),1)
           goto 10
         endif
       else if (revcom.eq.dotProd) then
*      perform the scalar product for the FGMRES iteration
*      work(colz) <-- work(colx) work(coly)
*
         call cgemv('C',n,nbscal,ONE,work(colx),n,
     &               work(coly),1,ZERO,work(colz),1)
         goto 10
       endif
*
*******************************
* dump the solution on a file
*******************************
*
      nout = 11
      open(nout,FILE='sol_Testfgmres',STATUS='unknown')
      if (icntl(5).eq.0) then
        write(nout,*) 'Orthogonalisation : MGS'
      elseif (icntl(5).eq.1) then
        write(nout,*) 'Orthogonalisation : IMGS'
      elseif (icntl(5).eq.2) then
        write(nout,*) 'Orthogonalisation : CGS'
      elseif (icntl(5).eq.3) then
        write(nout,*) 'Orthogonalisation : ICGS'
      endif
      write(nout,*) 'Restart : ', m
      write(nout,*) 'info(1) = ',info(1),'  info(2) = ',info(2)
      write(nout,*) 'rinfo = ',rinfo
      write(nout,*) 'Optimal workspace = ', info(3)
      write(nout,*) 'Solution : '
      do j=1,n
        write(nout,*) work(j)
      enddo
      write(nout,*) '   '
*
100    continue
*
      stop
      end

*************************************************************************
**                 TEST PROGRAMME FOR THE zfgmres CODE
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
      double complex  a(lda,lda+2), work(lwork), aux(lda+2)
      double precision  cntl(5), rinfo
      double precision  cntl2(5), rinfo2(2)
      real rn, rx, rc
*
      double complex ZERO, ONE
      parameter (ZERO = (0.0d0, 0.0d0), ONE = (1.0d0, 0.0d0))
*
* variables required by the parallel implementation
      include 'mpif.h'
      integer type, token, status(MPI_STATUS_SIZE)
      integer nproc, infompi, comm, me, nloc, iconf(2)
      integer nbcol, istart, jstart
      integer comm_inner, comm_outer
*
* MPI initialization
*
      call MPI_INIT(infompi)
      comm = MPI_COMM_WORLD
      call MPI_COMM_SIZE(comm,nproc,infompi)
      call MPI_COMM_RANK(comm,me,infompi)
*
* Duplicate the communicator to discreminate beetween inner and outer loop
*  exchanges
      call MPI_COMM_DUP(comm, comm_inner, infompi)
      call MPI_COMM_DUP(comm, comm_outer, infompi)
*
***************************************************************
** Generate the test matrix a and set the right-hand side
** in positions (n+1) to 2n of the array work.
** The right-hand side is chosen such that the exact solution
** is the vector of all ones.
***************************************************************
* The solution of the tridiagonal system is performed in parallel
* The matrix is decomposed by block of rows so that the matrix-vector
* can be easily performed.
* Each processor is in charged of a block of rows and stored the
* corresponding entries of the initial guess and the right hand-sides.
* We give below an example of the data distribution for a system of
* dimension 8 in complex arithmetic solved on 2 processors.
*
*                  A                             x  =  b
*
*
*       |  4  -2                                |   x1    b1 
*  P0   | -1+i 4  -2+i                          |   x2    b2  
*       |     -1+i   4  -2+i                    |   x3    b3      
*       |         -1+i   4  -2+i                |   x4    b4          
*  --------------------------------------------------------                 
*       |             -1+i   4  -2+i            |   x5    b5              
*       |                 -1+i   4  -2+i        |   x6    b6                
*  P1   |                     -1+i   4  -2+i    |   x7    b7               
*       |                         -1+i   4  -2+i|   x8    b8                          
*
*
* For the sake of simplicity each processor will have the same number of
* row denoted nloc. Consequently the size of the linear systems will be
* nloc times the number of processors
*
*
      if (me.eq.0) then
         write(*,*) '***********************************************'
         write(*,*) 'This code is an example of use of FGMRES'
         write(*,*) 'in double precision complex arithmetic'
         write(*,*) 'Results are written in output files'
         write(*,*)  'fort.41', ': log file of FGMRES iterations '
         write(*,*)  'fort.141',' : log file of inner
     & zgmres iterations '
         write(*,*) 'and ','sol_zTest','  : output of the computation.'
         write(*,*) '***********************************************'
         write(*,*)
        write(*,*) 'Local matrix size < ', lda
        read(*,*) nloc
        if (nloc.gt.lda) then
          write(*,*) 'You are asking for a too large matrix'
          goto 100
        endif
        write(*,*) ' Global matrix size ',nloc*nproc
*********************************
** Choose the restart parameter
*********************************
*
        write(*,*) 'Restart  <', ldstrt
        read(*,*) m
        iconf(1) = nloc
        iconf(2) = m
      endif
      call MPI_BCAST(iconf,2,MPI_INTEGER,0,comm,infompi)
      nloc = iconf(1)
      m    = iconf(2)
      n    = nloc*nproc
*
* Initialize the local matrix (nloc x (nloc+2) ) matrix
* only part of it might be used by the different processor
* depending on its rank
       do j = 1,nloc+2
          do i = 1,nloc
            a(i,j) = ZERO
          enddo
       enddo
*
       do i = 1,nloc
          a(i,i)   = (-1.0d0, 1.0d0)
          a(i,i+1) = (4.0d0,  0.0d0)
          a(i,i+2) = (-2.0d0, 1.0d0)
       enddo
*
*
* Intialise the column index of the first column 
* of the submatrix that will be involved in the mat-vec depending
* on the processor rank.
* Similarly initialize the index of the first entry of the local
* vector in the vector to be involved in the mat-vec.
*
*       jstart = 1 for all the processors but the first
*         ||       that does not have predecessor
*         \/
*                                           x1
*       |-1+i   4  -2+i             |       x2  <-- istart=2 for all the processors
*   A = |    -1+i   4  -2+i         |   x = x3               but the first
*       |        -1+i   4  -2+i     |       x4
*       |            -1+i   4  -2+i |       x5    
*                                           x6
*
       if (me.eq.0) then
        jstart = 2
        istart = 1
       else
        jstart = 1
        istart = 2
       endif
       nbCol = nloc+2
       if (me.eq.0) then
         nbCol = nbCol - 1
       endif
       if (me.eq.(nproc-1)) then
         nbCol = nbCol - 1
       endif
*
** Initialise the right hand side
      do j = 1,nloc+2
        aux(j) = ONE
      enddo
      call zgemv('N',nloc,nbCol,ONE,A(1,jstart),lda,aux,1,
     &            ZERO,work(nloc+1),1)
      do j = 1,nloc
        work(j) = ONE/2.0
      enddo
*
*
*******************************************************
** Initialize the control parameters to default value
*******************************************************
*
      call init_zfgmres(icntl,cntl)
      call init_zgmres(icntl2,cntl2)
*
*************************
*c Tune some parameters for zfgmres
*************************
*
* Save the convergence history standard output
      cntl(1) = 1.d-11
      if (me.eq.0) then
* orthogonalization
       print *,' FGMRES orthogonalization 0:MGS, 1:IMGS, 2:CGS, 3:ICGS'
       read(*,*) icntl(4)
       print *,' FGMRES initial guess 0: zero, 1: user supplied '
       read(*,*) icntl(5)
       print *,' Init residual at FGMRES restart 0:implicit, 1:explicit'
       read(*,*) icntl(7)
      endif
*
      call MPI_BCAST(icntl,7,MPI_INTEGER,0,comm,infompi)
      if (me.eq.0) then
       icntl(3) = 6
* orthogonalization
      else
       icntl(1) = 0
       icntl(2) = 0
       icntl(3) = 0
      endif
* Maximum number of iterations
      icntl(6) = 100
*
*************************
*c Tune some parameters for zgmres
*************************
*
* Tolerance
      cntl2(1) = 1.d-2
* warning output stream
      icntl2(2) = 0
      if (me.eq.0) then
* Save the convergence history in file 'fort.140'
        icntl2(3) = 141
      else
        icntl2(3) = 0
      endif
* Left preconditioning
      icntl2(4) = 1
* ICGS orthogonalization
      icntl2(5) = 3
*****************************************
** Reverse communication implementation
*****************************************
*
10     call drive_zfgmres(n,nloc,m,lwork,work,
     &         irc,icntl,cntl,info,rinfo)
       revcom = irc(1)
       colx   = irc(2)
       coly   = irc(3)
       colz   = irc(4)
       nbscal = irc(5)
*
       if (revcom.eq.matvec) then
* perform the matrix vector product
         call zcopy(nloc,work(colx),1,aux(istart),1)
* Send the entry of aux required to perform the parallel tridiagonal matrix-vector
* product
         if (me.ne.(nproc-1)) then
*  send the last entry of y local vector to the next processor that
*   needs this entry to perform its part of the matrix-vector product
            call MPI_SEND(aux(istart+nloc-1),1,MPI_DOUBLE_COMPLEX,
     &                me+1,2,comm_outer,infompi)
         endif
         if (me.ne.0) then
*  send the first entry of y local vector to the previous processor that
*   needs this entry to perform its part of the matrix-vector product
            call MPI_SEND(aux(istart),1,MPI_DOUBLE_COMPLEX,me-1,
     &                3,comm_outer,infompi)
         endif
* Receive the entry of aux required to perform the parallel tridiagonal matrix-vector
* product
         if (me.ne.(nproc-1)) then
*  receive the last entry of the vector to be involved in the mat-vec.
* this entry is computed by the next processor
            call MPI_RECV(aux(istart+nloc),1,MPI_DOUBLE_COMPLEX,
     &                me+1,3,comm_outer,status,infompi)
         endif
         if (me.ne.0) then
*  receive the first entry of the vector to be involved in the mat-vec.
* this entry is computed by the previous processor
            call MPI_RECV(aux(1),1,MPI_DOUBLE_COMPLEX,me-1,
     &                2,comm_outer,status,infompi)
         endif
*  Compute the local matrix-vector product
         call zgemv('N',nloc,nbCol,ONE,A(1,jstart),lda,aux,1,
     &            ZERO,work(colz),1)
*
         goto 10
*
       else if (revcom.eq.precondRight) then
* perform the right preconditioning for the zfgmres iteration
*
* Check if there is enough space left in the workspace to perfrom
* few steps of zgmres as right preconditioner
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
*
         if (m2.gt.0) then
* Perform at most two restarts in the inner zgmres
           icntl2(7) = 2*m2
* copy colx in the workspace (right hand side location) of the inner
* gmres iteration 
           call zcopy(nloc,work(colx),1,work(irc(6)+nloc),1)
 20        call drive_zgmres(n,nloc,m2,irc(7),
     &             work(irc(6)),irc2,icntl2,cntl2,info2,rinfo2) 
           revcom2 = irc2(1)
           colx2   = irc2(2) + irc(6) -1
           coly2   = irc2(3) + irc(6) -1
           colz2   = irc2(4) + irc(6) -1
           nbscal2 = irc2(5)
           if (revcom2.eq.matvec) then
* Perform the matrix vector product for the inner zgmres iteration
             call zcopy(nloc,work(colx2),1,aux(istart),1)
* Send the entry of aux required to perform the parallel tridiagonal matrix-vector
* product
             if (me.ne.(nproc-1)) then
*  send the last entry of y local vector to the next processor that
*   needs this entry to perform its part of the matrix-vector product
                call MPI_SEND(aux(istart+nloc-1),1,MPI_DOUBLE_COMPLEX,
     &                me+1,2,comm_inner,infompi)
             endif
             if (me.ne.0) then
*  send the first entry of y local vector to the previous processor that
*   needs this entry to perform its part of the matrix-vector product
                call MPI_SEND(aux(istart),1,MPI_DOUBLE_COMPLEX,me-1,
     &                3,comm_inner,infompi)
             endif
* Receive the entry of aux required to perform the parallel tridiagonal matrix-vector
* product
             if (me.ne.(nproc-1)) then
*  receive the last entry of the vector to be involved in the mat-vec.
* this entry is computed by the next processor
                call MPI_RECV(aux(istart+nloc),1,MPI_DOUBLE_COMPLEX,
     &                me+1,3,comm_inner,status,infompi)
             endif
             if (me.ne.0) then
*  receive the first entry of the vector to be involved in the mat-vec.
* this entry is computed by the previous processor
                call MPI_RECV(aux(1),1,MPI_DOUBLE_COMPLEX,me-1,
     &                2,comm_inner,status,infompi)
             endif
*  Compute the local matrix-vector product
             call zgemv('N',nloc,nbCol,ONE,A(1,jstart),lda,aux,1,
     &            ZERO,work(colz2),1)
             goto 20
           else if (revcom2.eq.precondRight) then
* perform the preconditioning for the inner zgmres iteration
             do i =0,nloc-1
               work(colz2+i) = work(colx2+i)/a(i+1,i+1)  
             enddo
             goto 20
           else if (revcom2.eq.precondleft) then
* perform the preconditioning for the inner zgmres iteration
             do i =0,nloc-1
               work(colz2+i) = work(colx2+i)/a(i+1,i+1)  
             enddo
             goto 20
           else if (revcom2.eq.dotProd) then
*      perform the scalar product for the inner zgmres iteration
*      work(colz) <-- work(colx) work(coly)
             call zgemv('C',nloc,nbscal2,ONE,work(colx2),nloc,
     &               work(coly2),1,ZERO,aux,1)
             call MPI_ALLREDUCE(aux,work(colz2),nbscal2,
     &          MPI_DOUBLE_COMPLEX, MPI_SUM,comm_inner,infompi)
             goto 20
           endif
           call zcopy(nloc,work(irc(6)),1,work(colz),1)
           goto 10
         else
* (m2.le.0)
           print *, ' Not enough space for inner gmres '
           call zcopy(nloc,work(colx),1,work(colz),1)
           goto 10
         endif
       else if (revcom.eq.dotProd) then
*      perform the scalar product for the zfgmres iteration
*      work(colz) <-- work(colx) work(coly)
*
         call zgemv('C',nloc,nbscal,ONE,work(colx),nloc,
     &               work(coly),1,ZERO,aux,1)
         call MPI_ALLREDUCE(aux,work(colz),nbscal,
     &          MPI_DOUBLE_COMPLEX, MPI_SUM,comm_outer,infompi)
         goto 10
       endif
*
*******************************
* dump the solution on a file
*******************************
*
      if (me.eq.0) then
        nout = 11
        open(nout,FILE='sol_zTest_Par',STATUS='unknown')
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
        do j=1,min(10,nloc)
          write(nout,*) j,work(j)
        enddo
        write(nout,*) '   '
      endif
*
100   continue
      call MPI_FINALIZE(infompi)
*
      stop
      end

c
c***********************************************************************
c
c                  EXHAUSTIVE DRIVER FOR PREQN
c
c                    * * * * * * * * * * * *
c
c     This driver exercises many of the options in the PREQN package.
c     Two tests are performed: i) Test I calls PREQN with 3 sets of 
c     valid parameters; ii) Test II calls PREQN with 3 sets of invalid
c     options. 
c
c     ------
c     Test I: We solve a sequence of 6 linear systems
c     ------
c                 A_i x = b_i,   i=1,...,6,          (1)
c
c     building a new preconditioner for each problem. This illustrates
c     the typical usage of the package.
c
c     The right hand side vectors $b_i$ are defined as in driver.f, 
c     and the 170x170 matrices A_i are stored, in dense format, in the 
c     files:
c
c         names(1) = 'data2.1'
c         names(2) = 'data2.2'
c         names(3) = 'data2.3'
c         names(4) = 'data2.4'
c         names(5) = 'data2.5'
c         names(6) = 'data2.6'
c
c     Other files required by Test I are:
c
c        pcg.f, extras.f, preqn.f, blas.f, data.dat
c
c     See driver1.f for an explanation of the first 4 files; the file
c     data.dat contains data for this test. 
c
c     -------
c     Test II: We try to solve the system 
c     -------
c                  A_6 x = b_6,                     (2)
c     
c     by calling PREQN with 3 sets of invalid parameters. Test II 
c     essentially exercises the routines that check errors. This test
c     requires the same files needed in Test I.
c
c
c***********************************************************************

      integer  nax, max
      parameter(nax = 170, max = 110)
      integer  n, nz, jptra(nax+1), indra(nax*10), iwp(2*max + 3), i,
     *         j, m, iop, liwp, lwp, info, maxit, iout, l, mssg, iopval
      double   precision adf(nax*nax), b(nax), x(nax), a(nax*10), 
     *         adiag(nax), 
     *         w(6*nax), wp( (max+1)*(4*nax+2) + 1 ), tolcg,
     *         anorm, dnorma, x0
      logical  build
      character names(6)*12, nfile*12
c
c------------
c     Test I: 
c------------
c     Open input file
c
      open(  8, file = 'data.dat',     status = 'unknown' )
c
c     Set the list of names
c
      names(1) = 'data2.1'
      names(2) = 'data2.2'
      names(3) = 'data2.3'
      names(4) = 'data2.4'
      names(5) = 'data2.5'
      names(6) = 'data2.6'
c
c     Read parameters for the PCG and PREQN routines
c      
      read(5,*) n, m, iop, maxit, tolcg, iout, mssg, x0
c
c------------------
c     Set the parameters for pcg.f and the PREQN package
c
c     l      number of linear systems to be solved
c     m      maximum number of pairs in the preconditioner
c     iop    specifies the saving strategy
c     lwp    length of real workspace for the preconditioner
c     liwp   length of integer workspace for the preconditioner
c-------------------
c
      l     = 6
      lwp   = (m+1)*(4*n+2) + 1
      liwp  = 2*m + 30
      build = .false.
c
c      Perform the following test for 3 different settings
c      of parameters
c
      write(6,930)
c
      do 300 iopval=1, 3
c    
c     Compute the first right hand side and define the initial point.
c     The initial point has the same value (x0) for all components
c
         call rhsg ( n, b )
c
         do 10 j=1, n
            x(j) = x0
 10      continue
c
c     Call the PCG routine sequentially.
c
         do 100 i=1, l
c
            nfile = names(i)
            open( 9, file = nfile, status = 'unknown' )         
            read(9,*) n
            read(9,*)( adf(j), j=1, n*n )
            close(9)
c
c     Set the non-zero structure of A_i. Compute ||A_i||. 
c
            call sparse ( adf, a, n, nz, adiag, jptra, indra  )
            anorm = dnorma ( n, a, adiag, jptra, indra, w )
c
c
c     Build the preconditioner after solving the first subproblem
c
            if ( i.ge.2 ) then
               build = .true.
            end if
c
            call pcg    ( n, nz,  a, adiag, jptra, indra, x, b, anorm,
     *                    w, 
     *                    m, iop, i, wp, lwp, iwp, liwp, 
     *                    build, info,
     *                    maxit, tolcg, iout, mssg )
            if ( info.ne.0 ) go to 200
c
c     Perturb the rhs vector, reset x to the initial point x0
c
            call perrhs ( n, b )
c
            do 90 j=1, n
               x(j) = x0
 90         continue
c
 100     continue    
c
 200     continue
c 
c     Sometimes it is appropriate to use the solution of system 
c     A_i x = b_1 as an initial approximation for the new system  
c     A_{i+1} x = b_{i+1}. An easy way to do this is to comment out
c     or remove the statements
c
c            do 90 j=1, n
c               x(j) = x0
c     90     continue
c 
c     in which x is reset to x0. 
c
c     Repeat the solution using now IOP=2
c     
         if ( iopval.eq.1 ) then
            iop = 2
            m   = 7 
            write (6,940)
         end if 
c
c     Define IOP=1 but set M=3 (should be an even number)
c
         if ( iopval.eq.2 ) then
            iop = 1
            m   = 3 
            write (6,950)
         end if       
c
 300  continue 
c
c-------------
c     Test II: 
c-------------
c
c     Call pcg with 3 sets of invalid parameters
c
      m     = 30
      iop   = 1
      i     = 1
      build = .false.
c
      do 400 iopval=1, 3
c
c     Set the length of array wp to a number less than (m+1)*(4*n+2)+1
c
         if ( iopval.eq.1) then
            lwp  = (m+1)*(4*n+2) 
            liwp = 2*m + 30 
            write (6,980)
         end if
c
c     Set the length of array iwp to a number less than 2*m + 30
c
         if ( iopval.eq.2) then
            lwp  = (m+1)*(4*n+2) + 1
            liwp = 2*m + 29  
            write (6,990)
         end if
c
c     Set IOP=0
c
         if ( iopval.eq.3) then
            iop  = 0
            liwp = 2*m + 30
            write (6,960)
         end if 
c
         call pcg    ( n, nz,  a, adiag, jptra, indra, x, b, anorm,
     *                 w, 
     *                 m, iop, i, wp, lwp, iwp, liwp, 
     *                 build, info,
     *                 maxit, tolcg, iout, mssg )
c
 400  continue
c
 930  format(///,2x,'<<<<<<     Set IOP = 1  and M = 30       >>>>>> ') 
 940  format(///,2x,'<<<<<<     Set IOP = 2  and M =  7       >>>>>> ')
 950  format(///,2x,'<<<<<<     Set IOP = 1  and M =  3       >>>>>> ')  
 960  format(///,2x,'<<<<<<           Set IOP = 0             >>>>>> ')       
 980  format(///,2x,'<<<<<<    Set LWP  < (m+1)*(4*n+2)+1     >>>>>> ')
 990  format(///,2x,'<<<<<<       Set LIWP < 2*m + 30         >>>>>> ')
c
      stop
      end    








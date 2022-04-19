      program driver
c     -----------------------------------------------------
c     Fortran driver for subroutine gqaps (sparse GRASP for
c     the QAP).
c
c     This file also contains 2 subroutines:
c     
c          readp  - read QAP instance     
c          outsol - write output report
c
c     -----------------------------------------------------
c     Parameters, variables and arrays needed for the         
c     subroutine gqaps.
c     -----------------------------------------------------  
      integer   nmax,nsq,in,out
      parameter (nmax = 256,nsq = nmax*nmax,in = 5,out = 6)
      integer   f(nsq),d(nsq),a(nmax),b(nmax)
      integer   srtf(nsq),srtif(nsq),srtd(nsq)
      integer   srtid(nsq),srtc(nsq),srtic(nsq)
      integer   indexd(nsq),indexf(nsq),cost(nsq)
      integer   fdind(nsq),perm(nmax)
      integer   cp(nsq),optb(nmax)
      integer   iter,maxitr,bestv,n,look4,seed,n2
      real      alpha,beta,sparsd,sparsf
c     ------------------------------------------------------------------
c     Variables used by driver (not needed by gqaps).
c     ------------------------------------------------------------------
      integer   iseed0
c     ------------------------------------------------------------------
c     Read problem data and gqaps parameters.
c     Write out summary of input.
c
c     Read problem dimension (n).
c     ------------------------------------------------------------------
      read (in,*) n
      if (n .gt. nmax) then
           write(out,10)
10         format(' error: n > nmax  (increase nmax in driver.f)')
           stop
      else
           n2 = n*n
c          -------------------------------------------------------------
c          Read matrices D and F and determine the sparsity of each.
c          -------------------------------------------------------------
           call readp(in,n,n2,sparsf,sparsd,f,d)
      endif
c     ------------------------------------------------------------------
c     Run sparse GRASP for QAP.
c     ------------------------------------------------------------------
c
c     Variables needed for input:
c       
c         n      : dimension of QAP
c                  integer
c         n2     : n * n
c         maxitr : maximum number of GRASP iterations
c                  integer
c         alpha  : GRASP construction phase parameter alpha (0.5)
c                  real  
c         beta   : GRASP construction phase parameter beta (0.1)
c                  real  
c         look4  : GRASP returns permutation if solution with cost
c                  less than or equal to look4 is found (look4 = -1
c                  causes GRASP to take maxitr iterations
c                  integer
c         seed   : seed for random number generator (270001)
c                  integer
c
c     Integer arrays needed for input:
c
c         f      : flow matrix (represented as row-by-row array of
c                  dimension NMAX*NMAX)
c         d      : distance matrix (represented as row-by-row array of
c                  dimension NMAX*NMAX)
c         
c     Integer arrays needed for work:
c   
c         a      : dimension NMAX
c         b      : dimension NMAX
c         optb   : dimension NMAX
c         srtf   : dimension NMAX*NMAX
c         srtif  : dimension NMAX*NMAX
c         srtd   : dimension NMAX*NMAX
c         srtid  : dimension NMAX*NMAX
c         srtc   : dimension NMAX*NMAX
c         srtic  : dimension NMAX*NMAX
c         indexd : dimension NMAX*NMAX
c         indexf : dimension NMAX*NMAX
c         cost   : dimension NMAX*NMAX
c         fdind  : dimension NMAX*NMAX
c         cp     : dimension NMAX*NMAX
c         
c     Integer array needed for output:
c   
c         perm   : permutation vector of best solution found 
c                  dimension NMAX
c
c     Integer variables needed for output:
c   
c         bestv  : cost of best assignment found
c         iter   : number of GRASP iterations taken  
c  
c     ------------------------------------------------------------------
c     Set GRASP parameters.
c     ------------------------------------------------------------------
      maxitr = 1000
      alpha  = 0.25
      beta   = 0.5
      look4  = -1
      seed   = 270001
      iseed0 = seed

c     ------------------------------------------------------------------
c     Find approximate solution to QAP.
c     ------------------------------------------------------------------
      call gqaps(n,n2,maxitr,alpha,beta,look4,seed,f,d,a,b,
     &           optb,srtf,srtif,srtd,srtid,srtc,srtic,indexd,
     &           indexf,cost,fdind,perm,bestv,iter,cp)
c     ------------------------------------------------------------------
c     Write output summary.
c     ------------------------------------------------------------------
      call outsol(out,n,iter,iseed0,bestv,maxitr,look4,alpha,beta,perm,
     &            sparsf,sparsd)
c     ------------------------------------------------------------------
      stop

c     end of program driver
      end



      subroutine readp(in,n,n2,sparsf,sparsd,f,d)
c     ------------------------------------------------------------------
c     readp: Read arrays d and f of dimension n2 and compute the         
c            sparsity of each. Sparsity is the ratio of the nondiagonal
c            nonzero elements to the zero elements. 
c     ------------------------------------------------------------------    
c     Passed input scalars:
c
c               in      - input device
c               n       - QAP dimension
c               n2      - n * n
c
c     ------------------------------------------------------------------
      integer in, n, n2
c     ------------------------------------------------------------------
c     Passed output scalars:
c
c               sparsd  - sparsity of array d
c               sparsf  - sparsity of array f
c
c     ------------------------------------------------------------------
      real    sparsd,sparsf
c     ------------------------------------------------------------------
c     Passed output arrays:
c
c               f       - flow matrix stored as 1-dimensional array, 
c                         row by row.
c               d       - distance matrix stored as 1-dimensinal array,
c                         row by row.
c
c     ------------------------------------------------------------------
      integer f(n2),d(n2)
c     ------------------------------------------------------------------      
c     Local scalars:
c
c               i       - index
c               j       - index
c               countd  - number of nondiagonal 0 entries in array d
c               countf  - number of nondiagonal 0 entries in array f
c
c     ------------------------------------------------------------------
      integer i,j,countd,countf

c     -------------------------------------------------------------------
c     Read matrices D and F.
c     -------------------------------------------------------------------
      do 10 i = 1,n
           read (in,*) (d((i - 1)*n + j), j = 1,n)
10    continue
      do 20 i = 1,n
           read (in,*) (f((i-1)*n + j), j = 1,n)
20    continue
c     ------------------------------------------------------------------
c     Compute the sparsity of matrices D and F.
c     ------------------------------------------------------------------
      countd = 0
      countf = 0
      do 30 i = 1,n2
         if (d(i) .eq. 0) then
            countd = countd + 1
         endif
         if (f(i) .eq. 0) then
            countf = countf + 1
         endif
30    continue
      countd = countd - n
      countf = countf - n
      sparsd = real(countd)/real(n2)
      sparsf = real(countf)/real(n2)
c     ------------------------------------------------------------------
      return

c     end of readp
      end
      
      
      
      subroutine outsol(out,n,iter,iseed0,opt,maxitr,look4,
     +                  alpha,beta,opta,sparsf,sparsd)
c     ------------------------------------------------------------------
c     Output best solution found and run statistics.
c     ------------------------------------------------------------------
c     Passed input scalars:
c
c          out    - output unit
c          n      - QAP dimension
c          iter   - number of GRASP iterations taken
c          iseed0 - initial random number generator seed
c          opt    - value of best permutation found
c          maxitr - maximum number GRASP iterations
c          look4  - value of permutation sought
c          alpha  - GRASP construction phase parameter
c          beta   - GRASP construction phase parameter
c          sparsf - sparsity of matrix f
c          sparsd - sparsity of matrix d
c          
c     ------------------------------------------------------------------
      integer out,n,iter,iseed0,opt,maxitr,look4
      real    alpha,beta,sparsf,sparsd
c     ------------------------------------------------------------------
c     Passed input array:
c
c          opta   - array of best permutation so far 
c
c     ------------------------------------------------------------------
      integer opta(n)
c     ------------------------------------------------------------------
c     Local scalar:
c
c          i      - do loop index
c
c     ------------------------------------------------------------------
      integer i
c     ------------------------------------------------------------------
      write(out,10)
10    format(
     +'----------------------------------------------------------')
      write(out,20)
20    format(
     +'G R A S P for Q A P---------------------------------------')
      write(out,30)
30    format('                                       ')

      write(out,40)
40    format(
     +'input-----------------------------------------------------')
      write(out,50) n
50    format(' dimension of qap                   : ',i20)
      write(out,60) alpha
60    format(' construction phase parameter alpha : ',f20.2)
      write(out,70) beta
70    format(' construction phase parameter beta  : ',f20.2)
      write(out,80) sparsf
80    format(' sparsity of matrix f               : ',f20.2)
      write(out,90) sparsd
90    format(' sparsity of matrix d               : ',f20.2)
      write(out,100) maxitr
100   format(' maximum number of grasp iterations : ',i20)
      write(out,110) look4
110   format(' look4                              : ',i20)
      write(out,120) iseed0
120   format(' initial seed                       : ',i20)
      write(out,30)

      write(out,130)
130   format(
     +'output----------------------------------------------------')
      write(out,140) iter
140   format(' grasp iterations                   : ',i20)
      write(out,150) opt
150   format(' cost of best permutation found     : ',i20)
      write(out,160) (opta(i),i = 1,n)
160   format(' best permutation found             : ',
     &       5i4 / (36x,': ',5i4))
      write(out,30)
      write(out,10)
c     ------------------------------------------------------------------
      return

c     end of outsol
      end


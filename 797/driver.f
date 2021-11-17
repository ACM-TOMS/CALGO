      program driver
c     ------------------------------------------------------------------
c
c     driver:  Driver program for gmpsg.f, Fortran subroutines for
c              a GRASP for graph planarization.
c
c              Inputs a graph with readp(), calls gmpsg() to find large
c              planar subgraph, and prints report with outsol().
c
c     authors: Mauricio G.C. Resende ( mgcr@research.att.com )
c              Celso C. Ribeiro      ( celso@inf.puc-rio.br  )
c
c     Reference: M.G.C. Resende & C.C. Ribeiro, "A GRASP for graph
c                planarization," Networks, vol.29, pp. 173- 189, 1997.
c
c     ------------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Input/ouput unit numbers
c
c          in     - input device unit
c          iout   - output device unit
c
c     ------------------------------------------------------------------
      integer in,             iout

      parameter ( in=5,iout=6 )

c     ------------------------------------------------------------------
c
c     Maximum dimensions arrays.
c
c          nmax   - max number of nodes
c          mmax   - max number of arcs
c          mmax2  - mmax * mmax
c          mmaxp1 - mmax + 1
c          mmaxt2 - mmax * 2
c          nmaxp1 - nmax + 1
c
c     ------------------------------------------------------------------
c
      integer nmax,           mmax,           mmax2,
     +        mmaxp1,         mmaxt2,         nmaxp1

      parameter ( nmax=100,mmax=1000,mmax2=mmax*mmax,
     +            mmaxt2=2*mmax,mmaxp1=mmax+1,nmaxp1=nmax+1 )

c     ------------------------------------------------------------------
c
c     Integer arrays:
c
c         adj(mmax2)     - adjacency matrix
c         arc1(mmax)     - node 1 of arc in input graph
c         arc1r(mmax)    - node 1 of arc in reduced graph
c         arc2(mmax)     - node 2 of arc in input graph
c         arc2r(mmax)    - node 2 of arc in reduced graph
c         bigw(mmax)     - set W in exact stable set algorithm
c         blue(mmax)     - set of blue arcs
c         currnt(mmax)   - current node in max stable algorithm
c         deg(nmax)      - node degrees
c         deg0(nmax)     - node degrees of input graph
c         extrem(mmaxt2) - node in adjcency list
c         iniadj(mmaxp1) - pointer to first node in adjacency list
c         iniux(mmaxp1)  - pointer to first interval in U(x) list
c         invpi(nmax)    - inverse of node permutation
c         invred(mmax)   - maps nodes in compressed graph to input
c         kset(mmax)     - set of stable nodes
c         lstblu(mmax)   - list of blue arcs that can be colored red
c         lstpal(mmax)   - list of pale arcs that can be colored blue
c         mark(mmax)     - mark in max stable set algorithm
c         nextux(mmax2)  - next element of U(x) list
c         numarc(mmaxt2) - arc number
c         nxtadj(mmax2)  - next node in adjacency list
c         optblu(mmax)   - set of best blue arcs
c         optprm(nmax)   - node permutation of best solution
c         optred(mmax)   - set of best red arcs
c         pale(mmax)     - set of pale arcs
c         pi(nmax)       - node permutation
c         point(m)       - POINTER in MAXSTABLE
c         ptrnod(nmaxp1) - pointer to first node in adjacency list
c         rcl(nmax)      - GRASP restricted candidade list
c         red(mmax)      - set of red arcs
c         sset(mmax2)    - set S(*)
c         stable(mmax)   - set of stable nodes (arcs in input graph)
c         stack(mmax2)   - stack to simulate recursive logic
c         stack2(mmax2)  - stack to simulate recursive logic
c         ux(mmax2)      - set UX for stable set algorithm
c         vrtex(m)       - Nodes of set U(x)
c         w(m)           - weight
c
c     ------------------------------------------------------------------
      integer adj(mmax2),     arc1(mmax),     arc1r(mmax),
     +        arc2(mmax),     arc2r(mmax),    bigw(mmax),
     +        blue(mmax),     currnt(mmax),   deg(nmax),
     +        deg0(nmax),     extrem(mmaxt2), iniadj(mmaxp1),
     +        iniux(mmaxp1),  invpi(nmax),    invred(mmax),
     +        kset(mmax),     lstblu(mmax),   lstpal(mmax),
     +        mark(mmax),     nextux(mmax2),  numarc(mmaxt2),
     +        nxtadj(mmax2),  optblu(mmax),   optprm(nmax),
     +        optred(mmax),   pale(mmax),     pi(nmax),
     +        point(mmax),    ptrnod(nmaxp1), rcl(nmax),
     +        red(mmax),      sset(mmax2),    stable(mmax),
     +        stack(mmax2),   stack2(mmax2),  ux(mmax2),
     +        vrtex(mmax),    w(mmax)

c     ------------------------------------------------------------------
c
c     Integer parameters:
c
c         errcnd         - error condition
c         gs             - size of planar subgraph
c         iter           - GRASP iteration counter
c         look4          - look for planar graph with look4 arcs
c         m              - number of arcs in input graph
c         maxgs          - size of largest planar subgraph found
c         maxitr         - maximum number of GRASP iterations
c         mred           - number of arcs in reduced graph
c         n              - number of nodes in input graph
c         nblue          - number of blue arcs
c         nopblu         - number of blue arcs in best planar subgraph
c         nopred         - number of red arcs in best planar subgraph
c         npale          - number of pale arcs
c         nred           - number of red arcs
c         nstbl          - number of nodes in stable set
c         optitr         - GRASP iteration best solution was found
c         prttyp         - type of print out (o,1,2)
c         seed           - seed for pseudo random number generator
c
c     ------------------------------------------------------------------
      integer errcnd,         gs,             iter,
     +        look4,          m,              maxgs,
     +        maxitr,         mred,           n,
     +        nblue,          nopblu,         nopred,
     +        npale,          nred,           nstbl,
     +        optitr,         prttyp,         seed

c     ------------------------------------------------------------------
c
c     Real parameter:
c
c         alpha          - RCL parameter
c
c     ------------------------------------------------------------------
      real    alpha

c     ------------------------------------------------------------------
c
c     Double precision parameter:
c
c         avggs          - sum of solutions found
c
c     ------------------------------------------------------------------
      double precision        avggs

c     ------------------------------------------------------------------
c     Input problem:
c     
c     n         - number of nodes
c     m         - number of arcs
c     arc1      - node 1 of arc
c     arc2      - node r of arc
c
c     ------------------------------------------------------------------
      call readp( arc1,   arc2,   errcnd, in,     m,      mmax,
     +            n,      nmax )

c     ------------------------------------------------------------------
c     If an error was encountered in input, print error type and
c     terminate execution.
c     ------------------------------------------------------------------
      if ( errcnd .gt. 0) then
           call errmsg( errcnd, iout )
           stop
      endif
      
c     ------------------------------------------------------------------
c     Set algorithm parameters:
c
c     alpha  - restricted candidate list parameter 
c              ( 0 .le. alpha .le. 1 ) real
c              ( alpha = 0: greedy construction )
c              ( alpha = 1: random construction )
c
c     look4  - look for a planar subgraph with at least look4 arcs
c              ( 1 .le. look4 .le. m ) integer
c
c     maxitr - maximum number of GRASP iterations
c              ( maxitr .ge. 1 ) integer
c
c     prttyp - type of output report 
c              ( prttyp = 0, 1, or 2 ) integer
c              ( prttyp = 0: silent run with no printing from gmpsg    )
c              ( prttyp = 1: print only solution improvements in gmpsg )
c              ( prttyp = 2: print size of solutions found in gmpsg    )
c
c     seed   - random number generator seed
c              ( 1 .le. seed .le. 2147483647 ) integer
c
c     ------------------------------------------------------------------
      alpha  = 0.1
      look4  = m
      maxitr = 2048
      prttyp = 1
      seed   = 270001

c     ------------------------------------------------------------------
c     Call main optimization module (gmpsg) to produce a large planar
c     subgraph of the input graph.
c     ------------------------------------------------------------------
      call gmpsg( adj,    alpha,  arc1,   arc1r,  arc2,   arc2r,
     +            avggs,  bigw,   blue,   currnt, deg,    deg0,   
     +            errcnd, extrem, gs,     iniadj, iniux,  invpi,
     +            invred, iout,   iter,   look4,  kset,   lstblu,
     +            lstpal, m,      mark,   maxgs,  maxitr, mmax,
     +            mred,   n,      nblue,  nextux, nmax,   nopblu,
     +            nopred, npale,  nred,   nstbl,  numarc, nxtadj,
     +            optblu, optitr, optprm, optred, pale,   pi,
     +            point,  prttyp, ptrnod, rcl,    red,    seed,
     +            sset,   stable, stack,  stack2, ux,     vrtex,
     +            w )

c     ------------------------------------------------------------------
c     Print out best solution found or error message.
c     ------------------------------------------------------------------
      call errmsg( errcnd, iout )
      if ( errcnd .eq. 0) then
           call outsol( avggs,  iout,   iter,   m,      maxgs,  n,
     +                  nopblu, nopred, optblu, optitr, optprm, optred )
      endif

c     ------------------------------------------------------------------
c     End of driver.
c     ------------------------------------------------------------------
      stop
      end




      subroutine readp( arc1,   arc2,   errcnd, in,     m,      mmax,
     +                  n,      nmax )
c     ------------------------------------------------------------------
c     readp: Input problem instance.
c     ------------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Input integer parameters:
c
c         in             - input device number
c         m              - number of arcs
c         mmax           - max value of dimension m
c         n              - number of nodes
c         nmax           - max value of dimension n
c
c     ------------------------------------------------------------------
      integer in,             m,              mmax,
     +        n,              nmax

c     ------------------------------------------------------------------
c
c     Integer working parameter:
c
c         edge           - arc
c         na             - arc counter
c
c     ------------------------------------------------------------------
      integer edge,           na

c     ------------------------------------------------------------------
c
c     Integer output parameter:
c
c         errcnd         - error condition
c
c     ------------------------------------------------------------------
      integer errcnd

c     ------------------------------------------------------------------
c
c     Integer output arrays:
c
c         arc1(m)        - node 1 of arc in input graph
c         arc2(m)        - node 2 of arc in input graph
c
c     ------------------------------------------------------------------
      integer arc1(mmax),     arc2(mmax)

c     ------------------------------------------------------------------
c     Input graph dimensions.
c     ------------------------------------------------------------------
      read (in,*) n,m
c     ------------------------------------------------------------------
c     Check validity of dimensions.
c     ------------------------------------------------------------------
      errcnd = 0
      if (n .gt. nmax) then
           errcnd=1
           return
      endif
      if (m .gt. mmax) then
           errcnd=2
           return
      endif

c     ------------------------------------------------------------------
c     Read arc (i,j) 
c     ------------------------------------------------------------------
      na=0
      do 10 edge=1,m
           read(in,*,end=20) arc1(edge), arc2(edge)
           na=na+1
10    continue
20    if (na .lt. m) then
           errcnd=9
           return
      endif
c     ------------------------------------------------------------------
c     End of readp.
c     ------------------------------------------------------------------
      end




      subroutine outsol( avggs,  iout,   iter,   m,      maxgs,  n,
     +                   nopblu, nopred, optblu, optitr, optprm, optred 
     +                 )
c     ------------------------------------------------------------------
c     outsol: Print report with best solution found.
c     ------------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Integer input parameters:
c
c         iout           - output device number
c         iter           - GRASP iteration counter
c         m              - number of arcs in input graph
c         maxgs          - size of largest planar subgraph found
c         n              - number of nodes in input graph
c         nopblu         - number of blue arcs in best planar subgraph
c         nopred         - number of red arcs in best planar subgraph
c         optitr         - GRASP iteration best solution was found
c
c     ------------------------------------------------------------------
      integer iout,           iter,           m,
     +        maxgs,          n,              nopblu,
     +        nopred,         optitr

c     ------------------------------------------------------------------
c
c     Double precision input parameter:
c
c         avggs          - sum of solutions found
c
c     ------------------------------------------------------------------
      double precision        avggs

c     ------------------------------------------------------------------
c
c     Integer work parameter:
c
c         i              - do loop parameter
c
c     ------------------------------------------------------------------
      integer i

c     ------------------------------------------------------------------
c
c     Integer input arrays:
c
c         optblu(m)      - set of best blue arcs
c         optprm(n)      - node permutation of best solution
c         optred(m)      - set of best red arcs
c
c     ------------------------------------------------------------------
      integer optblu(m),      optred(m),      optprm(n)

c     ------------------------------------------------------------------

      write(iout,10)
10    format(' ')
      write(iout,20)
20    format(
     +' GRASP solution------------------------------------------------',
     +'-----------')
      write(iout,10)

      write(iout,30) optitr
30    format(' largest planar subgraph found in iteration: ',i7)
      write(iout,10)

      write(iout,40) maxgs
40    format(' size of largest planar subgraph: ',i5)
      write(iout,50) avggs/iter
50    format(' avg  size of    planar subgraph: ',f7.1)
      write(iout,60) nopblu
60    format(' blue arcs in    planar subgraph: ',i5)
      write(iout,70) nopred
70    format(' red  arcs in    planar subgraph: ',i5)
      write(iout,10)

      write(iout,80) (optprm(i),i=1,n)
80    format(' vrtx permut:', 5i10)
      write(iout,10)

      write(iout,90) (optblu(i),i=1,nopblu)
90    format(' blue  arcs: ', 5i10)
      write(iout,10)

      write(iout,100) (optred(i),i=1,nopred)
100   format(' red   arcs: ', 5i10)
      write(iout,10)

      write(iout,110)
110   format(
     +' --------------------------------------------------------------',
     +'-----------')
      write(iout,10)

c     ------------------------------------------------------------------
c     End of outsol. 
c     ------------------------------------------------------------------
      return
      end

      subroutine gmpsg( adj,    alpha,  arc1,   arc1r,  arc2,   arc2r,
     +                  avggs,  bigw,   blue,   currnt, deg,    deg0,   
     +                  errcnd, extrem, gs,     iniadj, iniux,  invpi,  
     +                  invred, iout,   iter,   look4,  kset,   lstblu,
     +                  lstpal, m,      mark,   maxgs,  maxitr, mmax,
     +                  mred,   n,      nblue,  nextux, nmax,   nopblu,
     +                  nopred, npale,  nred,   nstbl,  numarc, nxtadj,
     +                  optblu, optitr, optprm, optred, pale,   pi,
     +                  point,  prttyp, ptrnod, rcl,    red,    seed,
     +                  sset,   stable, stack,  stack2, ux,     vrtex,
     +                  w )
c     ------------------------------------------------------------------
c     gmpsg.f: GRASP for graph planarization.
c
c     authors: Mauricio G.C. Resende ( mgcr@research.att.com )
c              Celso C. Ribeiro      ( celso@inf.puc-rio.br  )
c
c     Reference: M.G.C. Resende & C.C. Ribeiro, "A GRASP for graph
c                planarization," Networks, vol.29, pp. 173- 189, 1997.
c
c     ------------------------------------------------------------------
c
c     Subroutine gmpsg takes as input a graph ( m, n, arc1(*), arc2(*) )
c     and memory allocation parameters ( nmax, mmax ):
c
c           m       - number of arcs
c                     ( 1 .le. m .le. mmax )
c
c           n       - number of nodes
c                     ( 1 .le. n .le. nmax )
c
c           arc1(i) - node 1 of arc i
c                     ( 1 .le. arc1(i) .le. n )
c
c           arc2(i) - node 2 of arc i
c                     ( 1 .le. arc2(i) .le. n )
c
c     and parameters for running the GRASP ( alpha, look4, maxitr,
c     prttyp, seed ):
c
c           alpha   - restricted candidate list parameter 
c                     ( 0 .le. alpha .le. 1 ) real
c                     ( alpha = 0: greedy construction )
c                     ( alpha = 1: random construction )
c                     ( Rule of thumb for setting alpha: alpha=0.5 for
c                       small (n < 30) graphs; alpha=0.1 for larger
c                       graphs. )
c
c           look4   - look for a planar subgraph with at least look4 
c                     arcs
c                     ( 1 .le. look4 .le. m ) integer
c
c           maxitr  - maximum number of GRASP iterations
c                     ( maxitr .ge. 1 ) integer
c                     ( Rule of thumb for setting maxitr: larger
c                       values of maxitr increase the probability
c                       of finding larger planar subraphs but also
c                       increases running time. )
c
c           prttyp  - type of output report 
c                     ( prttyp = 0, 1, or 2 ) integer
c                     ( prttyp = 0: silent run with no printing from 
c                       gmpsg    )
c                     ( prttyp = 1: print only solution improvements 
c                       in gmpsg )
c                     ( prttyp = 2: print size of solutions found in
c                       gmpsg    )
c
c           seed    - random number generator seed
c                     ( 1 .le. seed .le. 2147483647 ) integer
c                     ( Rule of thumb for setting seed: Different seeds
c                       produce different runs (if alpha <  1).  If more
c                       than one processor is available, set a different 
c                       seed for each processor. )
c
c     Input data is checked and error condition (errcnd > 0) is returned
c     if an error in the input data is encountered.
c
c     If no error is present, errcnd = 0 is returned, along with a
c     large planar subgraph, represented by a node permutation
c     (optprm), a list of red arcs (optred) and a list of blue arcs
c     (optblu).  The number of red (blue) arcs is nopred (nopblu). 
c
c     The planar subgraph is drawn by placing the nodes in line in the
c     manner indicated by optprm and drawing blue arcs above the line
c     and red arcs below the line.  Red (blue) arcs do not intersect 
c     with respect to permutation optprm.
c
c     ------------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Integer input parameters:
c
c         iout           - output device number
c         look4          - look for planar subgraph of size look4
c         m              - number of arcs in input graph
c         maxitr         - maximum number of GRASP iterations
c         mmax           - max value of dimension m
c         n              - number of nodes in input graph
c         nmax           - max value of dimension n
c         prttyp         - type of reporting (0=silent,1=improvement,
c                          2=all iterations)
c         seed           - seed for pseudo random number generator
c
c     ------------------------------------------------------------------
      integer iout,           look4,          m,
     +        maxitr,         mmax,           n,
     +        nmax,           prttyp,         seed

c     ------------------------------------------------------------------
c
c     Integer work parameters:
c
c         gs             - size of planar subgraph
c         iter           - GRASP iteration counter
c         m2             - m * m
c         maxgs          - size of largest planar subgraph found
c         maxux          - size of set U(x)
c         mp1            - m + 1
c         mred           - number of arcs in reduced graph
c         mred*2         - mred * mred
c         mt2            - m * 2
c         nblue          - number of blue arcs
c         np1            - n + 1
c         npale          - number of pale arcs
c         nred           - number of red arcs
c
c     ------------------------------------------------------------------
      integer gs,             iter,           m2,
     +        maxgs,          maxux,          mp1,
     +        mred,           mred2,          mt2,
     +        nblue,          np1,            npale,
     +        nred
 
c     ------------------------------------------------------------------
c
c     Integer output parameters:
c
c         errcnd         - error condition
c         nopblu         - number of blue arcs in best planar subgraph
c         nopred         - number of red arcs in best planar subgraph
c         nstbl          - number of nodes in stable set
c         optitr         - GRASP iteration best solution was found
c
c     ------------------------------------------------------------------
      integer errcnd,         nopblu,         nopred,
     +        nstbl,          optitr

c     ------------------------------------------------------------------
c
c     real parameter
c
c         alpha          - RCL parameter
c
c     ------------------------------------------------------------------
      real                    alpha

c     ------------------------------------------------------------------
c
c     double precision parameter
c
c         avggs          - sum of solutions found
c
c     ------------------------------------------------------------------
      double precision        avggs

c     ------------------------------------------------------------------
c
c     Integer input arrays:
c
c         arc1(m)        - node 1 of arc in input graph
c         arc2(m)        - node 2 of arc in input graph
c
c     ------------------------------------------------------------------
      integer arc1(m),        arc2(m)

c     ------------------------------------------------------------------
c
c     Integer work arrays:
c
c         adj(m*m)       - adjacency matrix
c         arc1r(m)       - node 1 of arc in reduced graph
c         arc2r(m)       - node 2 of arc in reduced graph
c         bigw(m)        - set W in exact stable set algorithm
c         blue(m)        - set of blue arcs
c         currnt(m)      - current node in max stable algorithm
c         deg(n)         - node degrees
c         deg0(n)        - node degrees of input graph
c         extrem(m*2)    - node in adjcency list
c         iniadj(m+1)    - pointer to first node in adjacency list
c         iniux(m+1)     - pointer to first interval in U(x) list
c         invpi(n)       - inverse of node permutation
c         invred(m)      - maps nodes in compressed graph to input
c         kset(m)        - set of stable nodes
c         lstblu(m)      - list of blue arcs that can be colored red
c         lstpal(m)      - list of pale arcs that can be colored blue
c         mark(m)        - mark in max stable set algorithm
c         nextux(m*m)    - next element of U(x) list
c         numarc(m*2)    - arc number
c         nxtadj(m*m)    - next node in adjacency list
c         optprm(n)      - set of best pale arcs
c         optred(m)      - set of best red arcs
c         pale(m)        - set of pale arcs
c         pi(n)          - node permutation
c         point(m)       - POINTER in MAXSTABLE
c         ptrnod(n+1)    - pointer to first node in adjacency list
c         rcl(n)         - GRASP restricted candidade list
c         red(m)         - set of red arcs
c         sset(m*m)      - set S(*)
c         stable(m)      - set of stable nodes (arcs in input graph)
c         stack(m*m)     - stack to simulate recursive logic
c         stack2(m*m)    - stack to simulate recursive logic
c         vrtex(m)       - Nodes of set U(x)
c         w(m)           - weight
c         ux(m*m)        - set UX for stable set algorithm
c
c     ------------------------------------------------------------------
      integer adj(m*m),       arc1r(m),       arc2r(m),
     +        bigw(m),        blue(m),        currnt(m),
     +        deg(n),         deg0(n),        extrem(m*2),
     +        iniadj(m+1),    iniux(m+1),     invpi(n),
     +        invred(m),      kset(m),        lstblu(m),
     +        lstpal(m),      mark(m),        nextux(m*m),
     +        numarc(m*2),    nxtadj(m*m),    pale(m),
     +        pi(n),          point(m),       ptrnod(n+1),
     +        rcl(n),         red(m),         sset(m*m),
     +        stable(m),      stack(m*m),     stack2(m*m),
     +        ux(m*m),        vrtex(m),       w(m)

c     ------------------------------------------------------------------
c
c     Integer output arrays:
c
c         optblu(m)      - set of best blue arcs
c         optprm(n)      - node permutation of best solution
c         optred(m)      - set of best red arcs
c
c     ------------------------------------------------------------------
      integer optblu(m),      optprm(n),      optred(m)

c     ------------------------------------------------------------------
c     Initialize error condition.
c     ------------------------------------------------------------------
      errcnd=0

c     ------------------------------------------------------------------
c     Check if memory dimension is sufficient, if input graph
c     is consistent, and if algorithm parameters are allowable.
c     Return control to calling program if error is found.
c     ------------------------------------------------------------------
      call chkinp( alpha,  arc1,   arc2,   errcnd, look4,  m,
     +             maxitr, mmax,   n,      nmax,   prttyp, seed )

      if (errcnd.gt.0) return

c     ------------------------------------------------------------------
c     Initialize array dimensions.
c     ------------------------------------------------------------------
      m2=m*m
      np1=n+1
      mp1=m+1
      mt2=m*2

c     ------------------------------------------------------------------
c     Make graph data structures.
c     ------------------------------------------------------------------
      call mkadj( arc1,   arc2,   deg0,   extrem, m,      mt2,
     +            n,      np1,    numarc, ptrnod )

c     ------------------------------------------------------------------
c     Initialize size of best planar graph.
c
c           maxgs  - largest planar subgraph found
c           avggs  - total sum of arcs in planar subgraphs found.
c
c     ------------------------------------------------------------------
      maxgs=-1
      avggs=0.d0

c     ------------------------------------------------------------------
c     Loop for GRASP iterations.
c     ------------------------------------------------------------------
      do 50 iter=1,maxitr
c          -------------------------------------------------------------
c          GRASP construction phase to find node permutation with low
c          number of crossing arcs.
c          -------------------------------------------------------------
           call phase1( alpha,  deg,    deg0,   extrem, invpi,  mt2,
     +                  n,      np1,    pi,     ptrnod, rcl,    seed )

c          -------------------------------------------------------------
c          GRASP l)ocal search phase to find node permutation pi with low
c          number of crossing arcs.
c          -------------------------------------------------------------
           call phase2( arc1,   arc2,   invpi,  m,      mt2,    n, 
     +                  np1,    numarc, pi,     ptrnod )

c          -------------------------------------------------------------
c          Build graph of intersecting arcs (overlap graph).
c
c          Nodes in overlap graph are arcs of original graph. Arcs
c          in overlap graph connect pairs of nodes with corresponding
c          arcs in the original graph that cross in permutation pi.
c          -------------------------------------------------------------
           call mkux( arc1,   arc2,   iniux,  m,      m2,     maxux,
     +                n,      nextux, pi,     ux )

c          -------------------------------------------------------------
c          Find exact maximum stable set of overlap graph.  Stable nodes
c          are  arcs  of the  orginal graph  that do not  cross in 
c          permutation pi.  They can be colored blue and drawn above the
c          line of nodes.
c          -------------------------------------------------------------
           call mxstbl( adj,    arc1,   arc2,   bigw,   currnt, iniadj,
     +                  iniux,  kset,   m,      m2,     mark,   mp1,
     +                  n,      nextux, nstbl,  nxtadj, pi,     point,
     +                  sset,   stable, stack,  stack2, ux,
     +                  vrtex,  w )

c          -------------------------------------------------------------
c          Color stable set blue, remove blue arcs from original graph 
c          and compress original graph for next phase where arcs will be 
c          colored red.
c          -------------------------------------------------------------
           call mkblue ( arc1,   arc1r,  arc2,   arc2r,  blue,   invred,
     +                   m,      mred,   mred2,  nblue,  stable )

c          -------------------------------------------------------------
c          Build graph of intersecting arcs (overlap graph).
c
c          Nodes in overlap graph are arcs of compressed graph. Arcs
c          in overlap graph connect pairs of nodes with corresponding
c          arcs in the compressed graph that cross in permutation pi.
c          -------------------------------------------------------------
           call mkux( arc1,   arc2,   iniux,  mred,   mred2,  maxux,
     +                n,      nextux, pi,     ux )

c          -------------------------------------------------------------
c          Find exact maximum stable set of overlap graph.  Stable nodes
c          are  arcs  of the  compressed graph  that do not  cross in 
c          permutation pi.  They can be colored red and drawn below the
c          line of nodes.
c          -------------------------------------------------------------
           call mxstbl( adj,    arc1r,  arc2r,  bigw,   currnt, iniadj,
     +                  iniux,  kset,   mred,   mred2,  mark,   mp1,
     +                  n,      nextux, nstbl,  nxtadj, pi,     point,
     +                  sset,   stable, stack,  stack2, ux,     vrtex,
     +                  w )
           
c          -------------------------------------------------------------
c          Color stable set red and remaining arcs pale.
c
c          At this point, planar subgraph is made up of blue and red
c          arcs.  Pale arcs are excluded because in permutation pi
c          they cross with one or more red or blue arcs.
c          -------------------------------------------------------------
           call mkred( invred, m,      mred,   npale,  nred,   pale,
     +                 red,    stable )

c          -------------------------------------------------------------
c          Post processing on pale 2-coloring. 
c
c          Look for a blue arc that can be colored red and a pale arc
c          that can be colored blue.  If such arcs are present recolor
c          those arcs and increase size of planar subgraph.
c          -------------------------------------------------------------
           call postp( arc1,   arc2,   avggs,  blue,   gs,     lstblu,
     +                 lstpal, m,      n,      nblue,  npale,  nred,
     +                 pale,   pi,     red )

c          -------------------------------------------------------------
c          If solution found is best so far:
c          -------------------------------------------------------------
           if (gs.gt.maxgs) then
c               --------------------------------------------------------
c               Save it.
c               --------------------------------------------------------
                call savsol( blue,   gs,     invpi,  iter,    m,
     +                       maxgs,  n,      nblue,  nopblu,  nopred,
     +                       nred,   optblu, optitr, optprm,  optred,
     +                       red )

c               --------------------------------------------------------
c               Print solution update if not silent run.
c               --------------------------------------------------------
                if (prttyp .gt. 0) then
                     write(iout,30) iter,maxgs
30                   format(' itr: ',i10,'      size: ',i6,
     +                      ' <- improvement')
                endif
c               --------------------------------------------------------
c               If input graph is proven planar or planar subgraph of
c               size look4 is found, return control to calling program.
c               --------------------------------------------------------
                if ( maxgs.eq.m .or. maxgs.ge.look4 ) return
           else
c               --------------------------------------------------------
c               Print solution update if printing update every 
c               iteration.
c               --------------------------------------------------------
                if (prttyp .gt. 1) then
                     write(iout,40) iter,gs,maxgs
40                   format(' itr: ',i10,' size: ',i6,' best:',i6)
                endif
           endif
50    continue
c     ------------------------------------------------------------------
c     End of gmpsg.
c     ------------------------------------------------------------------
      return
      end




      subroutine chkinp( alpha,  arc1,   arc2,   errcnd, look4,  m,
     +                   maxitr, mmax,   n,      nmax,   prttyp, seed )
c     ------------------------------------------------------------------
c     chkinp: Check input problem instance.
c     ------------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Integer input parameters:
c
c         look4          - look for planar subgraph with look4 arcs
c         m              - number of arc
c         maxitr         - maximum number of iterations
c         mmax           - max value of dimension m
c         n              - number of nodes
c         nmax           - max value of dimension n
c         prttyp         - type of reporting (0=silent,1=improvement,
c                          2=all iterations)
c         seed           - seed for pseudo random number generator
c
c     ------------------------------------------------------------------

      integer look4,          m,              maxitr,
     +        mmax,           n,              nmax,
     +        prttyp,         seed

c     ------------------------------------------------------------------
c
c     Real input parameter:
c
c         alpha          - RCL parameter
c
c     ------------------------------------------------------------------
      real    alpha

c     ------------------------------------------------------------------
c
c     Integer work parameters:
c
c         edge           - arc
c         i              - node
c         j              - node
c
c     ------------------------------------------------------------------
      integer edge,           i,              j

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
c     Integer input arrays:
c
c         arc1(m)        - node 1 of arc in input graph
c         arc2(m)        - node 2 of arc in input graph
c
c     ------------------------------------------------------------------
      integer arc1(mmax),     arc2(mmax)

c     ------------------------------------------------------------------
c     Check n and m.
c     ------------------------------------------------------------------
      if (n.gt.nmax) then
           errcnd=1
           return
      endif
      if (m.gt.mmax) then
           errcnd=2
           return
      endif

c     ------------------------------------------------------------------
c     Scan arcs.
c     ------------------------------------------------------------------
      do 10 edge=1,m
           i=arc1(edge)
           j=arc2(edge)
c          -------------------------------------------------------------
c          Check arc nodes.
c          -------------------------------------------------------------
           if(i.lt.1 .or. i.gt.n. or. j.lt.1 .or. 
     +                    j.gt.n .or. i.eq.j) then
                errcnd=3
                return
           endif
10    continue

c     ------------------------------------------------------------------
c     Check GRASP parameters: alpha, look4, maxitr, seed.
c     ------------------------------------------------------------------
      if (alpha .lt. 0.0e0 .or. alpha .gt. 1.0e0) then
           errcnd=4
           return
      endif

      if (look4 .lt. 1 .or. look4 .gt. m) then
           errcnd=5
           return
      endif

      if (maxitr .lt. 1) then
           errcnd=6
           return
      endif

      if (.not.(prttyp.eq.0 .or. prttyp.eq.1 .or. prttyp.eq.2)) then
           errcnd=7
           return
      endif

      if (seed .lt. 1 .or. seed .gt. 2147483647) then
           errcnd=8
           return
      endif
c     ------------------------------------------------------------------
c     End of chkinp.
c     ------------------------------------------------------------------
      return
      end




      subroutine errmsg( errcnd, iout )
c     ------------------------------------------------------------------
c     errmsg: Prints error message
c     ------------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Integer input parameters:
c
c
c     ------------------------------------------------------------------
c
c         errcnd         - error condition
c         iout           - output device number
c
c     ------------------------------------------------------------------
      integer errcnd,         iout

      if (errcnd .eq. 0) then
           write(iout,10)
10         format(/' Execution terminated with no error. ')
           return
      endif

      if (errcnd .eq. 1) then
           write(iout,20)
20         format(/' Execution aborted.  Error in n (< 1 or > nmax) ')
           return
      endif

      if (errcnd .eq. 2) then
           write(iout,30)
30         format(/' Execution aborted.  Error in m (< 1 or > mmax) ')
           return
      endif

      if (errcnd .eq. 3) then
           write(iout,40)
40         format(/' Execution aborted. ',
     +            ' Error in input node number (< 1 or > n) ')
           return
      endif

      if (errcnd .eq. 4) then
           write(iout,50)
50         format(/' Execution aborted. ',
     +            ' Error in alpha (alpha < 0 or alpha > 1) ')
           return
      endif

      if (errcnd .eq. 5) then
           write(iout,60)
60         format(/' Execution aborted. ',
     +            ' Error in look4 (look4 < 1 or look4 > m) ')
           return
      endif

      if (errcnd .eq. 6) then
           write(iout,70)
70         format(/' Execution aborted.  Error in maxitr (< 1) ')
           return
      endif

      if (errcnd .eq. 7) then
           write(iout,80)
80         format(/' Execution aborted. ',
     +            ' Error in prttyp (not equal to 1, 2, or 3) ')
           return
      endif

      if (errcnd .eq. 8) then
           write(iout,90)
90         format(/' Execution aborted. ',
     +            ' Error in seed (< 1 or > 2147483647) ')
           return
      endif

      if (errcnd .eq. 9) then
           write(iout,100)
100        format(/' Execution aborted. ',
     +            ' Not enough arcs. ')
           return
      endif

c     ------------------------------------------------------------------
c     End of errmsg.
c     ------------------------------------------------------------------
      return
      end




      subroutine mkadj( arc1,   arc2,   deg0,   extrem, m,      mt2,
     +                  n,      np1,    numarc, ptrnod )
c     ------------------------------------------------------------------
c     mkadj: Build graph data structures.
c     ------------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Integer input parameters:
c
c         m              - number of arcs in input graph
c         mt2            - m * 2
c         n              - number of nodes in input graph
c         np1            - n + 1
c
c     ------------------------------------------------------------------
      integer m,              mt2,            n,
     +        np1

c     ------------------------------------------------------------------
c
c     Integer work parameters:
c
c         i              - do loop parameter
c         ilast          - pointer to last element
c         ipos           - pointer to element
c         u              - do loop parameter
c
c     ------------------------------------------------------------------
      integer i,              ilast,          inpos,
     +        u

c     ------------------------------------------------------------------
c
c     Integer input arrays:
c
c         arc1(m)        - node 1 of arc in input graph
c         arc2(m)        - node 2 of arc in input graph
c
c     ------------------------------------------------------------------
      integer arc1(m),        arc2(m)

c     ------------------------------------------------------------------
c
c     Integer output arrays:
c
c         extrem(m*2)    - node in adjcency list
c         deg0(n)        - node degrees of input graph
c         ptrnod(n+1)    - pointer to first node in adjacency list
c         numarc(m*2)    - arc number
c
c     ------------------------------------------------------------------
      integer extrem(mt2),    deg0(n),        ptrnod(np1),
     +        numarc(mt2)

c     ------------------------------------------------------------------
      inpos = 0
      ilast = 0
c     ------------------------------------------------------------------
c     For each node i:
c     ------------------------------------------------------------------
      do 20 i = 1,n
c          -------------------------------------------------------------
c          Scan each arc u:
c          -------------------------------------------------------------
           do 10 u = 1,m
c               --------------------------------------------------------
c               If arc u contains node i, add other node to 
c               list of adjacent nodes of i.
c               --------------------------------------------------------
                if (arc1(u).eq.i) then
                     inpos = inpos + 1
                     numarc(inpos) = u
                     extrem(inpos) = arc2(u)
                elseif (arc2(u).eq.i) then
                     inpos = inpos + 1
                     numarc(inpos) = u
                     extrem(inpos) = arc1(u)
                endif
10         continue
           ptrnod(i) = ilast + 1
           ilast = inpos
20    continue
      ptrnod(n+1) = ilast + 1

c     ------------------------------------------------------------------
c     Compute node degrees.
c     ------------------------------------------------------------------
      do 30 i = 1,n
          deg0(i) = ptrnod(i+1)-ptrnod(i)
30    continue

c     ------------------------------------------------------------------
c     End of mkadj.
c     ------------------------------------------------------------------
      return
      end




      subroutine phase1( alpha,  deg,    deg0,   extrem, invpi,  mt2,
     +                   n,      np1,    pi,     ptrnod, rcl,    seed )
c     ------------------------------------------------------------------
c     phase1: Construction phase of GRASP to build a permutation of 
c             nodes ordered in a GRASP fashion according to node 
c             degree.
c     ------------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Integer input parameters:
c
c         mt2            - m * 2
c         n              - number of nodes in input graph
c         np1            - n + 1
c
c     ------------------------------------------------------------------
      integer mt2,            n,              np1

c     ------------------------------------------------------------------
c
c     Integer input/output parameters:
c
c         seed           - seed for pseudo random number generator
c
c     ------------------------------------------------------------------
      integer seed

c     ------------------------------------------------------------------
c
c     Integer work parameters:
c
c         cutoff         - cut off to be in RCL
c         exnode         - extreme node
c         i              - do loop parameter
c         ii             - do loop parameter
c         maxdeg         - maximum degree
c         mindeg         - minimum degree
c         node           - node chosen from RCL
c         nselct         - index of RCL of node chosen from RCL
c         numrcl         - number of RCL elements
c         ptr            - pointer
c         xdeg           - degree
c
c     ------------------------------------------------------------------
      integer cutoff,         exnode,         i,
     +        ii,             maxdeg,         mindeg,
     +        node,           nselct,         numrcl,
     +        ptr,            xdeg

c     ------------------------------------------------------------------
c
c     Real input parameter:
c
c         alpha          - RCL parameter
c
c     ------------------------------------------------------------------
      real    alpha

c     ------------------------------------------------------------------
c
c     Real work parameter:
c
c         randp          - random number generator function
c         xrand          - probability returned by randp
c
c     ------------------------------------------------------------------
      real    randp,          xrand

c     ------------------------------------------------------------------
c
c     Integer input arrays:
c
c         deg0(n)        - node degrees of input graph
c         extrem(m*2)    - node in adjcency list
c         ptrnod(n+1)    - pointer to first node in adjacency list
c
c     ------------------------------------------------------------------
      integer deg0(n),        extrem(mt2),    ptrnod(np1)

c     ------------------------------------------------------------------
c
c     Integer work arrays:
c
c         deg(n)         - copy of node degrees of input graph
c         rcl(n)         - restricted candidate list (RCL)
c
c     ------------------------------------------------------------------
      integer deg(n),         rcl(n)

c     ------------------------------------------------------------------
c
c     Integer input/output arrays:
c
c         invpi(n)       - inverse of node permutation
c         pi(n)          - node permutation
c
c     ------------------------------------------------------------------
      integer invpi(n),       pi(n)

c     ------------------------------------------------------------------
c     Make temporary copy of node degrees.
c     ------------------------------------------------------------------
      call cpi4(n,deg0,deg)

c     ------------------------------------------------------------------
c     Find nodes with maximum and minimum degrees.
c     ------------------------------------------------------------------
      mindeg=deg(1)
      maxdeg=deg(1)
      do 10 i =2,n
           xdeg = deg(i)
           maxdeg = max( maxdeg, xdeg )
           mindeg = min( mindeg, xdeg )
10    continue
      
c     ------------------------------------------------------------------
c     Compute cutoff for restricted candidate list (RCL) and build RCL.
c     ------------------------------------------------------------------
      cutoff=alpha*(maxdeg-mindeg)+mindeg
      numrcl=0
      do 20 i =1,n
           if (deg(i).le.cutoff) then
                numrcl=numrcl+1
                rcl(numrcl)=i
           endif
20    continue

c     ------------------------------------------------------------------
c     Select at random from the best numrcl nodes.
c     ------------------------------------------------------------------
      xrand=randp(seed)
      nselct=1+seed/(2147483647/numrcl)
      node=rcl(nselct)

c     ------------------------------------------------------------------
c     Setup permutation and inverse permutation vectors.
c     ------------------------------------------------------------------
      pi(node)=1
      invpi(1)=node

c     ------------------------------------------------------------------
c     Adaptively change the degrees of remaining nodes to take
c     into account node put in permutation.
c     ------------------------------------------------------------------
      deg(node)=-1
      do 30 ptr=ptrnod(node),ptrnod(node+1)-1
           exnode=extrem(ptr)
           deg(exnode)=deg(exnode)-1
30    continue

c     ------------------------------------------------------------------
c     Put 2nd, 3rd, ... last elements in permutation. 
c     ------------------------------------------------------------------
      do 90 i=2,n-1
c          -------------------------------------------------------------
c          Find nodes with maximum and minimum degrees.
c          -------------------------------------------------------------
           mindeg=n+1
           maxdeg=-1
           do 40 ptr=ptrnod(node),ptrnod(node+1)-1
                exnode=extrem(ptr)
                xdeg=deg(exnode)
                if (xdeg.ge.0) then
                     maxdeg = max ( maxdeg, xdeg )
                     mindeg = min ( mindeg, xdeg )
                endif
40         continue
           if (maxdeg.ge.0) then
c               --------------------------------------------------------
c               If (i-1)st node in permutation has adjacent nodes: 
c               compute cutoff for restricted candidate list (RCL) and 
c               build RCL from set of adjacent nodes.
c               --------------------------------------------------------
                cutoff=alpha*(maxdeg-mindeg)+mindeg
                numrcl=0
                do 50 ptr=ptrnod(node),ptrnod(node+1)-1
                     exnode=extrem(ptr)
                     xdeg=deg(exnode)
                     if (xdeg.lt.0) goto 50
                     if (xdeg.le.cutoff) then
                          numrcl=numrcl+1
                          rcl(numrcl)=exnode
                     endif
50              continue
           else
c               --------------------------------------------------------
c               If (i-1)st node in permutation no adjacent nodes: 
c               compute cutoff for restricted candidate list (RCL) and 
c               build RCL from all remaining nodes.
c               --------------------------------------------------------
                mindeg=n+1
                maxdeg=-1
                do 60 ii =1,n
                     xdeg=deg(ii)
                     if (xdeg.ge.0) then
                          maxdeg = max( maxdeg, xdeg )
                          mindeg = min( mindeg, xdeg )
                     endif
60              continue
                cutoff=alpha*(maxdeg-mindeg)+mindeg
                numrcl=0
                do 70 ii=1,n
                     xdeg=deg(ii)
                     if (xdeg.lt.0) goto 70
                     if (xdeg.le.cutoff) then
                          numrcl=numrcl+1
                          rcl(numrcl)=ii
                     endif
70              continue
           endif
      
c          -------------------------------------------------------------
c          Select at random from the best numrcl nodes.
c          -------------------------------------------------------------
           xrand=randp(seed)
           nselct=1+seed/(2147483647/numrcl)
           node=rcl(nselct)
           pi(node)=i
           invpi(i)=node

c          -------------------------------------------------------------
c          Adaptively change the degrees of remaining nodes to take
c          into account node put in permutation.
c          -------------------------------------------------------------
           deg(node)=-1
           do 80 ptr=ptrnod(node),ptrnod(node+1)-1
                exnode=extrem(ptr)
                deg(exnode)=deg(exnode)-1
80         continue
90    continue

c     ------------------------------------------------------------------
c     Adjust pi and invpi.
c     ------------------------------------------------------------------
      do 100 i=1,n
           if (deg(i).ge.0) then
                pi(i)=n
                invpi(n)=i
                return
           endif
100   continue
c     ------------------------------------------------------------------
c     End of phase1.
c     ------------------------------------------------------------------
      return
      end




      real function randp(ix)
c     -----------------------------------------------------------------
c     randp: Portable pseudo-random number generator.
c            Reference: L. Schrage, "A More Portable Fortran
c            Random Number Generator", ACM Transactions on
c            Mathematical Software, Vol. 2, No. 2, (June, 1979).
c     -----------------------------------------------------------------

      implicit none
      integer a,p,ix,b15,b16,xhi,xalo,leftlo,fhi,k
      data a/16807/,b15/32768/,b16/65536/,p/2147483647/

      xhi=ix/b16
      xalo=(ix-xhi*b16)*a
      leftlo=xalo/b16
      fhi=xhi*a+leftlo
      k=fhi/b15
      ix=(((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
      if (ix.lt.0) ix=ix+p

      randp=float(ix)*4.656612875e-10

c     ------------------------------------------------------------------
c     End of phase1.
c     ------------------------------------------------------------------
      return
      end




      subroutine cpi4(n,x,y)
c     ------------------------------------------------------------------
c     cpi4: Copy integer array x of dimension n onto integer array y
c           of dimension n.
c     ------------------------------------------------------------------
c     Passed input parameter:
c
c          n     - dimension of arrays x and y.
c
c     ------------------------------------------------------------------
      implicit none
      integer n
c     ------------------------------------------------------------------
c     Passed input array:
c
c          x     - array
c     ------------------------------------------------------------------
      integer x(n)
c     ------------------------------------------------------------------
c     Passed output array:
c
c          y     - array
c     ------------------------------------------------------------------
      integer y(n)
c     ------------------------------------------------------------------
c     Local scalar:
c
c          i     - do loop counter
c     ------------------------------------------------------------------
      integer i
c     ------------------------------------------------------------------

      do 10 i=1,n
           y(i)=x(i)
10    continue
c     ------------------------------------------------------------------
c     End of cpi4.
c     ------------------------------------------------------------------
      return
      end




      subroutine phase2( arc1,   arc2,   invpi,  m,      mt2,    n, 
     +                   np1,    numarc, pi,     ptrnod )
c     ------------------------------------------------------------------
c     phase2: Local search phase of GRASP. Attempts to decrease the
c             number of crossing arcs, by swapping nodes that are 
c             adjacent in the permutation.
c     ------------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Integer input parameters:
c
c         m              - number of arcs in input graph
c         mt2            - m * 2
c         n              - number of nodes in input graph
c         np1            - n + 1
c
c     ------------------------------------------------------------------
      integer m,              mt2,            n,
     +        np1

c     ------------------------------------------------------------------
c
c     Integer work parameters:
c
c         bstcrs         - best number of crossing arcs 
c         crsedg         - function to compute number of crossing arcs
c                          for permutation pi.
c         goswap         - indicator to re-do 2-exchange
c         i              - do loop parameter
c         piv1           - pi(v1)
c         piv2           - pi(v2)
c         v1             - node
c         v2             - node
c         xg             - gain achieved by swap
c         xgain          - function to compute gain achieved by swap
c
c     ------------------------------------------------------------------
      integer bstcrs,         crsedg,         goswap,
     +        i,              piv1,           piv2,
     +        v1,             v2,             xg,
     +        xgain

c     ------------------------------------------------------------------
c
c     Integer input arrays:
c
c         arc1(m)        - node 1 of arc in input graph
c         arc2(m)        - node 2 of arc in input graph
c         invpi(n)       - inverse of node permutation
c         numarc(m*2)    - arc number
c         pi(n)          - node permutation
c         ptrnod(n+1)    - pointer to first node in adjacency list
c
c     ------------------------------------------------------------------
      integer arc1(m),        arc2(m),        invpi(n),
     +        numarc(mt2),    pi(n),          ptrnod(np1)

c     ------------------------------------------------------------------
c     Compute number of crossing arcs corresponding to permutation pi.
c     ------------------------------------------------------------------
      bstcrs=crsedg( arc1,   arc2,   m,      n,      pi )
c     ------------------------------------------------------------------
c     Loop until no improvement is possible swapping adjacent nodes
c     in permutation pi.
c     ------------------------------------------------------------------
10    goswap=0
      do 20 i=1,n-1
c          -------------------------------------------------------------
c          Nodes to swap are v1 and v2.
c          -------------------------------------------------------------
           v1=invpi(i)
           v2=invpi(i+1)

c          -------------------------------------------------------------
c          Compute gain achieved by swapping v1 and v2.
c          -------------------------------------------------------------
           xg=xgain( arc1,   arc2,   invpi,  m,      mt2,    n,
     +               np1,    numarc, pi,     ptrnod, v1,     v2 )

c          -------------------------------------------------------------
c          If gain is not achieved, unswap positions of v1 and v2.
c          Else, keep swap and update best number of crossing arcs.
c          -------------------------------------------------------------
           if (xg.ge.0) then
                piv1=pi(v1)
                piv2=pi(v2)
                invpi(piv1)=v2
                invpi(piv2)=v1
                pi(v1)=piv2
                pi(v2)=piv1
           else
                bstcrs=bstcrs+xg
                goswap=1
           endif

20    continue
c     ------------------------------------------------------------------
c     If at least one swap reduced number of crossing arcs, re-do
c     pairwise exchanges.
c     ------------------------------------------------------------------
      if (goswap.eq.1) goto 10

c     ------------------------------------------------------------------
c     End of phase2.
c     ------------------------------------------------------------------
      return
      end




      integer function crsedg( arc1,   arc2,   m,      n,      pi )
c     ------------------------------------------------------------------
c     crsedg: Computes the number of crossing arcs in node 
c             permutation defined by pi.
c     ------------------------------------------------------------------
      implicit none
c     ------------------------------------------------------------------
c
c     Integer input parameters:
c
c         m              - number of arcs in input graph
c         n              - number of nodes in input graph
c
c     ------------------------------------------------------------------
      integer m,              n

c     ------------------------------------------------------------------
c
c     Integer work parameters:
c
c         a              - node
c         b              - node
c         c              - node
c         d              - node
c         i              - do loop parameter
c         j              - do loop parameter
c         parc1i         - pi(arc1(i))
c         parc1j         - pi(arc1(j))
c         parc2i         - pi(arc2(i))
c         parc2j         - pi(arc2(j))
c
c     ------------------------------------------------------------------
      integer a,              b,              c,
     +        d,              i,              j,
     +        parc1i,         parc1j,         parc2i,
     +        parc2j

c     ------------------------------------------------------------------
c
c     Integer input arrays:
c
c         arc1(m)        - node 1 of arc in input graph
c         arc2(m)        - node 2 of arc in input graph
c         pi(n)          - node permutation
c
c     ------------------------------------------------------------------
      integer arc1(m),        arc2(m),        pi(n)

c     ------------------------------------------------------------------
c     Initialize counter.
c     ------------------------------------------------------------------
      crsedg=0
c     ------------------------------------------------------------------
c     Scan adjacent nodes in permutation.
c     ------------------------------------------------------------------
      do 20 i=1,m-1
          do 10 j=i+1,m
c               --------------------------------------------------------
c               Edges are (a,b) and (c,d).
c               --------------------------------------------------------
                parc1i=pi(arc1(i))
                parc2i=pi(arc2(i))
                parc1j=pi(arc1(j))
                parc2j=pi(arc2(j))
                a=min(parc1i,parc2i)
                b=max(parc1i,parc2i)
                c=min(parc1j,parc2j)
                d=max(parc1j,parc2j)

c               --------------------------------------------------------
c               If (a,b) and (c,d) cross, increment counter.
c               --------------------------------------------------------
                if ( (a.lt.c .and. c.lt.b .and. b.lt.d) .or.
     +               (c.lt.a .and. a.lt.d .and. d.lt.b) ) then
                     crsedg=crsedg+1
                endif
10        continue
20    continue
c     ------------------------------------------------------------------
c     End of crsedg.
c     ------------------------------------------------------------------
      return
      end




      subroutine swap( invpi,  n,      pi,     v1,     v2 )
c     ------------------------------------------------------------------
c     swap: Swaps v1 and v2 in permutation defined by pi and invpi.
c     ------------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Integer input parameters:
c
c         n              - number of nodes
c         v1             - first node for swapping
c         v2             - second node for swapping
c
c     ------------------------------------------------------------------
      integer n,              v1,             v2

c     ------------------------------------------------------------------
c     Integer work parameters:
c
c         piv1           - pi(v1)
c         piv2           - pi(v2)
c
c     ------------------------------------------------------------------
      integer piv1,           piv2

c     ------------------------------------------------------------------
c     Integer input arrays:
c
c         invpi(n)       - inverse of node permutation
c         pi(n)          - node permutation
c
c     ------------------------------------------------------------------
      integer invpi(n),       pi(n)

c     ------------------------------------------------------------------

      piv1=pi(v1)
      piv2=pi(v2)
      invpi(piv1)=v2
      invpi(piv2)=v1
      pi(v1)=piv2
      pi(v2)=piv1
c     ------------------------------------------------------------------
c     End of swap.
c     ------------------------------------------------------------------
      return
      end




      integer function xgain( arc1,   arc2,   invpi,  m,      mt2,    n,
     +                        np1,    numarc, pi,     ptrnod, v1,     v2
     +                      )
c     ------------------------------------------------------------------
c     xgain: Computes reduction in cross arcs by swapping v1 and v2 in 
c            permutation defined by pi and invpi.
c     ------------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Integer input parameters:
c
c         m              - number of arcs in input graph
c         mt2            - m * 2
c         n              - number of nodes in input graph
c         np1            - n + 1
c
c     ------------------------------------------------------------------
      integer m,              mt2,            n,
     +        np1

c     ------------------------------------------------------------------
c
c     Integer work parameters:
c
c         aftcrs         - number of crossings after swap
c         bfrcrs         - number of crossings before swap
c         piv1           - pi(v1)
c         piv2           - pi(v2)
c         u              - do loop parameter
c         v1             - node
c         v2             - node
c
c     ------------------------------------------------------------------
      integer aftcrs,         bfrcrs,         piv1,
     +        piv2,           u,              v1,
     +        v2

c     ------------------------------------------------------------------
c     Integer input arrays:
c
c         arc1(m)        - node 1 of arc in input graph
c         arc2(m)        - node 2 of arc in input graph
c         invpi(n)       - inverse of node permutation
c         numarc(m*2)    - arc number
c         pi(n)          - node permutation
c         ptrnod(n+1)    - pointer to first node in adjacency list
c
c     ------------------------------------------------------------------
      integer arc1(m),        arc2(m),        invpi(n),
     +        numarc(mt2),    pi(n),          ptrnod(np1)

c     ------------------------------------------------------------------

c     ------------------------------------------------------------------
c     Compute crossings before swap.
c     ------------------------------------------------------------------
      bfrcrs=0
      do 30 u=1,m
           call cntcrs( arc1,   arc2,   bfrcrs, m,      mt2,    n,
     +                   np1,   numarc, pi,     ptrnod, u,      v1)

           if (arc1(u).ne.v1.and.arc2(u).ne.v1) then
                call cntcrs( arc1,   arc2,   bfrcrs, m,      mt2,
     +                       n,      np1,    numarc, pi,     ptrnod,
     +                       u,      v2)
           endif
30    continue

c     ------------------------------------------------------------------
c     Swap nodes v1 and v2 in permutation (pi, invpi).
c     ------------------------------------------------------------------
      piv1=pi(v1)
      piv2=pi(v2)
      invpi(piv1)=v2
      invpi(piv2)=v1
      pi(v1)=piv2
      pi(v2)=piv1

c     ------------------------------------------------------------------
c     Compute crossings after swap.
c     ------------------------------------------------------------------
      aftcrs=0
      do 60 u=1,m
           call cntcrs( arc1,   arc2,   aftcrs, m,      mt2,    n,
     +                   np1,    numarc, pi,     ptrnod, u,      v1)

           call cntcrs( arc1,   arc2,   aftcrs, m,      mt2,    n,
     +                   np1,    numarc, pi,     ptrnod, u,      v2)
60    continue

c     ------------------------------------------------------------------
c     Compute the gain: number of crossings after swap minus number of
c     crossings before swap.
c     ------------------------------------------------------------------
      xgain=aftcrs-bfrcrs

c     ------------------------------------------------------------------
c     End of xgain.
c     ------------------------------------------------------------------
      return
      end




      subroutine cntcrs( arc1,   arc2,   crs,    m,      mt2,    n,
     +                   np1,    numarc, pi,     ptrnod, u,      v)
c     ------------------------------------------------------------------
c     cntcrs: Counts number of edge crossings.
c     ------------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Integer input parameters:
c
c         m              - number of arcs in input graph
c         mt2            - m * 2
c         n              - number of nodes in input graph
c         np1            - n + 1
c         u              - arc number
c         v              - arc number
c
c     ------------------------------------------------------------------
      integer m,              mt2,            n,
     +        np1,            u,              v

c     ------------------------------------------------------------------
c
c     Integer work parameters:
c
c         a              - node
c         arc            - arc
c         b              - node
c         c              - node
c         d              - node
c         parc1a         - pi(arc1(arc))
c         parc2a         - pi(arc2(arc))
c         parc1u         - pi(arc1(u))
c         parc2u         - pi(arc2(u))
c         ptr            - pointer
c
c     ------------------------------------------------------------------
      integer a,              arc,            b,
     +        c,              d,              parc1a,         
     +        parc2a,         parc1u,         parc2u,
     +        ptr

c     ------------------------------------------------------------------
c     Integer input arrays:
c
c         arc1(m)        - node 1 of arc in input graph
c         arc2(m)        - node 2 of arc in input graph
c         numarc(m*2)    - arc number
c         pi(n)          - node permutation
c         ptrnod(n+1)    - pointer to first node in adjacency list
c
c     ------------------------------------------------------------------
      integer arc1(m),        arc2(m),        numarc(mt2),
     +        pi(n),          ptrnod(np1)

c     ------------------------------------------------------------------
c     Integer input/output parameter:
c
c         crs            - number of crossings
c
c     ------------------------------------------------------------------
      integer crs

c     ------------------------------------------------------------------

      crs=0
      do 10 ptr=ptrnod(v),ptrnod(v+1)-1
           arc=numarc(ptr)
c          -------------------------------------------------------------
c          Edges are (a,b) and (c,d).
c          -------------------------------------------------------------
           parc1u=pi(arc1(u))
           parc2u=pi(arc2(u))
           parc1a=pi(arc1(arc))
           parc2a=pi(arc2(arc))
           a=min(parc1u,parc2u)
           b=max(parc1u,parc2u)
           c=min(parc1a,parc2a)
           d=max(parc1a,parc2a)

           if ((a.lt.c .and. c.lt.b .and. b.lt.d) .or.
     +         (c.lt.a .and. a.lt.d .and. d.lt.b)) then
c              ---------------------------------------------------------
c              Edges (a,b) and (c,d) cross.
c              ---------------------------------------------------------
               crs=crs+1
           endif
10    continue

c     ------------------------------------------------------------------
c     End of cntcrs.
c     ------------------------------------------------------------------
      return
      end




      subroutine mkux( arc1,   arc2,   iniux,  m,      m2,     maxux,
     +                 n,      nextux, pi,     ux )
c     ---------------------------------------------------------------
c     mkux:  Constructs the set U(x) of all intervals wholly
c            contained in interval x = (a,b), where a and b are
c            nodes of the permutation defined by pi.
c
c            The set is placed in a singly linked list =
c            {iniux,nextux,ux}.
c     ---------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Integer input parameters:
c
c         m              - number of arcs in input graph
c         m2             - m * m
c         n              - number of nodes in input graph
c
c     ------------------------------------------------------------------
      integer m,              m2,             n

c     ------------------------------------------------------------------
c
c     Integer work parameters:
c
c         a1             - arc1(i)
c         a2             - arc2(i)
c         i              - do loop parameter
c         j              - do loop parameter
c         lft            - node
c         lftj           - node
c         p1             - pi(a1)
c         p2             - pi(a2)
c         p1j            - pi(a1)
c         p2j            - pi(a2)
c         ptr            - pointer
c         rgt            - node
c         rgtj           - node
c
c     ------------------------------------------------------------------
      integer a1,             a2,             i,
     +        j,              lft,            lftj,
     +        p1,             p1j,            p2,
     +        p2j,            ptr,            rgt,
     +        rgtj

c     ------------------------------------------------------------------
c
c     Integer output parameter:
c
c         maxux          - number of elements in set U(x)
c
c     ------------------------------------------------------------------
      integer maxux

c     ------------------------------------------------------------------
c
c     Integer input arrays:
c
c         arc1(m)        - node 1 of arc in input graph
c         arc2(m)        - node 2 of arc in input graph
c         pi(n)          - node permutation
c
c     ------------------------------------------------------------------
      integer arc1(m),        arc2(m),        pi(n)

c     ------------------------------------------------------------------
c
c     Integer output arrays:
c
c         iniux(m+1)     - pointer to first interval in U(x) list
c         nextux(m*m)    - next element of U(x) list
c         ux(m*m)        - set UX for stable set algorithm
c
c     ------------------------------------------------------------------
      integer iniux(m+1),     nextux(m2),     ux(m2)

c     ------------------------------------------------------------------
c     If there are no arcs, do nothing.
c     ------------------------------------------------------------------
      if(m.eq.0) return

c     ------------------------------------------------------------------
c     Initialize linked list that represents the set U(x).
c     ------------------------------------------------------------------
      iniux(m+1)=m
      do 10 i=1,m
           iniux(i)=0
10    continue
      do 20 i=1,m
           ux(i)=i
           nextux(i)=i-1
20    continue

c     ------------------------------------------------------------------
c     Scan all arcs.
c     ------------------------------------------------------------------
      ptr=m+1
      do 40 i=1,m
c          -------------------------------------------------------------
c          Define interval x = (lft,rgt).
c          -------------------------------------------------------------
           a1=arc1(i)
           a2=arc2(i)
           p1=pi(a1)
           p2=pi(a2)
           lft=min(p1,p2)
           rgt=max(p1,p2)
c          -------------------------------------------------------------
c          Identify U(x), arcs wholly included in x.
c          -------------------------------------------------------------
           do 30 j=1,m
c               --------------------------------------------------------
c               Process only arcs that are diffent from those in x.
c               --------------------------------------------------------
                if (i.ne.j) then
c                    ---------------------------------------------------
c                    Define interval x = (lft,rgt).
c                    ---------------------------------------------------
                     a1=arc1(j)
                     a2=arc2(j)
                     p1j=pi(a1)
                     p2j=pi(a2)
                     lftj=min(p1j,p2j)
                     rgtj=max(p1j,p2j)

c                    ---------------------------------------------------
c                    If (lftj,rgtj) is inside (lft,rgt), add interval to
c                    set U(x).
c                    ---------------------------------------------------
                     if (lftj.ge.lft  .and.  rgtj.le.rgt) then
                          nextux(ptr)=iniux(i)
                          ux(ptr)=j
                          iniux(i)=ptr
                          ptr=ptr+1
                     endif
                endif
30         continue
40    continue
c     ------------------------------------------------------------------
c     Compute number of intervals in set U(x).
c     ------------------------------------------------------------------
      maxux=ptr-1
c     ------------------------------------------------------------------
c     End of mkux.
c     ------------------------------------------------------------------
      return
      end




      subroutine mxstbl( adj,    arc1,   arc2,   bigw,   currnt, iniadj,
     +                   iniux,  kset,   m,      m2,     mark,   mp1,
     +                   n,      nextux, nstbl,  nxtadj, pi,     point,
     +                   sset,   stable, stack,  stack2, ux,     vrtex,
     +                   w )
c     ------------------------------------------------------------------
c     mxstbl: Computes a maximum stable set in an overlap graph.
c
c             ref: Golumbic, M.C., "Algorithmic Graph Theory and 
c             Perfect Graphs", Chapters 5 and 11.
c       
c             Non-recursive implementation of:
c
c             procedure MAXSTABLE(UX)
c                begin
c             1.    if UX = {} return {};
c             2.    while (there exists x in UX with w(x) undefined) do
c             3.         S(x) = {x} UNION MAXSTABLE(U(x));
c             4.         w(x) = |S(x)|;
c             5.    endwhile
c             6.    T = MAXCLIQUE (UX);
c             7.    return UNION {v in T} S(v);
c                end
c
c     ------------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Integer input parameters:
c
c         m              - number of arcs in input graph
c         m2             - m * m
c         mp1            - m + 1
c         n              - number of nodes in input graph
c
c     ------------------------------------------------------------------
      integer m,              m2,             mp1,
     +        n              

c     ------------------------------------------------------------------
c
c     Integer work parameters:
c
c         flag           - flag
c         i              - do loop parameter and pointer 
c         ii             - do loop parameter
c         j              - do loop parameter and pointer 
c         jj             - do loop parameter
c         ku             - m*(uxx-1)
c         kv             - m*(v-1)
c         lft            - node
c         nkset          - number of elements in array kset
c         nstack         - number of elements in stack
c         nstbl          - number of stable nodes
c         nv             - node counter
c         ptra           - pointer used to scan adjacency list
c         rgt            - node
c         uxx            - ux(x)
c         v              - kset(i)
c         x              - pointer used to scan UX
c
c     ------------------------------------------------------------------
      integer flag,           i,              ii,
     +        j,              jj,             ku,
     +        kv,             lft,            nkset,
     +        nstack,         nstbl,          nv,
     +        ptra,           rgt,            uxx,
     +        v,              x


c     ------------------------------------------------------------------
c
c     Integer input arrays:
c
c         arc1(m)        - node 1 of arc in input graph
c         arc2(m)        - node 2 of arc in input graph
c         iniux(m+1)     - pointer to first interval in U(x) list
c         nextux(m*m)    - next element of U(x) list
c         ux(m*m)        - set UX for stable set algorithm
c
c     ------------------------------------------------------------------
      integer arc1(m),        arc2(m),        iniux(mp1),
     +        nextux(m2),     ux(m2)

c     ------------------------------------------------------------------
c
c     Integer work arrays:
c
c         adj(m*m)       - adjacency matrix
c         bigw(m)        - set W in exact stable set algorithm
c         currnt(m)      - current node in max stable algorithm
c         iniadj(m+1)    - pointer to first node in adjacency list
c         kset(m)        - set of stable nodes
c         mark(m)        - mark in max stable set algorithm
c         nxtadj(m*m)    - next node in adjacency list
c         pi(n)          - node permutation
c         point(m)       - array POINTER in procedure MAXSTABLE
c         sset(m*m)      - set S(*) of stable nodes (edges in graph)
c         stack(m*m)     - stack to simulate recursive logic
c         stack2(m*m)    - stack to simulate recursive logic
c         vrtex(m)       - nodes of set U(x)
c         w(m)           - weight
c
c     ------------------------------------------------------------------
      integer adj(m2),        bigw(m),        currnt(m),
     +        iniadj(mp1),    kset(m),        mark(m),
     +        nxtadj(m2),     pi(n),          point(m),
     +        sset(m2),       stack(m2),      stack2(m2),
     +        vrtex(m),       w(m) 

c     ------------------------------------------------------------------
c
c     Integer output arrays:
c
c         stable(m)      - set of stable nodes (arcs in input graph)
c
c     ------------------------------------------------------------------
      integer stable(m)

c     ------------------------------------------------------------------
c     If UX = EMPTYSET, then return EMPTYSET.     
c     ------------------------------------------------------------------
      if (m.eq.0) then
           nkset=0
           return
      endif

c     ------------------------------------------------------------------
c     Initialize w(*) and set S(*).
c     ------------------------------------------------------------------
      do 10 i=1,m
           w(i)=0
10    continue

      do 20 i=1,m*m
           sset(i)=0
20    continue

      nstack=0
      x = iniux(m+1)
c     ------------------------------------------------------------------
c     Initialize stack to simulate recursive logic: push 0 into stack.
c     ------------------------------------------------------------------
      nstack=nstack+1
      stack(nstack)=0

      flag = 2
c     ------------------------------------------------------------------
c     Control of recursive logic.
c     ------------------------------------------------------------------
30    if (x .eq. 0) then
           if (flag.eq.1) then
                goto 40
           else
                goto 50
           endif
      endif

c     ------------------------------------------------------------------
c     Search for next x in UX with w(x) undefined.
c     ------------------------------------------------------------------
      uxx=ux(x)
      if (w(uxx) .ne. 0) then
           x=nextux(x)
           flag=0
           goto 30
      endif

c     ------------------------------------------------------------------
c     Compute S(x) = {x} UNION MAXSTABLE(U(x)) and 
c     set w(x) = | S(x) | = 1.
c     ------------------------------------------------------------------
      sset(m*(uxx-1)+uxx)=1
      w(uxx)=1

c     ------------------------------------------------------------------
c     Push x into stack.
c     ------------------------------------------------------------------
      nstack=nstack+1
      stack(nstack)=x

c     ------------------------------------------------------------------
c     Get first element of UX(x).
c     ------------------------------------------------------------------
      x=iniux(uxx)
      flag=1
      goto 30

40    continue
c     ------------------------------------------------------------------
c     Pop x out of stack.
c     ------------------------------------------------------------------
      x=stack(nstack)
      nstack=nstack-1

c     ------------------------------------------------------------------
c     Get next element of UX(x).
c     ------------------------------------------------------------------
      x=nextux(x)
      if (x.ne.0) goto 30

50    continue
c     ------------------------------------------------------------------
c     Pop x out of stack.
c     ------------------------------------------------------------------
      x=stack(nstack)
      nstack=nstack-1

c     ------------------------------------------------------------------
c     Prepare data structures to be used in procedure MAXWEIGHT CLIQUE.
c     ------------------------------------------------------------------
      nv=0
      if (x.eq.0) then
           i=m
      else
           i=iniux(ux(x))
      endif

60    if (i.ne.0) then
           nv=nv+1
           vrtex(nv)=ux(i)
           i=nextux(i)
           goto 60
      endif

      do 70 i=1,m
           iniadj(i)=0
70    continue

      ptra=1
      do 90 ii=1,nv
           i=vrtex(ii)
           rgt=max(pi(arc1(i)),pi(arc2(i)))
           do 80 jj=1,nv
                j=vrtex(jj)
                if (i.eq.j) goto 80
                lft=min(pi(arc1(j)),pi(arc2(j)))
                if (lft.ge.rgt) then
                     nxtadj(ptra)=iniadj(i)
                     adj(ptra)=j
                     iniadj(i)=ptra
                     ptra=ptra+1
                endif
80         continue
           currnt(i)=iniadj(i)
90    continue

c     ------------------------------------------------------------------
c     Compute maximum weighted stable set by computing procedure
c     MAXWEIGHT CLIQUE(UX).
c     ------------------------------------------------------------------
      call maxclq( adj, bigw, currnt, iniadj, kset, m, m2, mark,
     +             mp1, nkset, nv, nxtadj, point, stack2, vrtex, w )

c     ------------------------------------------------------------------
c     If on top of search tree, finish up.
c     ------------------------------------------------------------------
      if (x.eq.0) goto 120

c     ------------------------------------------------------------------
c     If not on top of search tree, get next U(x) and compute
c     S(x) = {x} UNION MAXSTABLE(U(x)) and set w(x) = | S(x) |.
c     ------------------------------------------------------------------
      uxx=ux(x)
      ku=m*(uxx-1)
      do 110 i=1,nkset
           v=kset(i)
           kv=m*(v-1)
           do 100 j=1,m
                if (sset(ku+j).eq.0 .and. sset(kv+j).eq.1) then
                     w(uxx)=w(uxx)+1
                     sset(ku+j)=1
                endif
100        continue
110   continue

      flag=0
      x=nextux(x)
      goto 30

c     ------------------------------------------------------------------
c     Finish up.
c
c     Build UNION {s in T} S(v).  Stable set to be returned is stored
c     in stable(*).
c
c     stable (j) = 1, if node j is in stable set
c                = 0, if node j is not in stable set
c     ------------------------------------------------------------------
120   continue
      nstbl=0
      do 130 i=1,m
           stable(i)=0
130   continue
      do 150 i=1,nkset
           v=kset(i)
           kv=m*(v-1)
           do 140 j=1,m
                if (stable(j).eq.0) then
                     if (sset(kv+j).eq.1) then
                          stable(j)=1
                          nstbl=nstbl+1
                     endif
                endif
140        continue
150   continue

c     ------------------------------------------------------------------
c     End of mxstbl.
c     ------------------------------------------------------------------
      return
      end




      subroutine maxclq( adj,    bigw,   currnt, iniadj, kset,   m,
     +                   m2,     mark,   mp1,    nkset,  nv,     nxtadj,
     +                   point,  stack2, vrtex,  w )
c     ------------------------------------------------------------------
c     maxclq:  The MAXWEIGHT CLIQUE algorithm called by MAXSTABLE.
c
c             ref: Golumbic, M.C., "Algorithmic Graph Theory and 
c             Perfect Graphs", Chapter 5.
c
c             MAXWEIGHT CLIQUE (V)
c                begin
c              1.  for all (v in V) do
c              2.      if (v is unexplored) then
c              3.          EXPLORE(v);
c                  end for all;
c              4.  select y in V such that W(y) = max { W(v) | v in V };
c              5.  K = {y};
c              6.  y = POINTER(y);
c              7.  while (y not nil) do
c              8.     K = K UNION {y};
c              9.     y = POINTER(y);
c                  end while;
c             10.  return (K);  
c                 end
c
c     ------------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Integer input parameters:
c
c         m              - number of arcs in input graph
c         m2             - m * m
c         nv             - number of nodes
c
c     ------------------------------------------------------------------
      integer m,              m2,             mp1,
     +        nv             

c     ------------------------------------------------------------------
c
c     Integer work parameters:
c
c         i              - do loop parameter
c         v              - vrtex(i)
c         wmax           - maximum of bigw()
c         ymax           - index maximum element of bigw()
c
c     ------------------------------------------------------------------
      integer i,              v,              wmax,
     +        ymax

c     ------------------------------------------------------------------
c
c     Integer output parameters:
c
c         nkset          - size of kset()
c
c     ------------------------------------------------------------------
      integer nkset

c     ------------------------------------------------------------------
c
c     Integer input arrays:
c
c         adj(m*m)       - adjacency matrix
c         iniadj(m+1)    - pointer to first node in adjacency list
c         nxtadj(m*m)    - next node in adjacency list
c         vrtex(m)       - nodes of set U(x)
c         w(m)           - weight
c
c     ------------------------------------------------------------------
      integer adj(m2),        iniadj(mp1),     nxtadj(m2),
     +        vrtex(m),       w(m)

c     ------------------------------------------------------------------
c
c     Integer work arrays:
c
c         mark(m)        - mark in max stable set algorithm
c         point(m)       - point in max stable set algorithm
c         stack2(m*m)    -
c
c     ------------------------------------------------------------------
      integer mark(m),        point(m),        stack2(m2)

c     ------------------------------------------------------------------
c
c     Integer input/output arrays:
c
c         bigw(m)        - set W in exact stable set algorithm
c         currnt(m)      - current node in max stable algorithm
c
c     ------------------------------------------------------------------
      integer bigw(m),        currnt(m)      

c     ------------------------------------------------------------------
c
c     Integer output arrays:
c
c         kset(m)        - set of stable nodes
c
c     ------------------------------------------------------------------
      integer kset(m)

c     ------------------------------------------------------------------
c     Set clique size counter to zero and flag mark to indicate that
c     all nodes are unexplored.
c     ------------------------------------------------------------------
      nkset=0
      do 10 i=1,m
           mark(i)=0
10    continue

c     ------------------------------------------------------------------
c     Scan all nodes.
c     ------------------------------------------------------------------
      do 20 i=1,nv
           v=vrtex(i)
c          -------------------------------------------------------------
c          If node vrtex(i) has not yet been explored:
c          -------------------------------------------------------------
           if (mark(v).eq.0) then
c               --------------------------------------------------------
c               Explore it.
c               --------------------------------------------------------
                call explor( adj,    bigw,   currnt, iniadj, m,
     +                       m2,     mark,   mp1,    nxtadj, point,
     +                       stack2, v,      w )
           endif
20    continue
      
c     ------------------------------------------------------------------
c     Select y in V such that W(y) = max { W(v) | v in V }.
c     ------------------------------------------------------------------
      wmax=0
      do 30 i=1,nv
           v=vrtex(i)
           if (bigw(v) .gt. wmax) then
                wmax=bigw(v)
                ymax=v
           endif
30    continue

c     ------------------------------------------------------------------
c     While y is not NIL, do:
c     ------------------------------------------------------------------
40    continue
c          -------------------------------------------------------------
c          Add set {y} to set K.
c          -------------------------------------------------------------
           nkset=nkset+1
           kset(nkset)=ymax

c          -------------------------------------------------------------
c          y = POINTER(y).
c          -------------------------------------------------------------
           ymax=point(ymax)

           if (ymax.ne.0) goto 40
      return
c     ------------------------------------------------------------------
c     End of maxclq.
c     ------------------------------------------------------------------
      end




      subroutine explor( adj, bigw, currnt, iniadj, m, m2, mark, mp1,
     +                   nxtadj, point, stack2, v, w )
c     ------------------------------------------------------------------
c     explor: Subroutine EXPLORE called from MAXWEIGHT CLIQUE.
c
c             ref: Golumbic, M.C., "Algorithmic Graph Theory and 
c                  Perfect Graphs", Chapter 5.
c
c             recursive procedure EXPLORE(v)
c                begin     
c        1.         if ( Adj(v) = EMPTYSET ) then
c        2.              W(v) = w(v);
c        3.              POINTER(v) = NIL;
c        4.              return;
c                   end if;
c        5.         for all (x in Adj(v)) do
c        6.             if ( x is unexplored ) then
c        7.                  EXPLORE(x);
c                       end if;
c                   end for all;
c        8.         select y in Adj(v) such that W(y) =
c                          max { W(x) | x in Adj(v) };
c        9.         W(v) = w(v) + W(y);
c       10.         POINTER(v) = y;
c       11.         return
c                end
c
c     ------------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Integer input parameters:
c
c         m              - number of arcs in input graph
c         m2             - m * m
c         mp1            - m + 1
c         v              - node
c
c     ------------------------------------------------------------------
      integer m,              m2,             mp1,
     +        v

c     ------------------------------------------------------------------
c
c     Integer work parameters:
c
c         i              - do loop parameter
c         nstack         - number of stack elements
c         wmax           - largest element of array bigw()
c         x              - node
c         y              - node
c         ymax           - index of largest element of array bigw()        
c
c     ------------------------------------------------------------------
      integer i,              nstack,         wmax,
     +        x,              y,              ymax

c     ------------------------------------------------------------------
c
c     Integer input arrays:
c
c         adj(m*m)       - adjacency matrix
c         iniadj(m+1)    - pointer to first node in adjacency list
c         nxtadj(m*m)    - next node in adjacency list
c         w(m)           - weight
c
c     ------------------------------------------------------------------
      integer adj(m2),        iniadj(mp1),    nxtadj(m2),
     +        w(m)

c     ------------------------------------------------------------------
c
c     Integer input/output arrays:
c
c         bigw(m)        - set W in exact stable set algorithm
c         currnt(m)      - current node in max stable algorithm
c         mark(m)        - mark in max stable set algorithm
c         point(m)       - point in max stable set algorithm
c         stack2(m*m)    - stack used to simulate recursive logic
c
c     ------------------------------------------------------------------
      integer bigw(m),        currnt(m),      mark(m),
     +        point(m),       stack2(m2)

c     ------------------------------------------------------------------
c     ------------------------------------------------------------------
c     Initilize stack and mark v as explored.
c     ------------------------------------------------------------------
      nstack=0
      mark(v)=1

c     ------------------------------------------------------------------
c     Push v into stack2.
c     ------------------------------------------------------------------
      nstack=nstack+1
      stack2(nstack)=v

c     ------------------------------------------------------------------
c     If stack is empty, return control to MAXWEIGHT CLIQUE.
c     ------------------------------------------------------------------
10    if(nstack.eq.0) return

c     ------------------------------------------------------------------
c     Get node v from stack.
c     ------------------------------------------------------------------
      v=stack2(nstack)

c     ------------------------------------------------------------------
c     If v has no adjacent nodes:
c     ------------------------------------------------------------------
      if(iniadj(v).eq.0) then
c          -------------------------------------------------------------
c          Compute W(v) = w(v) and POINTER(v) = NIL.
c          -------------------------------------------------------------
           bigw(v) = w(v)
           point(v)=0
c          -------------------------------------------------------------
c          Pop v out of stack2.
c          -------------------------------------------------------------
           v=stack2(nstack)
           nstack=nstack-1
           go to 10
      else
c          -------------------------------------------------------------
c          Recursively explores nodes adjacent to v.
c          -------------------------------------------------------------
           y = currnt(v)
           if(y.ne.0) then
                x = adj(y)

c               --------------------------------------------------------
c               if node x is unexplored:
c               --------------------------------------------------------
                if(mark(x).eq.0) then
c                    ---------------------------------------------------
c                    Mark node x as explored.
c                    ---------------------------------------------------
                     mark(x) = 1

c                    ---------------------------------------------------
c                    Push v into stack2.
c                    ---------------------------------------------------
                     nstack=nstack+1
                     stack2(nstack)=x
                endif
                currnt(v) = nxtadj(currnt(v))
                go to 10
           else
c               --------------------------------------------------------
c               Select y in Adj(v) such that W(y) = 
c                              max { W(x) | x in Adj(v) }.
c               --------------------------------------------------------
                wmax=0
                i=iniadj(v)
20              if (i.eq.0) goto 30
                     if (bigw(adj(i)) .gt. wmax) then
                        ymax=adj(i)
                        wmax=bigw(ymax)
                     endif
                     i=nxtadj(i)
                     goto 20
30              continue
c               --------------------------------------------------------
c               Compute W(v) = w(v) + W(y) and set POINTER(v) = y.
c               --------------------------------------------------------
                bigw(v)=w(v)+wmax
                point(v)=ymax
c               --------------------------------------------------------
c               Pop v out of stack2.
c               --------------------------------------------------------
                v=stack2(nstack)
                nstack=nstack-1
                go to 10
           endif
      endif
c     ------------------------------------------------------------------
c     End of explor.
c     ------------------------------------------------------------------
      end




      subroutine mkblue ( arc1,   arc1r,  arc2,   arc2r,  blue,
     +                    invred, m,      mred,   mred2,  nblue,
     +                    stable )
c     -------------------------------------------------------------
c     mkblue: Color stable set blue, remove blue arcs from 
c             graph, and compress graph for next phase where 
c             arcs will be colored red.
c     -------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Integer input parameter:
c
c         m              - number of arcs in input graph
c
c     ------------------------------------------------------------------
      integer m

c     ------------------------------------------------------------------
c
c     Integer output parameters:
c
c         mred           - number of arcs in reduced graph
c         mred2          - mred * mred
c         nblue          - number of blue arcs
c
c     ------------------------------------------------------------------
      integer mred,           mred2,          nblue

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
c         arc1(m)        - node 1 of arc in input graph
c         arc2(m)        - node 2 of arc in input graph
c         stable(m)      - set of stable nodes (arcs in input graph)
c
c     ------------------------------------------------------------------
      integer arc1(m),        arc2(m),        stable(m)
      
c     ------------------------------------------------------------------
c
c     Integer output arrays:
c
c         arc1r(m)       - node 1 of arc in reduced graph
c         arc2r(m)       - node 2 of arc in reduced graph
c         blue(m)        - set of blue arcs
c         invred(m)      - maps nodes in compressed graph to input
c
c     ------------------------------------------------------------------
      integer arc1r(m),       arc2r(m),       blue(m),        
     +        invred(m)

c     -------------------------------------------------------------
c     Color arcs in stable set blue.
c     -------------------------------------------------------------
      nblue = 0
      do 10 i=1,m
           if (stable(i).eq.1) then
               nblue=nblue+1
               blue(nblue) = i
           endif
10    continue

c     -------------------------------------------------------------
c     Remove blue arcs from graph and compress graph for
c     next phase where arcs will be colored red.
c     -------------------------------------------------------------
      mred = 0
      do 20 i = 1,m
           if(stable(i).ne.1) then
                mred = mred + 1
                arc1r(mred) = arc1(i)
                arc2r(mred) = arc2(i)
                invred(mred)=i
           endif
20    continue
c     ------------------------------------------------------------------
c     Compute dimension for next call to mkux.
c     ------------------------------------------------------------------
      mred2=mred*mred

c     ------------------------------------------------------------------
c     End of mkblue.
c     ------------------------------------------------------------------
      return
      end




      subroutine mkred( invred, m,      mred,   npale,  nred,    pale, 
     +                  red,    stable )
c     -------------------------------------------------------------
c     mkred: Color stable set red and remaining arcs pale.
c     -------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Integer input parameters:
c
c         m              - number of arcs in input graph
c         mred           - number of arcs in reduced graph
c
c     ------------------------------------------------------------------
      integer m,              mred

c     ------------------------------------------------------------------
c
c     Integer output parameters:
c
c         npale          - number of pale arcs
c         nred           - number of red arcs
c
c     ------------------------------------------------------------------
      integer nred,           npale

c     ------------------------------------------------------------------
c
c     Integer work parameters:
c
c         i              - do loop parameter
c
c     ------------------------------------------------------------------
      integer i

c     ------------------------------------------------------------------
c
c     Integer input arrays:
c
c         invred(m)      - maps nodes in compressed graph to input
c         stable(m)      - set of stable nodes (arcs in input graph)
c
c     ------------------------------------------------------------------
      integer invred(m),      stable(m)
  
c     ------------------------------------------------------------------
c
c     Integer output arrays:
c
c         pale(m)        - set of pale arcs
c         red(m)         - set of red arcs
c
c     ------------------------------------------------------------------
      integer pale(m),        red(m)

c     ------------------------------------------------------------------
c     Initialize counters.
c     ------------------------------------------------------------------
      nred  = 0
      npale = 0

c     ------------------------------------------------------------------
c     Scan all arcs.
c     ------------------------------------------------------------------
      do 10 i=1,mred
c          -------------------------------------------------------------
c          If arc is in stable set ...
c          -------------------------------------------------------------
           if (stable(i).eq.1) then
c               --------------------------------------------------------
c               Color arc red.
c               --------------------------------------------------------
                nred=nred+1
                red(nred)=invred(i)
           else
c               --------------------------------------------------------
c               Color arc pale.
c               --------------------------------------------------------
                npale=npale+1
                pale(npale)=invred(i)
           endif
10    continue

c     ------------------------------------------------------------------
c     End of mkred.
c     ------------------------------------------------------------------
      return
      end




      subroutine postp( arc1,   arc2,   avggs,  blue,   gs,     lstblu,
     +                  lstpal, m,      n,      nblue,  npale,  nred,
     +                  pale,   pi,     red )
c     ------------------------------------------------------------------
c
c     postp: Post processing on pale 2-coloring.
c
c            Look for a blue arc that can be colored red and a pale arc
c            that can be colored blue.  If such arcs are present recolor
c            those arcs and increase size of planar subgraph.
c          
c     ------------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Integer input parameters:
c
c         m              - number of arcs in input graph
c         n              - number of nodes in input graph
c
c     ------------------------------------------------------------------
      integer m,              n

c     ------------------------------------------------------------------
c
c     Integer work parameters:
c
c         a              - node
c         b              - node
c         c              - node
c         d              - node
c         e              - node
c         f              - node
c         i              - do loop parameter
c         j              - do loop parameter
c         k              - do loop parameter
c         l              - do loop parameter
c         nltblu         - number of elements in list lstblu
c         nltpal         - number of elements in list lstpal
c         p1             - pi(arc1(*))
c         p2             - pi(arc2(*))
c         u              - pale(i)
c         v              - blue(j) or lstblu(j) or lstpal(j)
c         w              - red(k)
c
c     ------------------------------------------------------------------
      integer a,              b,              c,
     +        d,              e,              f,
     +        i,              j,              k,
     +        l,              nltblu,         nltpal,
     +        p1,             p2,             u,
     +        v,              w

c     ------------------------------------------------------------------
c
c     Integer input/output parameter:
c
c         gs             - size of planar subgraph
c
c     ------------------------------------------------------------------
      integer gs

c     ------------------------------------------------------------------
c
c     Integer output parameters:
c
c         nblue          - number of blue arcs
c         npale          - number of pale arcs
c         nred           - number of red arcs
c
c     ------------------------------------------------------------------
      integer nblue,          npale,          nred

c     ------------------------------------------------------------------
c
c     Integer input arrays:
c
c         arc1(m)        - node 1 of arc in input graph
c         arc2(m)        - node 2 of arc in input graph
c         pi(n)          - node permutation
c
c     ------------------------------------------------------------------
      integer arc1(m),        arc2(m),        pi(n)

c     ------------------------------------------------------------------
c
c     Integer work arrays:
c
c         lstblu(m)      - list of blue arcs that can be colored red
c         lstpal(m)      - list of pale arcs that can be colored blue
c
c     ------------------------------------------------------------------
      integer lstblu(m),      lstpal(m)

c     ------------------------------------------------------------------
c
c     Integer input/output arrays:
c
c         blue(m)        - set of blue arcs
c         pale(m)        - set of pale arcs
c         red(m)         - set of red arcs
c
c     ------------------------------------------------------------------
      integer blue(m),        pale(m),        red(m)

c     ------------------------------------------------------------------
c
c     Double precision input/output parameter:
c
c         avggs          - sum of solutions found
c
c     ------------------------------------------------------------------
      double precision avggs

c     ------------------------------------------------------------------
c     Scan pale arcs.
c     ------------------------------------------------------------------
10    nltpal=0
      do 70 i=1,npale
           u=pale(i)
c          -------------------------------------------------------------
c          Scan blue arcs.
c          -------------------------------------------------------------
           nltblu=0
           do 30 j=1,nblue
                v=blue(j)
                p1=pi(arc1(u))
                p2=pi(arc2(u))
                a=min(p1,p2)
                b=max(p1,p2)
                p1=pi(arc1(v))
                p2=pi(arc2(v))
                c=min(p1,p2)
                d=max(p1,p2)
                if ((a.lt.c .and. c.lt.b .and. b.lt.d) .or.
     +              (c.lt.a .and. a.lt.d .and. d.lt.b)) then
c                    ---------------------------------------------------
c                    If pale and blue arcs intersect, scan red arcs.
c                    ---------------------------------------------------
                     do 20 k=1,nred
                          w=red(k)
                          p1=pi(arc1(w))
                          p2=pi(arc2(w))
                          e=min(p1,p2)
                          f=max(p1,p2)
c                         ----------------------------------------------
c                         If red and blue arcs intersect, nothing done,
c                         since blue cannot be colored red to make space 
c                         for pale arc.
c                         ----------------------------------------------
                          if ((c.lt.e .and. e.lt.d .and. d.lt.f) .or.
     +                        (e.lt.c .and. c.lt.f .and. f.lt.d)) then
                               goto 70
                          endif
20                   continue
                     nltblu=nltblu+1
                     lstblu(nltblu)=v
                endif
30         continue
c          -------------------------------------------------------------
c          Blue arcs in lstblu can be colored red and pale arc colored 
c          blue, thus increasing size of 2-coloring by 1.
c
c          Make blue arcs red and pale arc blue.
c          -------------------------------------------------------------
           do 60 j=1,nltblu
                v=lstblu(j)
                k=1
40              if (blue(k).eq.v) then
                     nred=nred+1
                     red(nred)=v
                     do 50 l=k,nblue-1
                          blue(l)=blue(l+1)
50                   continue
                     nblue=nblue-1
                     goto 60
                else
                     k=k+1
                     goto 40
                endif
60         continue 
           nblue=nblue+1
           blue(nblue)=u
           nltpal=nltpal+1
           lstpal(nltpal)=u
70    continue
      if (nltpal.ne.0) then
           do 100 j=1,nltpal
                v=lstpal(j)
                k=1
80              if (pale(k).eq.v) then
                     do 90 l=k,npale-1
                          pale(l)=pale(l+1)
90                   continue
                     npale=npale-1
                     goto 100
                else
                     k=k+1
                     goto 80
                endif
100        continue 
           goto 10
      else
c          -------------------------------------------------------------
c          Compute size after post optimization.  
c          Accumulate for averaging.
c          -------------------------------------------------------------
           gs = nblue + nred
           avggs=avggs+gs
           return
      endif
c     ------------------------------------------------------------------
c     End of postp.
c     ------------------------------------------------------------------
      end




      subroutine savsol( blue,   gs,     invpi,  iter,   m,      maxgs,
     +                   n,      nblue,  nopblu, nopred, nred,   optblu,
     +                   optitr, optprm, optred, red )
c     ------------------------------------------------------------------
c     savsol: Saves best solution found so far in maxgs, optblu, optred,
c             nopblu, nopred, optprm.
c     ------------------------------------------------------------------
      implicit none

c     ------------------------------------------------------------------
c
c     Integer input parameters:
c
c         gs             - size of planar subgraph
c         iter           - GRASP iteration counter
c         m              - number of arcs in input graph
c         n              - number of nodes in input graph
c         nblue          - number of blue arcs
c         nred           - number of red arcs
c
c     ------------------------------------------------------------------
      integer gs,             iter,           m,
     +        n,              nblue,          nred
      
c     ------------------------------------------------------------------
c
c     Integer output parameters:
c
c         maxgs          - size of largest planar subgraph found
c         nopblu         - number of blue arcs in best planar subgraph
c         optitr         - GRASP iteration best solution was found
c         nopred         - number of red arcs in best planar subgraph
c
c     ------------------------------------------------------------------
      integer maxgs,          nopblu,         nopred,
     +        optitr

c     ------------------------------------------------------------------
c
c     Integer input arrays:
c
c         blue(m)        - set of blue arcs
c         invpi(n)       - inverse of node permutation
c         red(m)         - set of red arcs
c
c     ------------------------------------------------------------------
      integer blue(m),        invpi(n),       red(m)

c     ------------------------------------------------------------------
c
c     Integer output arrays:
c
c         optblu(m)      - set of best blue arcs
c         optprm(n)      - node permutation of best solution
c         optred(m)      - set of best red arcs
c
c     ------------------------------------------------------------------
      integer optblu(m),      optprm(n),      optred(m)

c     ------------------------------------------------------------------
c     Copy parameters.
c     ------------------------------------------------------------------
      nopblu=nblue
      nopred=nred
      maxgs=gs
      optitr=iter

c     ------------------------------------------------------------------
c     Copy arrays.
c     ------------------------------------------------------------------
      call cpi4(nblue,blue,optblu)
      call cpi4(nred,red,optred)
      call cpi4(n,invpi,optprm)

c     ------------------------------------------------------------------
c     End of savsol.
c     ------------------------------------------------------------------
      return
      end

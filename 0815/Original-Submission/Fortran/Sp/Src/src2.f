      subroutine gfvs( alpha,  bcset,  bestz,  cset,   cset0,  din,
     +                 din0,   dout,   dout0,  errcnd, fin0,   fout0,
     +                 gitrb,  lin,    lin0,   look4,  lout,   lout0,
     +                 maxe,   maxitr, maxv,   max2e,  ne,     nin,
     +                 nin0,   nout,   nout0,  nrcl,   nv,     out,
     +                 pin,    pin0,   pout,   pout0,  prttyp, rcl,
     +                 seedb,  seedd,  vertex, vin,    vin0,   vout,
     +                 vout0,  vtx1,   vtx2 )
c     ------------------------------------------------------------------
c
c     gfvs: A package of FORTRAN subroutines for finding feedback vertex
c           sets in a graph using GRASP.
c
c     authors:       Paola Festa             [paofes@udsab.dia.unisa.it]
c                    Mauricio G.C. Resende   [mgcr@research.att.com]
c                    Panos P.M. Pardalos     [pardalos@ufl.edu]
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          alpha  - GRASP RCL parameter
c          look4  - size of cutset sought
c          max2e  - array dimension = 2*maxe
c          maxe   - array dimension
c          maxitr - maximum number of GRASP iterations
c          maxv   - array dimension
c          ne     - number of arcs in original graph
c          nv     - number of vertices in original graph 
c          out    - Fortran output device
c          prttyp - iteration summary report
c                   = 0 - silent
c                   = 1 - show improvements only
c                   = 2 - show all iterations
c
c     ------------------------------------------------------------------
      real       alpha
      integer    look4,  max2e,  maxe,   maxitr, maxv,   ne,
     +           nv,     prttyp, out

c     ------------------------------------------------------------------
c
c     Passed input/output scalars:
c
c          errcnd - error condition
c          fin    - number of incoming arcs (copy)
c          fin0   - number of incoming arcs (original)
c          fout   - number of outgoing arcs (copy)
c          fout0  - number of outgoing arcs (original)
c          ncset  - size of cutset currently found
c          nrcl   - number of elements in RCL
c          seed   - seed for pseudo random number generator
c          seedd  - seed for pseudo random number generator
c
c     ------------------------------------------------------------------
      integer          errcnd, fin,    fin0,   fout,   fout0,  ncset,
     +                 nrcl,   seed
      double precision seedd

c     ------------------------------------------------------------------
c
c     Passed output scalars:
c
c          bestz  - number of vertices in best cutset found
c          gitrb  - number of iterations to find best cutset
c          seedb  - seed at start of best iteration 
c
c     ------------------------------------------------------------------
      integer    bestz,  gitrb,  seedb

c     ------------------------------------------------------------------
c
c     Passed input/output arrays:
c
c           bcset   - best cutset found
c           cset    - cutset currently found (copy)
c           cset0   - cutset currently found
c           din     - in-degree of vertex (copy)
c           din0    - in-degree of vertex
c           dout    - out-degree of vertex (copy)
c           dout0   - out-degree of vertex
c           lin     - pointer to end of in-adjacency list (copy)
c           lin0    - pointer to end of in-adjacency list (original)
c           lout    - pointer to end of out-adjacency list (copy)
c           lout0   - pointer to end of out-adjacency list (original)
c           nin     - pointer to next element of in-adjacency
c                     list (copy)
c           nin0    - pointer to next element of in-adjacency
c                     list (original)
c           nout    - pointer to next element of out-adjacency
c                     list (copy)
c           nout0   - pointer to next element of out-adjacency
c                     list (original)
c           pin     - pointer to start of in-adjacency list (copy)
c           pin0    - pointer to start of in-adjacency list (original)
c           pout    - pointer to start of out-adjacency list (copy)
c           pout0   - pointer to start of out-adjacency list (original)
c           rcl     - Restricted Candidate List
c           vertex  - existence of vertex in the reduced graph
c           vin     - vertex in in-adjacency list (copy)
c           vin0    - vertex in in-adjacency list (original)
c           vout    - vertex in out-adjacency list (copy)
c           vout0   - vertex in out-adjacency list (original)
c           vtx1    - vertex 1 of arc
c           vtx2    - vertex 2 of arc
c
c     ------------------------------------------------------------------
      integer    bcset(maxv),    cset(maxv),     cset0(maxv),
     +           din(maxv),      din0(maxv),     dout(maxv),
     +           dout0(maxv),    lin(max2e),     lin0(max2e),
     +           lout(max2e),    lout0(max2e),   nin(max2e),
     +           nin0(max2e),    nout(max2e),    nout0(max2e),
     +           pin(maxv),      pin0(maxv),     pout(maxv),
     +           pout0(maxv),    rcl(maxv),      vertex(maxv),
     +           vin(max2e),     vin0(max2e),    vout(max2e),
     +           vout0(max2e),   vtx1(maxe),     vtx2(maxe)

c     ------------------------------------------------------------------
c
c     Local scalars:
c
c          alpha0 - GRASP RCL parameter (copy)
c          gitr   - do loop index
c          k      - do loop index
c          n2e    - array dimension = 2*ne
c          randp  - portable pseudo-random number generator
c          seed0  - seed for pseudo random number generator (copy)
c
c     ------------------------------------------------------------------
      real       alpha0, randp
      integer    gitr,   k,      n2e,    seed0

c     ------------------------------------------------------------------
c     Checks parameters and problem dimension.
c     ------------------------------------------------------------------
      call chkfvs( errcnd, look4,  maxe,   maxitr, maxv,   ne,
     +             nv,     prttyp, seed,   seedd,  vtx1,   vtx2 )
      
      if (errcnd .gt.0) return

c     ------------------------------------------------------------------
c     Compute array dimension.
c     ------------------------------------------------------------------
      n2e=2*ne

c     ------------------------------------------------------------------
c     Make data structures.
c     ------------------------------------------------------------------
      call mkds( din0,   dout0,  fin0,   fout0,  lin0,   lout0, n2e,
     +           ne,     nin0,   nout0,  nv,    pin0,    pout0, vertex,
     +           vin0,   vout0,  vtx1,   vtx2 )

c     ------------------------------------------------------------------
c     Warm up random number generator.
c     ------------------------------------------------------------------
      call warmup( seed )

c     ------------------------------------------------------------------
c     Do maxitr GRASP iterations.
c
c     Initialize value of best solution found.
c     ------------------------------------------------------------------
      bestz=nv+1

c     ------------------------------------------------------------------
c     Copy original alpha value.
c     ------------------------------------------------------------------
      alpha0=alpha

      if (prttyp .ge. 1) then
           write(out,10)
10        format(//,' GRASP for Feedback Set Problem-----------------',
     +           //)
      endif

      do 50 gitr=1,maxitr
c          -------------------------------------------------------------
c          Save current seed in seed0.
c          -------------------------------------------------------------
           seed0=seed

c          -------------------------------------------------------------
c          Copy original graph onto work data structure.
c          -------------------------------------------------------------
           call cpds( din,    din0,   dout,   dout0,  fin,    fin0,
     +                fout,   fout0,  lin,    lin0,   lout,   lout0,
     +                n2e,    ne,     nin,    nin0,   nout,   nout0,
     +                nv,     pin,    pin0,   pout,   pout0,  vertex,
     +                vin,    vin0,   vout,   vout0 )

c          -------------------------------------------------------------
c          Generate random alpha if alpha is not fixed.          
c          -------------------------------------------------------------
           if (alpha0 .lt. 0) then
                alpha=randp( seed )
           endif

c          -------------------------------------------------------------
c          Do GRASP construction phase.
c          -------------------------------------------------------------
           call build( alpha,  cset,   din,    dout,   lin,    lout,
     +                 ncset,  n2e,    nin,    nout,   nrcl,   nv,
     +                 pin,    pout,   rcl,    seed,   vertex, vin,
     +                 vout )

c          -------------------------------------------------------------
c          Do GRASP local search phase.
c          -------------------------------------------------------------
           call local( cset,   cset0,  din,    din0,   dout,   dout0,
     +                 fin0,   fout0,  lin,    lin0,   lout,   lout0,
     +                 ncset,  n2e,    ne,     nin,    nin0,   nout,
     +                 nout0,  nv,     pin,    pin0,   pout,   pout0,
     +                 vertex, vin,    vin0,   vout,   vout0 )

c          -------------------------------------------------------------
c          Check if the local search has found a better solution.
c          -------------------------------------------------------------
           if (ncset.lt.bestz) then
c               --------------------------------------------------------
c               Update the variables related to the best solution 
c               currently found.
c               --------------------------------------------------------
                bestz=ncset
                seedb=seed0
                gitrb=gitr
                do 20 k=1,bestz
                     bcset(k)=cset(k)
20              continue
c               --------------------------------------------------------
c               The output option parameter is set in order to print out 
c               the solution improvements. 
c               --------------------------------------------------------
                if (prttyp .ge. 1) then
                     write(out,30) gitr, bestz
30                   format('     itr = ',i8,
     +                      '     cutset size = ',i6,' ***')
                endif
c               --------------------------------------------------------
c               Check if the cardinality of the cutset found is less 
c               than the size sought. 
c               --------------------------------------------------------
                if (ncset.le.look4) then
                     errcnd=0
                     return
                endif
           else
c               --------------------------------------------------------
c               The output option parameter is set in order to print out 
c               the solution found in each GRASP iteration. 
c               --------------------------------------------------------
                if (prttyp .eq. 2) then
                     write(out,40) gitr,ncset
40                   format('     itr = ',i8,
     +                      '     cutset size = ',i6)
                endif
           endif
50    continue

c     ------------------------------------------------------------------
c     Execution terminated successfully.
c     ------------------------------------------------------------------
      errcnd=0
      return

c     ------------------------------------------------------------------
c     End of subroutine gfvs.
c     ------------------------------------------------------------------
      end





      subroutine adde( index,  lxxx,   n2e,    nv,     nxxx,   pxxx,
     +                 u,      vxxx,   w )
c     ------------------------------------------------------------------
c
c     adde: Add arc(index) pointed (u,w) into adjacency list of u.
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          index  - arc index
c          n2e    - array dimension = 2*ne
c          nv     - array dimension
c          u      - vertex index
c          w      - vertex index

c     ------------------------------------------------------------------
      integer    index,  nv,     n2e,    u,      w

c     ------------------------------------------------------------------
c
c     Passed input/output arrays:
c
c          lxxx   - pointer to end of xxx-adjacency list
c          nxxx   - pointer to next element of xxx-adjacency list 
c          pxxx   - pointer to start of xxx-adjacency list 
c          vxxx   - vertex in xxx-adjacency list 
c
c     ------------------------------------------------------------------
      integer    lxxx(n2e),      nxxx(n2e),      pxxx(nv),
     +           vxxx(n2e)

c     ------------------------------------------------------------------
c     Add w into the out-adjacent list of u.
c     ------------------------------------------------------------------
      vxxx(index)=w

c     ------------------------------------------------------------------
c     Set the next arc to the arc index equal to the first outgoing arc
c     from the vertex u.
c     ------------------------------------------------------------------
      nxxx(index)=pxxx(u)

c     ------------------------------------------------------------------
c     If the first outgoing arc is the last element of list, set last
c     of next to element.
c     ------------------------------------------------------------------
      if (pxxx(u).ne.0) lxxx(pxxx(u))=index

c     ------------------------------------------------------------------
c     Set the arc index to be the first and not last outgoing arc 
c     from u.
c     ------------------------------------------------------------------
      pxxx(u)=index
      lxxx(index)=0
      return

c     ------------------------------------------------------------------
c     End of subroutine adde.
c     ------------------------------------------------------------------
      end





      subroutine build( alpha,  cset,   din,    dout,   lin,    lout,
     +                  ncset,  n2e,    nin,    nout,   nrcl,   nv,
     +                  pin,    pout,   rcl,    seed,   vertex, vin,
     +                  vout )
c     ------------------------------------------------------------------
c
c     build: GRASP construction phase (phase 1)
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          alpha  - GRASP RCL parameter
c          n2e    - array dimension = 2*ne
c          nv     - number of nodes in reduced graph
c
c     ------------------------------------------------------------------
      real       alpha
      integer    n2e,    nv

c     ------------------------------------------------------------------
c
c     Passed input/output scalars:
c
c          ncset  - size of cutset found
c          nrcl   - number of elements in RCL
c          seed   - seed for pseudo random number generator
c
c     ------------------------------------------------------------------
      integer    ncset,  nrcl,   seed

c     ------------------------------------------------------------------
c
c     Passed input arrays:
c
c          din    - in-degree of vertex (copy)
c          dout   - out-degree of vertex (copy)
c          lin    - pointer to end of in-adjacency list (copy)
c          lout   - pointer to end of out-adjacency list (copy)
c          nin    - pointer to next element of in-adjacency
c                   list (copy)
c          nout   - pointer to next element of out-adjacency
c                   list (copy)
c          pin    - pointer to start of in-adjacency list (copy)
c          pout   - pointer to start of out-adjacency list (copy)
c          vertex - existence of vertex in the reduced graph
c          vin    - vertex in in-adjacency list (copy)
c          vout   - vertex in out-adjacency list (copy)
c
c     ------------------------------------------------------------------
      integer    din(nv),        dout(nv),       lin(n2e),
     +           lout(n2e),      nin(n2e),       nout(n2e),
     +           pin(nv),        pout(nv),       vertex(nv),
     +           vin(n2e),       vout(n2e)

c     ------------------------------------------------------------------
c
c     Passed output arrays:
c
c          cset   - set of feedback vertices
c
c     ------------------------------------------------------------------
      integer    cset(nv)

c     ------------------------------------------------------------------
c
c     Passed working arrays:
c
c          rcl    - Restricted Candidate List
c
c     ------------------------------------------------------------------
      integer    rcl(nv)

c     ------------------------------------------------------------------
c
c     Local scalars:
c
c          nvleft - number of unselected vertices
c          v      - vertex
c
c     ------------------------------------------------------------------
      integer    nvleft, v

c     ------------------------------------------------------------------
c     Initialize number of vertices left (nvleft) and number of vertices
c     in cut set (ncset).
c     ------------------------------------------------------------------
      nvleft=nv
      ncset=0

c     ------------------------------------------------------------------
c     Apply initial graph reduction. 
c     ------------------------------------------------------------------
      call redg( cset,   din,    dout,   lin,    lout,   ncset,  n2e,
     +           nin,    nout,   nv,     nvleft, pin,    pout,   vertex,
     +           vin,    vout )

c     ------------------------------------------------------------------
c     Begin construction loop
c     ------------------------------------------------------------------
10    continue
           if (nvleft.eq.0) goto 20
c          -------------------------------------------------------------
c          If there are still vertices to be considered,
c          build restricted candidate list (RCL).
c          -------------------------------------------------------------
           call mkrcl( alpha,  din,    dout,   nrcl,   nv,     rcl,
     +                 vertex )

c          -------------------------------------------------------------
c          Pick a candidate from the RCL at random and place v into
c          cut set.
c          -------------------------------------------------------------
           call pickv( cset,   ncset,  nrcl,   nv,     rcl,    seed,
     +                 v )

c          -------------------------------------------------------------
c          Reduce graph by deleting v and its incident arcs from graph.
c          -------------------------------------------------------------
           call delv( din,    dout,   lin,    lout,   n2e,    nin,
     +                nout,   nv,     nvleft, pin,    pout,   v,  
     +                vertex, vin,    vout )

c          -------------------------------------------------------------
c          Apply graph reductions.
c          -------------------------------------------------------------
           call redg( cset,   din,    dout,   lin,    lout,   ncset,
     +                n2e,    nin,    nout,   nv,     nvleft, pin,
     +                pout,   vertex, vin,    vout )

c          -------------------------------------------------------------
c          End loop, and go to check if there still are vertices to be 
c          removed.
c          -------------------------------------------------------------
           goto 10
20    continue
      return

c     ------------------------------------------------------------------
c     End of subroutine build
c     ------------------------------------------------------------------
      end





      subroutine chkfvs( errcnd, look4,  maxe,   maxitr, maxv,   ne,    
     +                   nv,     prttyp, seed,   seedd,  vtx1,   vtx2 )
c     ------------------------------------------------------------------
c
c     chkfvs: Checks values of parameters and problem dimension for
c             validity.  Feedback vertex set case.
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          look4  - size of cutset sought
c          maxe   - array dimension
c          maxitr - maximum number of GRASP iterations
c          maxv   - array dimension
c          ne     - number of arcs in original graph
c          nv     - number of vertices in original graph
c          prttyp - iteration summary report
c                   = 0 - silent
c                   = 1 - show improvements only
c                   = 2 - show all iterations
c          seedd  - seed for pseudo random number generator
c
c     ------------------------------------------------------------------
      integer          look4,  maxe,   maxitr, maxv,   ne,     nv,
     +                 prttyp
      double precision seedd

c     ------------------------------------------------------------------
c
c     Passed output scalars:
c
c          errcnd - error condition     
c          seed   - seed for pseudo random number generator
c
c     ------------------------------------------------------------------
      integer    errcnd, seed

c     ------------------------------------------------------------------
c
c     Passed input arrays:
c
c          vtx1  - vertex 1 of arc
c          vtx2  - vertex 2 of arc
c
c     ------------------------------------------------------------------
      integer    vtx1(maxe),     vtx2(maxe)

c     ------------------------------------------------------------------
c
c     Local scalars:
c
c          i      - do loop index
c
c     ------------------------------------------------------------------
      integer    i

c     ------------------------------------------------------------------
c     Check if the number of vertices of the input graph exceeds maxv.
c     ------------------------------------------------------------------
      if (nv.gt.maxv) then
           errcnd=1
           return
      endif

c     ------------------------------------------------------------------
c     Check if the number of arcs of the input graph exceeds maxe.
c     ------------------------------------------------------------------
      if (ne.gt.maxe) then
           errcnd=2
           return
      endif

c     ------------------------------------------------------------------
c     Check if the vertices indices are well-defined.
c     ------------------------------------------------------------------
      do 10 i=1,ne
           if (vtx1(i).lt.1 .or. vtx1(i).gt.nv) then
                errcnd=3
                return
           endif
           if (vtx2(i).lt.1 .or. vtx2(i).gt.nv) then
                errcnd=3
                return
           endif
10    continue

c     ------------------------------------------------------------------
c     Check if the stopping parameter look4 is well-defined.
c     ------------------------------------------------------------------
      if (look4.lt.0 .or. look4.gt.nv) then
           errcnd=4
           return
      endif

c     ------------------------------------------------------------------
c     Check if the maximum number of GRASP iteration is well-defined.
c     ------------------------------------------------------------------
      if (maxitr.lt.1) then
           errcnd=5
           return
      endif

c     ------------------------------------------------------------------
c     Check if the output option parameter is well-defined.
c     ------------------------------------------------------------------
      if ((prttyp.ne.0) .and. (prttyp.ne.1) .and. (prttyp.ne.2) ) then
           errcnd=6
           return
      endif

c     ------------------------------------------------------------------
c     Check if the pseudo random number generator seed is well-defined.
c     ------------------------------------------------------------------
      if (seedd.lt.1 .or. seedd.gt.2147483647) then
           errcnd=7
           return
      else
           seed = idint(seedd)
      endif
      return

c     ------------------------------------------------------------------
c     End of subroutine chkfvs.
c     ------------------------------------------------------------------
      end





      subroutine cpds( din,    din0,   dout,   dout0,  fin,    fin0,
     +                 fout,   fout0,  lin,    lin0,   lout,   lout0,
     +                 n2e,    ne,     nin,    nin0,   nout,   nout0,
     +                 nv,     pin,    pin0,   pout,   pout0,  vertex,
     +                 vin,    vin0,   vout,   vout0 )
c     ------------------------------------------------------------------
c
c     cpds: Copy original data structure into work data structure.
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          fin0   - number of incoming arcs (original)
c          fout0  - number of outgoing arcs (original)
c          n2e    - array dimension = 2*ne
c          ne     - number of arcs in graph
c          nv     - number of nodes in graph
c
c     ------------------------------------------------------------------
      integer    fin0,   fout0,  n2e,    ne,     nv

c     ------------------------------------------------------------------
c
c     Passed input arrays:
c
c         din0    - in-degree of vertex
c         dout0   - out-degree of vertex
c         lin0    - pointer to end of in-adjacency list (original)
c         lout0   - pointer to end of out-adjacency list (original)
c         nin0    - pointer to next element of in-adjacency
c                   list (copy)
c         nout0   - pointer to next element of out-adjacency
c                   list (original)
c         pin0    - pointer to start of in-adjacency list (original)
c         pout0   - pointer to start of out-adjacency list (original)
c         vertex  - existence of vertex in the reduced graph
c         vin0    - vertex in in-adjacency list (original)
c         vout0   - vertex in out-adjacency list (original)
c
c     ------------------------------------------------------------------
      integer    din0(nv),       dout0(nv),      lin0(n2e),
     +           lout0(n2e),     nin0(n2e),      nout0(n2e),
     +           pin0(nv),       pout0(nv),      vertex(nv),
     +           vin0(n2e),      vout0(n2e)

c     ------------------------------------------------------------------
c
c     Passed output scalars:
c
c         fin     - number of incoming arcs (copy)
c         fout    - number of outgoing arcs (copy)
c
c     ------------------------------------------------------------------
      integer    fin,    fout

c     ------------------------------------------------------------------
c
c     Passed output arrays:
c
c         din    - in-degree of vertex (copy)
c         dout   - out-degree of vertex (copy)
c         lin    - pointer to end of in-adjacency list (copy)
c         lout   - pointer to end of out-adjacency list (copy)
c         nin    - pointer to next element of in-adjacency
c                  list (copy)
c         nout   - pointer to next element of out-adjacency
c                  list (copy)
c         pin    - pointer to start of in-adjacency list (copy)
c         pout   - pointer to start of out-adjacency list (copy)
c         vin    - vertex in in-adjacency list (copy)
c         vout   - vertex in out-adjacency list (copy)
c
c     ------------------------------------------------------------------
      integer    din(nv),        dout(nv),       lin(n2e),
     +           lout(n2e),      nin(n2e),       nout(n2e),
     +           pin(nv),        pout(nv),       vin(n2e),
     +           vout(n2e)

c     ------------------------------------------------------------------
c
c     Local scalars:
c    
c          i     - do loop index
c    
c     ------------------------------------------------------------------
      integer    i

c     ------------------------------------------------------------------
c     cpi4( n,x,y ): Copy the n-dimensional array x into 
c                    the n-dimensional array y
c     ------------------------------------------------------------------
      n2e=2*ne
      call cpi4( nv,     din0,   din  )
      call cpi4( nv,     dout0,  dout )
      call cpi4( n2e,    lin0,   lin  )
      call cpi4( n2e,    lout0,  lout )
      call cpi4( n2e,    nin0,   nin  )
      call cpi4( n2e,    nout0,  nout )
      call cpi4( nv,     pin0,   pin  )
      call cpi4( nv,     pout0,  pout )
      call cpi4( n2e,    vin0,   vin  )
      call cpi4( n2e,    vout0,  vout )

c     ------------------------------------------------------------------
c     Initialize all nodes to existing in reduced graph.
c     ------------------------------------------------------------------
      do 10 i=1,nv
           vertex(i)=1
10    continue

c     ------------------------------------------------------------------
c     Copy number of incoming and outgoing arcs.
c     ------------------------------------------------------------------
      fin=fin0
      fout=fout0
      return

c     ------------------------------------------------------------------
c     End of subroutine cpds.
c     ------------------------------------------------------------------
      end





      subroutine cpi4( n,      x,      y )
c     ------------------------------------------------------------------
c
c     cpi4:  Copy the n-dimensional array x into the n-dimensional 
c            array y
c
c     ------------------------------------------------------------------
c
c     Passed input scalar:
c
c          n   - array dimension
c
c     ------------------------------------------------------------------
      integer    n
c     ------------------------------------------------------------------
c
c     Passed input array:
c    
c          x   - array to be copied 
c
c     ------------------------------------------------------------------
      integer    x(n)

c     ------------------------------------------------------------------
c
c     Passed output array:
c    
c          y   - copy array  
c
c     ------------------------------------------------------------------
      integer    y(n)

c     ------------------------------------------------------------------
c
c     Local scalar:
c    
c          i      - do loop index
c    
c     ------------------------------------------------------------------
      integer    i

c     ------------------------------------------------------------------
      do 10 i=1,n
           y(i)=x(i)
10    continue
      return

c     ------------------------------------------------------------------
c     End of subroutine cpi4.
c     ------------------------------------------------------------------
      end





      subroutine dele( lxxx,   n2e,    nv,     nxxx,   pxxx,   q,      
     +                 u )
c     ------------------------------------------------------------------
c
c     dele: Delete arc pointed to by q from adjacency list of u.
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          n2e    - array dimension = 2*ne
c          nv     - array dimension
c          q      - pointer to arc
c          u      - vertex index
c
c     ------------------------------------------------------------------
      integer    n2e,    nv,     q,      u

c     ------------------------------------------------------------------
c
c     Passed input/output arrays:
c
c          lxxx   - pointer to end of xxx-adjacency list
c          nxxx   - pointer to next element of xxx-adjacency list
c          pxxx   - pointer to start of xxx-adjacency list
c
c     ------------------------------------------------------------------
      integer    lxxx(n2e),      nxxx(n2e),      pxxx(nv)

c     ------------------------------------------------------------------
c
c     Local scalars:
c
c          lq     - pointer to vertex
c          nq     - pointer to vertex
c
c     ------------------------------------------------------------------
      integer    lq,     nq

c     ------------------------------------------------------------------
c     Initialize nq and lq to be next and previous element pointed to by
c     q in adjacency list of u.
c     ------------------------------------------------------------------
      nq=nxxx(q)
      lq=lxxx(q)

c     ------------------------------------------------------------------
c     Repeat while list is not scanned.
c     ------------------------------------------------------------------
      if (nq.eq.0) then
           if (lq.ne.0) then
c               --------------------------------------------------------
c               Last element of list, not first.
c               --------------------------------------------------------
                nxxx(lq)=0
           else
c               --------------------------------------------------------
c               Only element of list.
c               --------------------------------------------------------
                pxxx(u)=0
           endif
      else
           if (lq.ne.0) then
c               --------------------------------------------------------
c               Not last element of list, not first.
c               --------------------------------------------------------
                nxxx(lq)=nq
                lxxx(nq)=lq
           else
c               --------------------------------------------------------
c               First, but not last element of list.
c               --------------------------------------------------------
                pxxx(u)=nq
                lxxx(nq)=0
           endif
      endif
      return

c     ------------------------------------------------------------------
c     End of subroutine dele.
c     ------------------------------------------------------------------
      end





      subroutine delv( din,    dout,   lin,    lout,   n2e,    nin,
     +                 nout,   nv,     nvleft, pin,    pout,   v,
     +                 vertex, vin,    vout )
c     ------------------------------------------------------------------
c
c     delv:  Remove v and all its incident arcs from graph.
c            Adapt greedy function by changing degrees.
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          n2e    - array dimension = 2*ne
c          nv     - array dimension
c          nvleft - number of unselected vertices
c          v      - vertex
c
c     ------------------------------------------------------------------
      integer    n2e,    nv,     nvleft, v

c     ------------------------------------------------------------------ 
c
c     Passed input/output arrays:
c
c          din     - in-degree of vertex (copy)
c          dout    - out-degree of vertex (copy)
c          lin     - pointer to end of in-adjacency list (copy)
c          lout    - pointer to end of out-adjacency list (copy)
c          nin     - pointer to next element of in-adjacency
c                    list (copy)
c          nout    - pointer to next element of out-adjacency
c                    list (copy)
c          pin     - pointer to start of in-adjacency list (copy)
c          pout    - pointer to start of out-adjacency list (copy)
c          vertex  - existence of vertex in the reduced graph
c          vin     - vertex in in-adjacency list (copy)
c          vout    - vertex in out-adjacency list (copy)
c
c     ------------------------------------------------------------------
      integer    din(nv),        dout(nv),       lin(n2e),
     +           lout(n2e),      nin(n2e),       nout(n2e),
     +           pin(nv),        pout(nv),       vertex(nv),
     +           vin(n2e),       vout(n2e)

c     ------------------------------------------------------------------
c
c     Local scalars:
c
c          p      - pointer to arc
c          q      - pointer to arc
c          u      - vertex
c          w      - vertex
c
c     ------------------------------------------------------------------
      integer    p,      q,      u,      w

c     ------------------------------------------------------------------
c     Remove from graph all arcs incoming into vertex v.
c     
c     Let p be the first incoming arc in v.
c     ------------------------------------------------------------------
      p=pin(v)

c     ------------------------------------------------------------------
c     Test if v has no incoming arc.
c     ------------------------------------------------------------------
10    if (p.eq.0) goto 40

c          -------------------------------------------------------------
c          Let u be the tail of the arc p.
c          -------------------------------------------------------------
           u=vin(p)

c          -------------------------------------------------------------
c          If u does not belong to the reduced graph, then let p be 
c          the next incoming arc.
c          -------------------------------------------------------------
           if (vertex(u).eq.0) then
                p=nin(p)
                goto 10
           endif 

c          -------------------------------------------------------------
c          If p is not a self-loop arc, then let q be the first 
c          outgoing arc from u and let w be its head.
c          -------------------------------------------------------------
           if (u.ne.v) then
                q=pout(u)
20              w=vout(q)

c               --------------------------------------------------------
c               Test if w is vertex v.
c               -------------------------------------------------------- 
                if (w.eq.v) then
c                   ----------------------------------------------------
c                   Delete arc pointed to by q from adjacency list of u.
c                   ----------------------------------------------------
                    call dele( lout,   n2e,    nv,     nout,   pout,
     +                         q,      u )

c                   ----------------------------------------------------
c                   Decrease the out-degree of vertex u.
c                   ----------------------------------------------------
                    dout(u)=dout(u)-1
                    goto 30
                else
c                   ----------------------------------------------------
c                   Let q be the next outgoing arc from u.
c                   ----------------------------------------------------
                    q=nout(q)
                    goto 20
                endif
           endif
30         continue

c          -------------------------------------------------------------
c          Let p be the next arc incoming into v.
c          -------------------------------------------------------------
           p=nin(p)
           goto 10
40    continue
      
c     ------------------------------------------------------------------
c     Vertex v has no incoming arcs.
c     ------------------------------------------------------------------
      pin(v)=0

c     ------------------------------------------------------------------
c     Remove from graph all outgoing arcs from vertex v.
c     
c     Let p be the first outgoing arc from v.
c     ------------------------------------------------------------------
      p=pout(v)

c     ------------------------------------------------------------------
c     Test if v has no outgoing arcs.
c     ------------------------------------------------------------------
50    if (p.eq.0) goto 80

c          -------------------------------------------------------------
c          Let u be the head of arc p.
c          -------------------------------------------------------------
           u=vout(p)

c          -------------------------------------------------------------
c          If u does not belong to the reduced graph, then let p be 
c          the next outgoing arc.
c          -------------------------------------------------------------
           if (vertex(u).eq.0) then
                p=nout(p)
                goto 50
           endif 

c          -------------------------------------------------------------
c          If p is not a self-loop arc, then let q be the first incoming 
c          arc of u and w be its tail.
c          -------------------------------------------------------------
           if (u.ne.v) then
                q=pin(u)
60              w=vin(q)

c               --------------------------------------------------------
c               Test if w is vertex v.
c               --------------------------------------------------------
                if (w.eq.v) then
c                   ----------------------------------------------------
c                   Delete arc pointed to by q from adjacency list of u.
c                   ----------------------------------------------------
                    call dele( lin,    n2e,    nv,     nin,    pin,
     +                         q,      u )

c                   ----------------------------------------------------
c                   Decrease the in-degree of vertex u.
c                   ----------------------------------------------------
                    din(u)=din(u)-1
                    goto 70
                else
c                   ----------------------------------------------------
c                   Let q be the next incoming arc of u.
c                   ----------------------------------------------------
                    q=nin(q)
                    goto 60
                endif
           endif
70         continue
c          -------------------------------------------------------------
c          Let p be the next outgoing arc from v.
c          -------------------------------------------------------------
           p=nout(p)
           goto 50
80    continue

c     ------------------------------------------------------------------
c     Vertex v has no outgoing arcs.
c     ------------------------------------------------------------------
      pout(v)=0

c     ------------------------------------------------------------------
c     Remove vertex v from the graph.
c     ------------------------------------------------------------------
      vertex(v)=0
      nvleft=nvleft-1
      return

c     ------------------------------------------------------------------
c     End of subroutine delv.
c     ------------------------------------------------------------------
      end





      subroutine dflset( alpha,  look4,  maxitr, prttyp, seedd )
c     ------------------------------------------------------------------
c
c     dflset: Sets default values for parameters.
c
c     ------------------------------------------------------------------
c
c     Passed output scalars:
c
c          alpha  - RCL parameter
c          look4  - size of cutset sought
c          maxitr - maximum number of GRASP iterations
c          prttyp - iteration summary report
c                   = 0 - silent
c                   = 1 - show improvements only
c                   = 2 - show all iterations
c          seedd  - seed for pseudo random number generator
c 
c     ------------------------------------------------------------------
      real             alpha
      integer          look4,  maxitr, prttyp
      double precision seedd

c     ------------------------------------------------------------------
      alpha =-1.0
      seedd =270001.d0
      maxitr=2048
      prttyp=1
      look4 =0
c     ------------------------------------------------------------------
      return

c     ------------------------------------------------------------------
c     End of subroutine dflset.
c     ------------------------------------------------------------------
      end





      subroutine drop( cset,   i,      ncset0, nv )
c     ------------------------------------------------------------
c
c     drop: Drops i-th element from the current solution set.
c
c     ------------------------------------------------------------
c
c     Passed input scalars:
c
c          i      - vertex index
c          ncset0 - size of cutset found (copy)
c          nv     - number of nodes in graph
c
c     ------------------------------------------------------------
      integer    i,      ncset0, nv

c     ------------------------------------------------------------
c
c     Passed working arrays:
c
c          cset0  - cutset currently found
c 
c     ------------------------------------------------------------
      integer    cset(nv)

c     ------------------------------------------------------------
c
c     Local scalars:
c
c          j      - do loop index
c
c     ------------------------------------------------------------
      integer    j

c     ------------------------------------------------------------
c     Check if the i-th element is the last inserted in the 
c     current solution set.
c     ------------------------------------------------------------
      if (i.eq.ncset0) then
c          -------------------------------------------------------
c          Decrease the cardinality of the current solution set.
c          -------------------------------------------------------
           ncset0=ncset0-1
      else
c          -------------------------------------------------------
c          Shift one position left the last (ncset0-1-i) elements 
c          of the current solution set.
c          ------------------------------------------------------- 
           do 10 j=i,ncset0-1
                cset(j)=cset(j+1)
10         continue

c          -------------------------------------------------------
c          Decrease the cardinality of the current solution set.
c          -------------------------------------------------------
           ncset0=ncset0-1
      endif
      return

c     ------------------------------------------------------------------
c     End of subroutine drop.
c     ------------------------------------------------------------------
      end





      subroutine errmsg( errcnd, out )
c     ------------------------------------------------------------------
c
c     errmsg: Prints out error message corresponding to error condition
c             given by errcnd (returned from gfvs).
c
c     ------------------------------------------------------------------
c
c     Passed input scalar:
c
c       errcnd - error condition number
c
c          successful termination
c
c                 = 0    - found target look4 in fewer than maxitr
c                          iterations
c
c          unsuccessful termination
c
c                 = 1    - nv too big
c                 = 2    - ne too big 
c                 = 3    - input vertex label invalid 
c                 = 4    - cutset size too small or too big
c                 = 5    - maxitr too small 
c                 = 6    - prttyp invalid 
c                 = 7    - invalid seed 
c
c
c          out    - Fortran output device
c
c     ------------------------------------------------------------------
      integer    errcnd, out 

c     ------------------------------------------------------------------
      if (errcnd .eq. 1) then
           write(out,10)
10         format(' error: number of vertices too big. Increase maxv.')
      endif

      if (errcnd .eq. 2) then 
           write(out,20)
20         format(' error: number of arcs too big. Increase maxe.')
      endif

      if (errcnd .eq. 3) then
           write(out,30)
30         format(' error: node number is not valid.')
      endif

      if (errcnd .eq. 4) then
           write(out,40)
40         format(' error: value of look4 is not valid.')
      endif

      if (errcnd .eq. 5) then
           write(out,50)
50         format(' error: value of maxitr is not valid.')
      endif

      if (errcnd .eq. 6) then
           write(out,60)
60         format(' error: value of prttyp is not valid.')
      endif

      if(errcnd .eq. 7) then
           write(out,70)
70         format(' error: value of seed is not valid.')
      endif

      if(errcnd .eq. 0) then
           write(out,80)
80         format(//,'     Execution terminated with no error.')
      endif
      return

c     ------------------------------------------------------------------
c     End of subroutine errmsg.
c     ------------------------------------------------------------------
      end





      subroutine local( cset,   cset0,  din,    din0,   dout,   dout0,
     +                  fin0,   fout0,  lin,    lin0,   lout,   lout0,
     +                  ncset,  n2e,    ne,     nin,    nin0,   nout,
     +                  nout0,  nv,     pin,    pin0,   pout,   pout0,
     +                  vertex, vin,    vin0,   vout,   vout0 )

c     ------------------------------------------------------------------
c
c     local: A local search procedure to improve the solution obtained 
c	     by the GRASP construction phase.
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          n2e    - array dimension = 2*ne
c          ne     - number of arcs in reduced graph
c          nv     - number of nodes in reduced graph
c
c     ------------------------------------------------------------------
      integer    n2e,    ne,     nv

c     ------------------------------------------------------------------
c
c     Passed input/output scalars:
c
c          fin    - number of incoming arcs (copy)
c          fin0   - number of incoming arcs (original)
c          fout   - number of outgoing arcs (copy)
c          fout0  - number of outgoing arcs (original)
c          ncset  - size of cutset found
c          ncset0 - size of cutset found (copy)
c
c     ------------------------------------------------------------------
      integer    fin,    fin0,   fout,   fout0,  ncset,  ncset0

c     ------------------------------------------------------------------
c
c     Passed input arrays:
c
c          din    - in-degree of vertex (copy)
c          dout   - out-degree of vertex (copy)
c          lin    - pointer to end of in-adjacency list (copy)
c          lout   - pointer to end of out-adjacency list (copy)
c          nin    - pointer to next element of in-adjacency
c                   list (copy)
c          nout   - pointer to next element of out-adjacency
c                   list (copy)
c          pin    - pointer to start of in-adjacency list (copy)
c          pout   - pointer to start of out-adjacency list (copy)
c          vertex - existence of vertex in the reduced graph
c          vin    - vertex in in-adjacency list (copy)
c          vout   - vertex in out-adjacency list (copy)
c
c     ------------------------------------------------------------------
      integer    din(nv),        dout(nv),       lin(n2e),
     +           lout(n2e),      nin(n2e),       nout(n2e),
     +           pin(nv),        pout(nv),       vertex(nv),
     +           vin(n2e),       vout(n2e)

c     ------------------------------------------------------------------
c
c     Passed output arrays:
c
c          cset  - set of feedback vertices
c
c     ------------------------------------------------------------------
      integer    cset(nv)

c     ------------------------------------------------------------------
c
c     Passed working arrays:
c
c          cset0  - cutset currently found (copy)
c          din0   - in-degree of vertex
c          dout0  - out-degree of vertex
c          lin0   - pointer to end of in-adjacency list (original)
c          lout0  - pointer to end of out-adjacency list (original)
c          nin0   - pointer to next element of in-adjacency
c                   list (copy)
c          nout0  - pointer to next element of out-adjacency
c                   list (original)
c          pin0   - pointer to start of in-adjacency list (original)
c          pout0  - pointer to start of out-adjacency list (original)
c          vin0   - vertex in in-adjacency list (original)
c          vout0  - vertex in out-adjacency list (original)
c
c     ------------------------------------------------------------------
      integer    cset0(nv),      din0(nv),       dout0(nv),
     +           lin0(n2e),      lout0(n2e),     nin0(n2e),
     +           nout0(n2e),     pin0(nv),       pout0(nv),
     +           vin0(n2e),      vout0(n2e)

c     ------------------------------------------------------------------
c
c     Local scalars:
c
c          i      - do loop index
c          ib     - start of i
c          j      - do loop index
c          nrcl   - number of elements in RCL
c          nvleft - number of unselected vertices
c          v      - vertex
c
c     ------------------------------------------------------------------
      integer    i, ib,     j,      nvleft, v

c     ------------------------------------------------------------------
c     Make copy of the current solution.
c     ------------------------------------------------------------------
      do 10 i=1,ncset
           cset0(i)=cset(i)
10    continue
      ncset0=ncset
      
      do 20 v=1,nv
           vertex(v)=1
20    continue

      ib = 1
30    continue
           do 50 i=ib,ncset0
                nvleft=nv
                ncset=0
c               ------------------------------------------------------
c               Copy original graph onto work data structure.
c               ------------------------------------------------------
                call cpds( din,    din0,   dout,   dout0,  fin,
     +                     fin0,   fout,   fout0,  lin,    lin0,
     +                     lout,   lout0,  n2e,    ne,     nin,
     +                     nin0,   nout,   nout0,  nv,     pin,
     +                     pin0,   pout,   pout0,  vertex, vin,
     +                     vin0,   vout,   vout0 )

                do 40 j=1,ncset0
c                    -------------------------------------------------
c                    Produce reduced graph.
c                    -------------------------------------------------
                     if(j.ne.i) then
                          call delv( din,    dout,   lin,    lout,
     +                               n2e,    nin,    nout,   nv,
     +                               nvleft, pin,    pout,   cset0(j),
     +                               vertex, vin,    vout )
                     endif
40              continue
c               --------------------------------------------------------
c               Apply reduction on the reduced graph.
c               --------------------------------------------------------
                call redg( cset,   din,    dout,   lin,    lout,
     +                     ncset,  n2e,    nin,    nout,   nv,
     +                     nvleft, pin,    pout,   vertex, vin,
     +                     vout )

                if((nvleft+ncset).eq.0) then
c                    ---------------------------------------------------
c                    Drop i-th element from the current solution set
c                    cset0.
c                    ---------------------------------------------------
                     call drop( cset0,  i,      ncset0, nv )
                     ib = i
                     goto 30
                endif
50         continue

c     ------------------------------------------------------------------
c     Copy the found solution onto original variables.
c     ------------------------------------------------------------------
      ncset=ncset0
      call cpi4( ncset0, cset0,  cset )

      return

c     ------------------------------------------------------------------
c     End of subroutine local.
c     ------------------------------------------------------------------
      end





      subroutine mkds( din0,   dout0,  fin0,   fout0,  lin0,   lout0,
     +                 n2e,    ne,     nin0,   nout0,  nv,     pin0,
     +                 pout0,  vertex, vin0,   vout0,  vtx1,   vtx2 )
c     ------------------------------------------------------------------
c
c     mkds: Builds adjacency data structure for graph, computes degrees.
c
c     ------------------------------------------------------------------
c
c     Passed input parameters:
c
c          n2e    - array dimension = 2*ne
c          ne     - number of arcs in graph
c          nv     - number of nodes in graph
c
c     ------------------------------------------------------------------
      integer    n2e,    ne,     nv

c     ------------------------------------------------------------------
c
c     Passed input arrays:
c    
c          vtx1  - node 1 of arc
c          vtx2  - node 2 of arc
c
c     ------------------------------------------------------------------
      integer    vtx1(ne),       vtx2(ne)

c     ------------------------------------------------------------------
c
c     Passed output arrays:
c    
c          din0   - in-degree of vertex (original) 
c          dout0  - out-degree of vertex (original)
c          lin0   - pointer to end of in-adjacency list (original)
c          lout0  - pointer to end of out-adjacency list (original)
c          nin0   - pointer to next vertex of in-adjacency 
c                   list (original)
c          nout0  - pointer to next vertex of out-adjacency 
c                   list (original)
c          pin0   - pointer to first vertex of in-adjacency 
c                   list (original)
c          pout0  - pointer to first vertex of out-adjacency 
c                   list (original)
c          vertex - existence of vertex in the reduced graph
c          vin0   - vertex in in-adjacency list (original)
c          vout0  - vertex in out-adjacency list (original) 
c
c     ------------------------------------------------------------------
      integer    din0(nv),       dout0(nv),      lin0(n2e),
     +           lout0(n2e),     nin0(n2e),      nout0(n2e),
     +           pin0(nv),       pout0(nv),      vertex(nv),
     +           vin0(n2e),      vout0(n2e)

c     ------------------------------------------------------------------
c
c     Local scalars:
c    
c          i      - do loop index
c          v1     - node 1
c          v2     - node 2
c          fin0   - number of incoming arcs (original)
c          fout0  - number of outgoing arcs (original)
c    
c     -----------------------------------------------------------------
      integer    i,      v1,     v2,     fin0,   fout0

c     ------------------------------------------------------------------
c     Initialize pointers to list of adjacent vertices and degrees.
c     ------------------------------------------------------------------
      do 10 i=1,nv
          pin0(i)=0
          pout0(i)=0
          din0(i)=0
          dout0(i)=0
          vertex(i)=1
10    continue

      do 20 i=1,n2e
          nin0(i)=i+1
          nout0(i)=i+1
          lin0(i)=0
          lout0(i)=0
          vin0(i)=0
          vout0(i)=0
20    continue

      nin0(n2e)=0
      nout0(n2e)=0
      fin0=1
      fout0=1

      do 30 i=1,ne
c          -------------------------------------------------------------
c          Get directed edge = (v1,v2).
c          -------------------------------------------------------------
           v1=vtx1(i)
           v2=vtx2(i)
c          -------------------------------------------------------------
c          Add v2 to v1's out-adjacency list.
c          -------------------------------------------------------------
           vout0(fout0)=v2
           nout0(fout0)=pout0(v1)
           if (pout0(v1).ne.0) then
                lout0(pout0(v1))=fout0
           endif
           lout0(fout0)=0
           pout0(v1)=fout0
           dout0(v1)=dout0(v1)+1
           fout0=fout0+1
c          -------------------------------------------------------------
c          Add v1 to v2's in-adjacency list.
c          -------------------------------------------------------------
           vin0(fin0)=v1
           nin0(fin0)=pin0(v2)
           if (pin0(v2).ne.0) then
                lin0(pin0(v2))=fin0
           endif
           lin0(fin0)=0
           pin0(v2)=fin0
           din0(v2)=din0(v2)+1
           fin0=fin0+1
30    continue
      return

c     ------------------------------------------------------------------
c     End of subroutine mkds.
c     ------------------------------------------------------------------
      end





      subroutine mkrcl( alpha,  din,    dout,   nrcl,   nv,     rcl,
     +                  vertex )
c     ------------------------------------------------------------------
c
c     mkrcl: Build the restricted candidate list (RCL).  
c
c            All vertices v with 
c
c            in-degree(v) * out-degree(v) >= g + alpha * (G - g)
c
c            are put in the RCL, where
c
c            g = min (over all v) in-degree(v) * out-degree(v),
c            G = max (over all v) in-degree(v) * out-degree(v),      
c            and
c            alpha is a real number such that 0 < alpha < 1.
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          alpha  - GRASP RCL parameter
c          nv     - number of nodes in graph
c
c     ------------------------------------------------------------------
      real       alpha
      integer    nv

c     ------------------------------------------------------------------
c
c     Passed input arrays:
c   
c          din    - in-degree of vertex (copy)
c          dout   - out-degree of vertex (copy) 
c          vertex - existence of vertex in the reduced graph
c
c     ------------------------------------------------------------------
      integer    din(nv),        dout(nv),       vertex(nv)

c     ------------------------------------------------------------------
c
c     Passed output scalars:
c
c          nrcl   - number of elements in RCL
c
c     ------------------------------------------------------------------
      integer    nrcl

c     ------------------------------------------------------------------
c
c     Passed output array:
c
c          rcl   - Restricted Candidate List
c
c     ------------------------------------------------------------------
      integer    rcl(nv)

c     ------------------------------------------------------------------
c
c     Local scalars:
c
c          coff  - cutoff for RCL
c          gf    - value of greedy function
c          maxgf - maximum value of greedy function
c          mingf - minimum value of greedy function 
c          v     - vertex
c
c     ------------------------------------------------------------------
      integer    gf,     mingf,  maxgf,  coff,   v 

c     ------------------------------------------------------------------
c     Find the maximum value of in-degree(v) * out-degree(v).
c     ------------------------------------------------------------------
      maxgf=0
      mingf=nv*nv

      do 10 v=1,nv
c          -------------------------------------------------------------
c          Consider only existing vertices.
c          -------------------------------------------------------------
           if (vertex(v).eq.1) then
                gf = din(v) *  dout(v)
                maxgf=max(maxgf,  gf)
                mingf=min(mingf,  gf)
           endif
10    continue

c     ------------------------------------------------------------------
c     Compute cut off.
c     ------------------------------------------------------------------
      coff=mingf+ifix(alpha*(maxgf-mingf))

c     ------------------------------------------------------------------
c     Put vertex v into RCL if in-degree(v) * out-degree(v) >= cut off.
c     ------------------------------------------------------------------
      nrcl=0
      do 20 v=1,nv
           if (vertex(v).eq.1) then
                gf=din(v) * dout(v)
                if (gf.ge.coff) then
                     nrcl=nrcl+1
                     rcl(nrcl)=v
                endif
           endif
20    continue
c     ------------------------------------------------------------------
      return

c     ------------------------------------------------------------------
c     End of subroutine mkrcl.
c     ------------------------------------------------------------------
      end





      subroutine outvtx( bcset,  gitrb,  maxv,   ncset,  out,    seedb )
c     ------------------------------------------------------------------
c
c     outvtx: Prints out best feedback vertex set found.
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          gitrb  - iteration best solution was found
c          maxv   - dimension of cutset array
c          ncset  - number of vertices in cutset
c          out    - output device number
c          seedb  - seed at start of GRASP iteration best solution 
c
c     ------------------------------------------------------------------
      integer    gitrb,  maxv,   ncset,  out,    seedb

c     ------------------------------------------------------------------
c
c     Passed input array:
c
c          bcset   - array if vertices of cutset
c
c     ------------------------------------------------------------------
      integer    bcset(maxv)

c     ------------------------------------------------------------------
c
c     Local scalar:
c
c          i      - do loop index
c
c     ------------------------------------------------------------------
      integer    i

c     ------------------------------------------------------------------
c     Write size of minimum cutset found.
c     ------------------------------------------------------------------
      write(out,10)
10    format(//,' GRASP solution---------------------------------',//)
      write(out,20) ncset
20    format('  size of feedback vertex cutset: ',i14,//)
      write(out,30) gitrb
30    format('     iteration best cutset found: ',i14,//)
      write(out,40) seedb
40    format(' seed at start of best iteration: ',i14,//)
c     ------------------------------------------------------------------
c     Write vertices of minimum cutset found.
c     ------------------------------------------------------------------
      write(out,50)
50    format(' smallest vertex cutset: ',/)
      do 70 i=1,ncset
           write(out,60) bcset(i)
60         format(' vertex: ',i8)
70    continue
      write(out,80)
80    format(//,' -----------------------------------------------',//)
      return

c     ------------------------------------------------------------------
c     End of subroutine outvtx.
c     ------------------------------------------------------------------
      end





      subroutine pickv( cset,   ncset,  nrcl,   nv,     rcl,    seed,
     +                  v )
c     ------------------------------------------------------------------
c
c     pickv: Pick a vertex at random from elements of Restricted
c            Candidate List (RCL). 
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          nv     - dimension array
c          seed   - seed for pseudo random number generator
c
c     ------------------------------------------------------------------
      integer    nv,     seed

c     ------------------------------------------------------------------
c
c     Passed input/output scalar:
c
c          ncset  - size of cutset currently found
c          nrcl   - number of vertices currently in RCL
c
c     ------------------------------------------------------------------
      integer    ncset,  nrcl

c     ------------------------------------------------------------------
c
c     Passed input array:
c
c          rcl    - Restricted Candidate List
c
c     ------------------------------------------------------------------
      integer    rcl(nv)

c     ------------------------------------------------------------------ 
c
c     Passed input/output array:
c
c          cset   - cutset currently found 
c
c     ------------------------------------------------------------------
      integer    cset(nv)

c     ------------------------------------------------------------------
c
c     Local scalars and functions:
c
c          nselct - pointer to vertex
c          randp  - portable pseudo-random number generator
c          xrand  - output of randp
c          v      - vertex
c
c     ------------------------------------------------------------------
      integer    nselct, v
      real       xrand,  randp

c     ------------------------------------------------------------------
c     Pick v at random.
c     ------------------------------------------------------------------
      xrand=randp( seed )
      nselct=1+seed/(2147483647/nrcl)
      v=rcl(nselct)

c     ------------------------------------------------------------------
c     Put v in cut set.
c     ------------------------------------------------------------------
      ncset=ncset+1
      cset(ncset)=v

      return
 
c     ------------------------------------------------------------------
c     End of subroutine pickv.
c     ------------------------------------------------------------------
      end





      real function randp( ix )
c     -----------------------------------------------------------------
c
c     randp: Portable pseudo-random number generator.
c            Reference: L. Schrage, "A More Portable Fortran
c            Random Number Generator", ACM Transactions on
c            Mathematical Software, Vol. 2, No. 2, (June, 1979).
c
c     -----------------------------------------------------------------
      integer    a,      p,      ix,     b15,    b16,    xhi,
     +           xalo,   leftlo, fhi,    k
      data       a/16807/, b15/32768/, b16/65536/, p/2147483647/

      xhi=ix/b16
      xalo=(ix-xhi*b16)*a
      leftlo=xalo/b16
      fhi=xhi*a+leftlo
      k=fhi/b15
      ix=(((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
      if (ix.lt.0) ix=ix+p

      randp=float(ix)*4.656612875e-10

      return

c     ------------------------------------------------------------------
c     End of subroutine randp.
c     ------------------------------------------------------------------
      end





      subroutine readp( errcnd, in,     ne,     nv,     maxe,   maxv,
     +                  vtx1,   vtx2 )
c     ------------------------------------------------------------------
c
c     readp: Inputs problem data.
c
c     ------------------------------------------------------------------
c
c     Passed input parameters:
c
c          in     - input device number
c          maxe   - array dimension
c          maxv   - array dimension
c
c     ------------------------------------------------------------------
      integer    in,     maxe,   maxv

c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          nv     - number of nodes
c          ne     - number of arcs
c
c     ------------------------------------------------------------------
      integer    ne,     nv

c     ------------------------------------------------------------------
c     
c     Passed output scalars:
c
c          errcnd - error condition
c
c     ------------------------------------------------------------------
      integer    errcnd

c     ------------------------------------------------------------------
c
c     Passed output arrays:
c
c          vtx1  - node 1 of arc
c          vtx2  - node 2 of arc
c
c     ------------------------------------------------------------------
      integer    vtx1(maxe),     vtx2(maxe)

c     ------------------------------------------------------------------
c
c     Local parameters:
c
c          i     - do loop parameter
c
c     ------------------------------------------------------------------
      integer    i

c     ------------------------------------------------------------------
c     Initilaize error condition.
c     ------------------------------------------------------------------
      errcnd=0

c     ------------------------------------------------------------------
c     Read graph dimensions.
c     ------------------------------------------------------------------
      read (in,*) nv,ne

c     ------------------------------------------------------------------
c     Check if the number of vertices of the input graph exceeds maxv
c     ------------------------------------------------------------------
      if (nv .gt.maxv) then
           errcnd=1
           return
      endif

c     ------------------------------------------------------------------
c     Check if the number of arcs of the input graph exceeds maxe
c     ------------------------------------------------------------------
      if (ne .gt.maxe) then
           errcnd=2
           return
      endif

c     ------------------------------------------------------------------
c     Read arcs.
c     ------------------------------------------------------------------
      do 10 i=1,ne
c          -------------------------------------------------------------
c          Read directed edge = (v1,v2).
c          -------------------------------------------------------------
           read(in,*) vtx1(i),vtx2(i)
10    continue
      return
 
c     ------------------------------------------------------------------
c     End of subroutine readp.
c     ------------------------------------------------------------------
      end





      subroutine redg( cset,   din,    dout,   lin,    lout,   ncset,
     +                 n2e,    nin,    nout,   nv,     nvleft, pin,
     +                 pout,   vertex, vin,    vout )
c     ------------------------------------------------------------------

c     redg:  Apply reductions to graph.  Keep reapplying until no
c            further reduction is possible.
c
c            Adapt greedy function by changing degrees.

c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          n2e    - array dimension = 2*ne
c          ncset  - size of cutset found
c          nv     - number of nodes in graph
c          nvleft - number of unselected vertices
c          red    - indicator of reduction
c
c     ------------------------------------------------------------------
      integer    n2e,    ncset,  nv,     nvleft, red

c     ------------------------------------------------------------------
c
c     Passed input/output arrays:
c
c         cset    - cutset currently found 
c         din     - in-degree of vertex (copy)
c         dout    - out-degree of vertex (copy)
c         lin     - pointer to end of in-adjacency list (copy)
c         lout    - pointer to end of out-adjacency list (copy)
c         nin     - pointer to next element of in-adjacency
c                   list (copy)
c         nout    - pointer to next element of out-adjacency
c                   list (copy)
c         pin     - pointer to start of in-adjacency list (copy)
c         pout    - pointer to start of out-adjacency list (copy)
c         vertex  - existence of vertex in the reduced graph 
c         vin     - vertex in in-adjacency list (copy)
c         vout    - vertex in out-adjacency list (copy)
c
c     ------------------------------------------------------------------
      integer    cset(nv),       din(nv),        dout(nv),
     +           lin(n2e),       lout(n2e),      nin(n2e),
     +           nout(n2e),      pin(nv),        pout(nv),
     +           vertex(nv),     vin(n2e),       vout(n2e)

c     ------------------------------------------------------------------
c     Begin reduction loop.
c     ------------------------------------------------------------------
10    continue
           red=0
c          -------------------------------------------------------------
c          Apply reduction type 3:
c
c          If, for any vertex v, (v,v) is in E, then
c
c               v in CUTSET and V = V - {v}, 
c
c               E = E - OUT_EDGES(v) - IN_EDGES(v)
c          -------------------------------------------------------------
           call rgloop( cset,   din,    dout,   lin,    lout,   ncset,
     +                  n2e,    nin,    nout,   nv,     nvleft, pin,
     +                  pout,   red,    vertex, vin,    vout )

c          -------------------------------------------------------------
c          If reduction was made, restart.
c          -------------------------------------------------------------
           if (red .eq. 1) goto 10

c          -------------------------------------------------------------
c          Apply graph reduction type 1:
c
c          If, for any vertex v, out_degree(v) = 0 or 
c                                 in_degree(v) = 0, then
c
c          V = V - {v} and E = E - IN_EDGES(v) - OUT_EDGES(v)
c          -------------------------------------------------------------
           call rgdeg0( din,    dout,   lin,    lout,   n2e,    nin,
     +                  nout,   nv,     nvleft, pin,    pout,   red,
     +                  vertex, vin,    vout )

c          -------------------------------------------------------------
c          If reduction was made, restart.
c          -------------------------------------------------------------
           if (red .eq. 1) goto 10

c          -------------------------------------------------------------
c          Apply reduction type 2:
c
c          If, for any vertex v, in_degree(v) = 1, and (u,v) in E, then
c
c          V = V - {v}, 
c          E = E - IN_EDGES(v) - OUT_EDGES(v) + (u, OUT_VERTICES(v))
c          -------------------------------------------------------------
           call rgd1( din,    dout,   lin,    lout,   n2e,    nin,
     +                nout,   nv,     nvleft, pin,    pout,   red,
     +                vertex, vin,    vout )

c          -------------------------------------------------------------
c          If reduction was made, restart.
c          -------------------------------------------------------------
           if (red .eq. 1) goto 10

c          -------------------------------------------------------------
c          Apply reduction type 4:
c
c          If, for any vertex v, out_degree(v) = 1, and (v,u) in E, then
c
c          V = V - {v}, 
c          E = E - OUT_EDGES(v) - IN_EDGES(v) + (IN_VERTICES(v),u)
c          -------------------------------------------------------------
           call rgd1( dout,   din,    lout,   lin,    n2e,    nout,
     +                nin,    nv,     nvleft, pout,   pin,    red,
     +                vertex, vout,   vin )

c          -------------------------------------------------------------
c          If reduction was made, restart.
c          -------------------------------------------------------------
           if (red .eq. 1) goto 10

c     ------------------------------------------------------------------
c     If no more reductions are possible, return.
c     ------------------------------------------------------------------
      return
c
c     ------------------------------------------------------------------
c     End of subroutine redg.
c     ------------------------------------------------------------------
      end





      subroutine rgd1( dxx,    dxxx,   lxx,    lxxx,   n2e,    nxx,
     +                 nxxx,   nv,     nvleft, pxx,    pxxx,   red,
     +                 vertex, vxx,    vxxx )
c     ------------------------------------------------------------------
c
c     rgd1: Apply reduction type 2 and type 3:
c
c     If, for any vertex v, in_degree(v) = 1, and (u,v) in E, then
c
c     v not_in CUTSET and V = V - {v}, 
c
c     E = E - IN_EDGES(v) - OUT_EDGES(v) + (u, OUT_VERTICES(v))
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          n2e    - array dimension = 2*ne
c          nv     - number of nodes in graph
c          nvleft - number of unselected vertices
c          red    - indicator of reduction
c
c     ------------------------------------------------------------------
      integer    n2e,    nv,     nvleft, red

c     ------------------------------------------------------------------
c
c     Passed input/output arrays:
c
c         dxx     - in-degree of vertex (copy)
c         dxxx    - out-degree of vertex (copy)
c         lxx     - pointer to end of in-adjacency list (copy)
c         lxxx    - pointer to end of out-adjacency list (copy)
c         nxx     - pointer to next element of in-adjacency
c                   list (copy)
c         nxxx    - pointer to next element of out-adjacency
c                   list (copy)
c         pxx     - pointer to start of in-adjacency list (copy)
c         pxxx    - pointer to start of out-adjacency list (copy)
c         vertex  - existence of vertex in the reduced graph
c         vxx     - vertex in in-adjacency list (copy)
c         vxxx    - vertex in out-adjacency list (copy)
c
c     ------------------------------------------------------------------
      integer    dxx(nv),        dxxx(nv),       lxx(n2e),  
     +           lxxx(n2e),      nxx(n2e),       nxxx(n2e),  
     +           pxx(nv),        pxxx(nv),       vertex(nv),  
     +           vxx(n2e),       vxxx(n2e)

c     ------------------------------------------------------------------
c
c     Local scalars:
c
c         me     - existence of arc in the reduced graph
c         p      - pointer to arc
c         r      - pointer to arc
c         u      - pointer to vertex
c         v      - vertex index
c         x      - pointer to vertex
c         y      - pointer to arc
c         w      - pointer to arc
c
c     ------------------------------------------------------------------
      integer    me,     p,      r,      u,      v,      x,
     +           y,      w

c     ------------------------------------------------------------------
      do 90 v=1,nv
           if (vertex(v).ne.0) then
                if (dxx(v).eq.1) then
c                    ---------------------------------------------------
c                    If v has a self loop, then skip rgdin.
c                    ---------------------------------------------------
                     if(vxx(pxx(v)).eq.v) then
                          red=1
                          goto 90
                     endif    

c                    ---------------------------------------------------
c                    u is the vertex such that (u,v) or (v,u) in E.
c                    ---------------------------------------------------
                     u=vxx(pxx(v))

c                    ---------------------------------------------------
c                    Delete v from the xxx_adj of u.
c                    ---------------------------------------------------
                     y=pxxx(u)
10                   if (y.eq.0) goto 20 
                          if (vxxx(y).eq.v) then
                               call dele( lxxx,   n2e,    nv,     nxxx,
     +                                    pxxx,   y,      u )
                               dxxx(u)=dxxx(u)-1
                               goto 20
                          else
                               y=nxxx(y)
                               goto 10
                          endif   
20                   continue      
c                    ---------------------------------------------------
c                    Start to scan the xxx_adj of v.
c                    ---------------------------------------------------
                     p=pxxx(v)

30                   if(p.eq.0) goto 80
c                    ---------------------------------------------------
c                    Temporarily save nxxx(p) in r.
c                    ---------------------------------------------------
                     r=nxxx(p)
                     x=vxxx(p)
c                    ---------------------------------------------------
c                    If vxxx(p)=u , then  u has a self loop.
c                    ---------------------------------------------------
                     if(x.eq.u) then
c                         ----------------------------------------------
c                         Delete edge p from the xxx_adj of v.
c                         ----------------------------------------------
                          call dele( lxxx,   n2e,    nv,     nxxx,
     +                               pxxx,   p,      v )
                          vxx(p)=u

c                         ----------------------------------------------
c                         Add edge p to the xxx_adj of u.
c                         ----------------------------------------------
                          call adde( p,      lxxx,   n2e,    nv,  
     +                               nxxx,   pxxx,   u,      vxxx,
     +                               u )
                          dxxx(u)=dxxx(u)+1
                     else 
c                         ----------------------------------------------
c                         Otherwise, check whether edge (u,x) or (x,u) 
c                         exists, and set me =1 if so.
c                         ----------------------------------------------
                          y=pxxx(u)
                          me=0
40                        if (y.eq.0) goto 50
                          if (vxxx(y).eq.x) then
                               me=1
                               goto 50
                          else
                               y=nxxx(y)
                               goto 40
                          endif
50                        continue
c                         ----------------------------------------------
c                         If (u,x) or (x,u) exists, delete edge p from 
c                         xx_adj of x.
c                         ----------------------------------------------
                          if (me.eq.1) then
                               call dele( lxx,    n2e,    nv,     nxx,
     +                                    pxx,    p,      x )
                               dxx(x)=dxx(x)-1
                          else
c                              -----------------------------------------
c                              Otherwise, add (u,x) or (x,u) to the 
c                              xxx_adj of x.
c                              -----------------------------------------
                               w=pxx(x)
60                             if(w.eq.0) goto 70
                               if (vxx(w).eq.v) then
                                    vxx(w)=u
                                    call dele( lxxx,   n2e,    nv,   
     +                                         nxxx,   pxxx,   w, 
     +                                         v )
                                    dxxx(v)=dxxx(v)-1
                                    call adde( w,      lxxx,   n2e,  
     +                                         nv,     nxxx,   pxxx, 
     +                                         u,      vxxx,   x )
                                    dxxx(u)=dxxx(u)+1
                               else
                                    w=nxx(w)
                                    goto 60
                               endif    
70                             continue
                          endif
                     endif   

c                    ---------------------------------------------------
c                    Let p=nxxx(p).
c                    ---------------------------------------------------
                     p=r
                     goto 30
80                   continue 
                     red=1
                     vertex(v)=0
                     nvleft=nvleft-1
                endif
           endif
90    continue
      return
 
c     ------------------------------------------------------------------
c     End of subroutine rdg1.
c     ------------------------------------------------------------------
      end





      subroutine rgdeg0( din,    dout,   lin,    lout,   n2e,    nin,
     +                   nout,   nv,     nvleft, pin,    pout,   red,
     +                   vertex, vin,    vout )
c     ------------------------------------------------------------------
c
c     rgdeg0: Apply graph reduction type 1:
c
c          If, for any vertex v, out_degree(v) = 0 or in_degree(v) = 0,
c
c          then 
c             v not_in CUTSET and 
c
c                 V = V - {v} 
c                 E = E - INCIDENT_EDGES(v)
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          n2e    - array dimension = 2*ne
c          nv     - number of nodes in graph
c          nvleft - number of unselected vertices
c          red    - indicator of reduction
c
c     ------------------------------------------------------------------
      integer    n2e,    nv,     nvleft, red

c     ------------------------------------------------------------------
c
c     Passed input/output arrays:
c
c         din     - in-degree of vertex (copy)
c         dout    - out-degree of vertex (copy)
c         lin     - pointer to end of in-adjacency list (copy)
c         lout    - pointer to end of out-adjacency list (copy)
c         nin     - pointer to next element of in-adjacency
c                   list (copy)
c         nout    - pointer to next element of out-adjacency
c                   list (copy)
c         pin     - pointer to start of in-adjacency list (copy)
c         pout    - pointer to start of out-adjacency list (copy)
c         vertex  - existence of vertex in the reduced graph
c         vin     - vertex in in-adjacency list (copy)
c         vout    - vertex in out-adjacency list (copy)
c
c     ------------------------------------------------------------------
      integer    din(nv),        dout(nv),       lin(n2e),  
     +           lout(n2e),      nin(n2e),       nout(n2e),  
     +           pin(nv),        pout(nv),       vertex(nv),   
     +           vin(n2e),       vout(n2e)

c     ------------------------------------------------------------------
c
c     Local scalars:
c
c         dinv   - in-degree of vertex
c         doutv  - out-degree of vertex
c         v      - vertex index
c
c     ------------------------------------------------------------------
      integer    dinv,   doutv,  v

c     ------------------------------------------------------------------
c     Scan all vertices of GRAPH.
c     ------------------------------------------------------------------
      do 10 v=1,nv
c          -------------------------------------------------------------
c          If vertex v is in GRAPH, check if reduction is possible.
c          -------------------------------------------------------------
           if (vertex(v).ne.0) then
c               --------------------------------------------------------
c               If in-degree or out-degree of v is 0, apply reduction.
c               --------------------------------------------------------
                dinv=din(v)
                doutv=dout(v)
                if (dinv.eq.0 .or. doutv.eq.0) then
                     red=1
c                    ---------------------------------------------------
c                    Delete v and all edges incident to v from GRAPH.
c                    ---------------------------------------------------
                     if ((dinv+doutv) .eq. 0) then
                          vertex(v)=0
                          nvleft=nvleft-1
                     else
                          call delv( din,    dout,   lin,    lout,
     +                               n2e,    nin,    nout,   nv,  
     +                               nvleft, pin,    pout,   v,
     +                               vertex, vin,    vout )
                     endif
                endif
           endif
10    continue
c     ------------------------------------------------------------------
      return

c     ------------------------------------------------------------------
c     End of subroutine rgdeg0.
c     ------------------------------------------------------------------
      end





      subroutine rgloop( cset,   din,    dout,   lin,    lout,   ncset,  
     +                   n2e,    nin,    nout,   nv,     nvleft, pin,
     +                   pout,   red,    vertex, vin,    vout )
c     ------------------------------------------------------------------
c
c     rgloop: Apply reduction type 3:
c
c          If, for any vertex v, (v,v) is in E, then
c
c          v in CUTSET and V = V - {v}, 
c
c          E = E - OUT_EDGES(v) - IN_EDGES(v)
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          n2e    - array dimension = 2*ne
c          ncset  - size of cutset found
c          nv     - number of nodes in graph
c          nvleft - number of unselected vertices
c          red    - reduction indicator
c
c     ------------------------------------------------------------------
      integer    n2e,    ncset,  nv,     nvleft, red

c     ------------------------------------------------------------------
c
c     Passed input/output arrays:
c
c         cset    - cutset currently found
c         din     - in-degree of vertex (copy)
c         dout    - out-degree of vertex (copy)
c         lin     - pointer to end of in-adjacency list (copy)
c         lout    - pointer to end of out-adjacency list (copy)
c         nin     - pointer to next element of in-adjacency
c                   list (copy)
c         nout    - pointer to next element of out-adjacency
c                   list (copy)
c         pin     - pointer to start of in-adjacency list (copy)
c         pout    - pointer to start of out-adjacency list (copy)
c         vertex  - existence of vertex in the reduced graph
c         vin     - vertex in in-adjacency list (copy)
c         vout    - vertex in out-adjacency list (copy)
c
c     ------------------------------------------------------------------
      integer    cset(nv),       din(nv),        dout(nv),
     +           lin(n2e),       lout(n2e),      nin(n2e),  
     +           nout(n2e),      pin(nv),        pout(nv),
     +           vertex(nv),     vin(n2e),       vout(n2e)

c     ------------------------------------------------------------------
c
c     Local scalars:
c
c         p      - pointer to a vertex
c         v      - vertex index
c
c     ------------------------------------------------------------------
      integer    p,      v

c     ------------------------------------------------------------------
c     Scan all vertices.
c     ------------------------------------------------------------------
      do 30 v=1,nv
c          -------------------------------------------------------------
c          If vertex v is still in graph, try to apply reduction to it.
c          -------------------------------------------------------------
           if (vertex(v).ne.0) then
c               --------------------------------------------------------
c               Scan all edges into v.
c               --------------------------------------------------------
                p=pin(v)

c               --------------------------------------------------------
c               Check if there are still more edges incident to v.
c               --------------------------------------------------------
10              if (p.eq.0) goto 20
c                    ---------------------------------------------------
c                    If there still edges incident to v,
c                    check if the edge is (v,v).
c                    ---------------------------------------------------
                     if (vin(p).eq.v) then
c                         ----------------------------------------------
c                         If edge is (v,v), add v to CUTSET.
c                         ----------------------------------------------
                          ncset=ncset+1
                          cset(ncset)=v
                          red=1

c                         ----------------------------------------------
c                         Delete v from GRAPH, along with all edges
c                         incident to it.
c                         ----------------------------------------------
                          call delv( din,    dout,   lin,    lout,
     +                               n2e,    nin,    nout,   nv,   
     +                               nvleft, pin,    pout,   v,
     +                               vertex, vin,    vout )

                          goto 30
                     endif
c                    ---------------------------------------------------
c                    Get next edge incident to v.
c                    ---------------------------------------------------
                     p=nin(p)
                     goto 10
20              continue
           endif
30    continue
c     ------------------------------------------------------------------
      return
 
c     ------------------------------------------------------------------
c     End of subroutine rgloop.
c     ------------------------------------------------------------------
      end





      subroutine warmup( seed )
c     ------------------------------------------------------------------
c
c     warmup: Warms up random number generator.
c
c     ------------------------------------------------------------------
c
c     Passed input/output scalars:
c
c          seed   - seed for pseudo random number generator
c
c     ------------------------------------------------------------------
      integer    seed
c     ------------------------------------------------------------------
c
c     Local parameter:
c
c          wapitr - max number of warmup iterations
c
c     ------------------------------------------------------------------
      integer    wupitr
      parameter  (wupitr=100)

c     ------------------------------------------------------------------
c
c     Local scalars and functions:
c
c          i      - do loop index 
c          ix     - output of the function randp
c          randp  - portable pseudo-random number generator
c
c     ------------------------------------------------------------------
      integer    i
      real       ix,     randp

c     ------------------------------------------------------------------
c     Call random number generator wupitr times to warm it up.
c     ------------------------------------------------------------------
      do 10 i=1,wupitr
           ix=randp( seed )
10    continue
      return
c
c     ------------------------------------------------------------------
c     End of subroutine warmup.
c     ------------------------------------------------------------------
      end

      program dgfas
c     ------------------------------------------------------------------
c
c     dgfas:         Driver for gfas, a GRASP for finding approximate
c                    solutions to the Feedback Arc Set Problem
c
c
c     authors:       Paola Festa             [paofes@udsab.dia.unisa.it]
c                    Mauricio G.C. Resende   [mgcr@research.att.com]
c                    Panos P.M. Pardalos     [pardalos@ufl.edu]
c
c     ------------------------------------------------------------------
c
c     Dimensions of maximum graph size:
c
c          maxe0 - max number of arcs (original graph)
c          maxe  - max number of arcs (tranformed graph)
c          maxv  - max number of vertices (tranformed graph)
c
c     ------------------------------------------------------------------
      integer    maxe0
      parameter  (maxe0=1500)
      integer    maxe
      parameter  (maxe=maxe0*maxe0)
      integer    maxv
      parameter  (maxv=maxe0)

c     ------------------------------------------------------------------
c
c     Other array dimension.
c
c          max2e  - max number of vertices
c
c     ------------------------------------------------------------------
      integer    max2e
      parameter  (max2e=2*maxe)
c     ------------------------------------------------------------------
c
c     Input/output device numbers passed to errmsg, gfvas, outarc, and
c     readp:
c
c          in     - input device
c          out    - output device
c
c     ------------------------------------------------------------------
      integer    in,     out
      parameter  (in=5,out=6)

c     ------------------------------------------------------------------
c
c     Scalars passed to gfas:
c
c          alpha  - GRASP RCL parameter
c          bestz  - size of best cutset found
c          errcnd - error condition
c          fin0   - number of incoming arcs (original)
c          fout0  - number of outgoing arcs (original)
c          gitrb  - iteration cut best set was found
c          look4  - size of cutset sought
c          maxitr - maximum number of GRASP iterations
c          ne     - number of arcs in original graph
c          nrcl   - number of vertices in RCL
c          nv     - number of nodes in original graph
c          out    - Fortran output device
c          prttyp - iteration summary report
c                   = 0 - silent
c                   = 1 - show improvements only
c                   = 2 - show all iterations
c          seedb  - seed at start of best iteration
c          seedd  - seed for pseudo random number generator
c
c     ------------------------------------------------------------------
      real             alpha
      integer          bestz,  errcnd, fin0,   fout0,  gitrb,  look4,
     +                 maxitr, ne,     nrcl,   nv,     prttyp, seedb
      double precision seedd

c     ------------------------------------------------------------------
c
c     Arrays passed to gfas:
c          
c           adjv    - out-adjacency list
c           bcset   - best cutset found
c           caset1  - node 1 of arc in cutset found for the feedback arc
c                     problem
c           caset2  - node 2 of arc in cutset found for the feedback arc
c                     problem
c           cset    - cutset currently found (copy)
c           cset0   - cutset currently found
c           din     - in-degree of vertex (copy)
c           din0    - in-degree of vertex
c           dout    - out-degree of vertex (copy)
c           dout0   - out-degree of vertex
c           first   - pointer to start of out-adjacency list
c           ind     - index of vertex in the translated problem
c           last    - pointer to end of out-adjacency list
c           lin     - pointer to end of in-adjacency list (copy)
c           lin0    - pointer to end of in-adjacency list (original)
c           lout    - pointer to end of out-adjacency list (copy)
c           lout0   - pointer to end of out-adjacency list (original)
c           next    - pointer to next element of out-adjacency
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
c           vertex  - existence of vertex in the graph
c           vin     - vertex in in-adjacency list (copy)
c           vin0    - vertex in in-adjacency list (original)
c           vmap1   - node 1 of arc in the translated problem
c           vmap2   - node 2 of arc in the translated problem
c           vout    - vertex in out-adjacency list (copy)
c           vout0   - vertex in out-adjacency list (original)
c           vtx1    - vertex 1 of arc
c           vtx2    - vertex 2 of arc
c          
c     ------------------------------------------------------------------
      integer    adjv(max2e),      bcset(maxe),      caset1(maxe),
     +           caset2(maxe),     cset(maxv),       cset0(maxe),
     +           din(maxv),        din0(maxv),       dout(maxv),
     +           dout0(maxv),      first(maxv),      ind(maxv*maxv), 
     +           last(max2e),      lin(max2e),       lin0(max2e),
     +           lout(max2e),      lout0(max2e),     next(max2e),
     +           nin(max2e),       nin0(max2e),      nout(max2e),
     +           nout0(max2e),     pin(maxv),        pin0(maxv),
     +           pout(maxv),       pout0(maxv),      rcl(maxv),
     +           vertex(maxv),     vin(max2e),       vin0(max2e),
     +           vmap1(maxv*maxv), vmap2(maxv*maxv), vout(max2e),
     +           vout0(max2e),     vtx1(maxe),       vtx2(maxe)

c     ------------------------------------------------------------------
c     Sets default values for parameters.
c     ------------------------------------------------------------------
      call dflset( alpha,  look4,  maxitr, prttyp, seedd )

c     ------------------------------------------------------------------
c     Read input graph.
c     ------------------------------------------------------------------
      call readp( errcnd, in,     ne,     nv,     maxe,   maxv,   vtx1,
     +            vtx2 )

c     ------------------------------------------------------------------
c     Solve problem approximately using GRASP.
c     ------------------------------------------------------------------
      if (errcnd.eq.0) then
           call gfas( adjv,   alpha,  bcset,  bestz,  caset1, caset2,
     +                cset,   cset0,  din,    din0,   dout,   dout0,
     +                errcnd, fin0,   first,  fout0,  gitrb,  ind,
     +                last,   lin,    lin0,   look4,  lout,   lout0,
     +                max2e,  maxe,   maxitr, maxv,   ne,     next,
     +                nin,    nin0,   nout,   nout0,  nrcl,   nv,
     +                pin,    pin0,   pout,   pout0,  prttyp, out,
     +                rcl,    seedb,  seedd,  vertex, vin,    vin0,
     +                vmap1,  vmap2,  vout,   vout0,  vtx1,   vtx2 )

           call errmsg( errcnd, out )
c          -------------------------------------------------------------
c          Print the cut set.
c          -------------------------------------------------------------
           if (errcnd .eq. 0) then
                call outarc( caset1, caset2, gitrb,  maxe,   bestz,
     +                       out,    seedb )
           endif
      else
           call errmsg( errcnd, out )
      endif
      stop

c     ------------------------------------------------------------------
c     End driver for gfas.
c     ------------------------------------------------------------------
      end

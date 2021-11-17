      program dgfvs
c     ------------------------------------------------------------------
c
c     dggvs:         Driver for gfvs, a GRASP for finding approximate
c                    solutions to the Feedback Vertex Set Problem 
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
c          maxe  - max number of arcs 
c          maxv  - max number of vertices
c
c     ------------------------------------------------------------------
      integer    maxe,   maxv
      parameter  (maxe=1500*1500,maxv=1500)

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
c     Input/output device numbers passed to errmsg, gfvs, outvtx, and
c     readp: 
c
c          in     - Fortran input device
c          out    - Fortran output device
c
c     ------------------------------------------------------------------
      integer   in,      out
      parameter (in=5,out=6)
c     ------------------------------------------------------------------
c
c     Scalars passed to gfvs:
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
c     Arrays passed to gfvs:
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
c           vertex  - existence of vertex in the graph
c           vin     - vertex in in-adjacency list (copy)
c           vin0    - vertex in in-adjacency list (original)
c           vout    - vertex in out-adjacency list (copy)
c           vout0   - vertex in out-adjacency list (original)
c           vtx1    - vertex 1 of arc
c           vtx2    - vertex 2 of arc

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
           call gfvs( alpha,  bcset,  bestz,  cset,   cset0,  din,
     +                din0,   dout,   dout0,  errcnd, fin0,   fout0,
     +                gitrb,  lin,    lin0,   look4,  lout,   lout0,
     +                maxe,   maxitr, maxv,   max2e,  ne,     nin,
     +                nin0,   nout,   nout0,  nrcl,   nv,     out,
     +                pin,    pin0,   pout,   pout0,  prttyp, rcl,
     +                seedb,  seedd,  vertex, vin,    vin0,   vout,
     +                vout0,  vtx1,   vtx2 )

           call errmsg( errcnd, out )
c          -------------------------------------------------------------
c          Print the cut set.
c          -------------------------------------------------------------
           if (errcnd .eq. 0) then
                call outvtx( bcset,  gitrb,  maxv,   bestz,  out,    
     +                       seedb )
           endif
      else
         call errmsg( errcnd, out )
      endif
      stop

c     ------------------------------------------------------------------
c     End driver for gfvs.
c     ------------------------------------------------------------------
      end

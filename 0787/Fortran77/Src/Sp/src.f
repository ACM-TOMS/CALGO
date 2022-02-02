      subroutine gmis(adj,adjn,adjn0,alpha,beta,biter,bkset,btup,deg,
     +                deg0,errcnd,fdmiq,fdmq,fwdmap,incdij,iniad,
     +                iniad0,invmap,iq,ittype,kset,kset0,look4,lstype,
     +                lswhen,max2a,max2n2,maxa,maxitr,maxn,maxn2,maxt,
     +                maxtd,narcs0,nkset0,nlowdg,nnode0,node1,node2,out,
     +                ptrn,ptrn0,q,rcl,seed,tupdim,tupind,tuple,tuplmt,
     +                tuplst,vtup)
c     ------------------------------------------------------------------
c
c     gmis: A package of FORTRAN subroutines for finding large 
c           independent sets in a graph using GRASP.
c
c     authors: Mauricio G.C. Resende [mgcr@research.att.com]
c              Thomas A. Feo         [feo@emx.utexas.edu]
c              Stuart H. Smith       [stuart@alk.com]
c
c     The algorithm implemented in this package is based on the paper:
c     Thomas A. Feo, Mauricio G.C. Resende, and S.H. Smith, "A greedy 
c     randomized adaptive search procedure for maximum independent 
c     set," Operations Research, vol. 42, pp. 860--878, 1994.
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          alpha  - GRASP RCL parameter
c          beta   - tuple list size restriction parameter 
c          ittype - iteration summary report
c                   = 0 - silent
c                   = 1 - show improvements only
c                   = 2 - show all iterations
c          look4  - size of independent set sought
c          lstype - type of local search
c          lswhen - indicates when local search is done
c          max2a  - array dimension = 2*maxa
c          max2n2 - array dimension = 2*maxn2
c          maxa   - array dimension
c          maxitr - maximum number of GRASP iterations
c          maxn   - array dimension
c          maxn2  - array dimension
c          maxt   - maximum number of tuples in freedom queue
c          maxtd  - array dimension
c          narcs0 - number of arcs in original graph
c          nlowdg - number of low degree vertices used to make tuples 
c          nnode0 - number of nodes in original graph
c          out    - Fortran output device
c          tupdim - dimension of tuple vector
c          tuplmt - maximum number of tuples allowed
c
c     ------------------------------------------------------------------
      real     alpha,beta
      integer  ittype,look4,lstype,lswhen,max2a,max2n2,maxa,maxitr,
     +         maxn,maxn2,maxt,maxtd,narcs0,nlowdg,nnode0,out,tupdim,
     +         tuplmt
c     ------------------------------------------------------------------
c
c     Passed input/output scalars:
c
c          errcnd - error condition
c          seed   - seed for pseudo random number generator
c
c     ------------------------------------------------------------------
      integer  errcnd,seed
c     ------------------------------------------------------------------
c
c     Passed working scalars:
c
c          gitr   - number of grasp iterations taken
c          narcs  - number of arcs in reduced graph
c          nkset  - number of vertices in independent set of reduced
c                   graph
c          nnode  - number of vertices in reduced graph
c          ntuple - number of tuples produced by subroutine mktup
c
c     ------------------------------------------------------------------
      integer gitr,narcs,nkset,nnode,ntuple
c     ------------------------------------------------------------------
c
c     Passed output scalars:
c
c          biter  - number of iterations to find best set
c          btup   - number of tuples to find best set
c          nkset0 - number of vertices in independent set
c
c     ------------------------------------------------------------------
      integer biter,btup,nkset0
c     ------------------------------------------------------------------
c
c     Passed input arrays:
c
c          node1  - vertex 1 of egde
c          node2  - vertex 2 of edge
c
c     ------------------------------------------------------------------
      integer   node1(maxa),node2(maxa)
c     ------------------------------------------------------------------
c
c     Passed working arrays:
c
c          adj    - indicator array of adjacent nodes
c          adjn   - vertex in adjacency list (copy)
c          adjn0  - vertex in adjacency list (original)
c          bkset  - set of vertices in independent set
c          deg    - degree of vertex (copy)
c          deg0   - degree of vertex (original)
c          fdmiq  - freedom heap indices 
c          fdmq   - freedom heap values 
c          fwdmap - node i in original graph is fwdmap(i) in reduced
c                   graph
c          incdij - graph incidence matrix
c          iniad  - pointer to start of adjacency list (copy)
c          iniad0 - pointer to start of adjacency list (original)
c          invmap - node i in reduced graph is invmap(i) in original 
c                   graph, passed to mkrgrf, grspss
c          iq     - heap index array
c          kset   - set of independent nodes
c          ptrn   - pointer to next element of adjacency list (copy) 
c          ptrn0  - pointer to next element of adjacency list (original)
c          q      - heap value array
c          rcl    - restricted candidate list 
c          tupind - index of tuple 
c          tuple  - list of tuples
c          tuplst - list of tuples
c          vtup   - tuple vertices 
c
c     ------------------------------------------------------------------
      integer   adj(maxn2),adjn(max2a),adjn0(max2a),bkset(maxn),
     +          deg(maxn),deg0(maxn),fdmiq(maxt),fdmq(maxt),
     +          fwdmap(maxn),incdij(maxn2),iniad(maxn),iniad0(maxn),
     +          invmap(maxn),iq(maxn2),kset(maxn),ptrn(max2a),
     +          ptrn0(max2a),q(maxn2),rcl(maxn),tupind(maxtd),
     +          tuple(max2n2),tuplst(max2n2),vtup(maxtd)
c     ------------------------------------------------------------------
c
c     Passed output array:
c
c          kset0  - best set of independent nodes, in original graph
c
c     ------------------------------------------------------------------
      integer   kset0(maxn)
c     ------------------------------------------------------------------
c
c     Local scalars:
c
c          biter0 - number of iterations in grspss to find best set
c          i      - do loop index
c          k      - do loop index
c          maxtup - maximum number of tuples to be scanned
c          nq     - number of elements in heap
c          tmtm1  - tupdim*(tpl-1)
c          tpl    - do loop index (tuple number)
c
c     ------------------------------------------------------------------
      integer  biter0,i,k,maxtup,nq,tmtm1,tpl

c     ------------------------------------------------------------------
c     Check if parameters are set correctly.
c     ------------------------------------------------------------------
      call chkerr(alpha,beta,errcnd,look4,lstype,lswhen,maxitr,
     +             maxtd,nlowdg,nnode0,seed,tupdim)
      if (errcnd .ge. 100) return
      if (ittype.ge.1) then
           write(out,10)
10         format(
     +' iteration summary---------------------------------------------')
      endif

c     ------------------------------------------------------------------
c     Create adjacency data structure.
c     ------------------------------------------------------------------
      call mkds(adjn0,deg0,incdij,iniad0,max2a,maxa,maxn,maxn2,narcs0,
     +          nnode0,node1,node2,ptrn0)
c     ------------------------------------------------------------------
c     Create list of tuples for processing reduced graph
c 
c            G - VTUP - adjv(VTUP). 
c
c     Array invmap is passed as a work array named lowdeg in mktup. 
c     ------------------------------------------------------------------
      if (tupdim .ne. 0) then
           call mktup(adj,adjn0,beta,deg0,errcnd,fdmiq,fdmq,incdij,
     +          iniad0,iq,invmap,max2a,max2n2,maxn,maxn2,maxt,maxtd,
     +          nlowdg,nnode0,ntuple,ptrn0,q,tupdim,tupind,tuple,
     +          tuplst,vtup)
      endif

c     ------------------------------------------------------------------
c     Initialize nkset: number of vertices in max indep set.
c
c     Limit number of tuples to be examined to at most tuplmt.
c     ------------------------------------------------------------------
      nkset0=0
      if (tupdim .gt.0) then
           maxtup=min0(tuplmt,ntuple)
      else
           maxtup=1
      endif
c     ------------------------------------------------------------------
c     For each tuple (vtup), run GRASP on reduced graph
c     G - vtup - adjv(vtup).
c     ------------------------------------------------------------------
      maxitr=maxitr/maxtup
      do 50 tpl=1,maxtup

c          -------------------------------------------------------------
c          Identify vertices (vtup) of tuple.
c          ------------------------------------------------------------
           tmtm1=tupdim*(tpl-1)
           do 20 k=1,tupdim
                vtup(k)=tuplst(tmtm1+k)
20         continue

c          -------------------------------------------------------------
c          Make the reduced graph G - vtup - adjv(vtup)
c
c          The graph is placed in data structure (nnode,narcs,iniad,
c          ptrn,adjn,deg).
c
c          Node i in reduced  graph is invmap(i) in original graph. 
c          Node i in original graph is fwdmap(i) in reduced  graph. 
c
c          Array tuple is passed as a work array named vthere in mkrgrf.
c          -------------------------------------------------------------
           call mkrgrf(adjn,adjn0,deg,fwdmap,iniad,iniad0,invmap,
     +                 max2a,max2n2,maxn,maxtd,narcs,nnode,nnode0,ptrn,
     +                 ptrn0,tupdim,tuple,vtup)

c          -------------------------------------------------------------
c          Run maxitr GRASP iterations on reduced graph.
c          Largest independent set found is of size nkset and its
c          vertices (of reduced graph) are returned in kset.
c
c          Array tuple is passed as a work array that is called 
c          degv in grspss.
c
c          Array q is passed as a work array that is called 
c          vset in grspss.
c          -------------------------------------------------------------
           call grspss(adjn,alpha,biter0,bkset,deg,tuple,gitr,incdij,
     +                 iniad,invmap,ittype,kset,look4,lstype,lswhen,
     +                 max2a,max2n2,maxn,maxn2,maxitr,nkset,nkset0,
     +                 nnode,nnode0,out,ptrn,rcl,seed,tpl,tupdim,q)

c          -------------------------------------------------------------
c          If the independent set is the largest found so far,
c          unpack it and save it in array kset0.
c          The size of the set is nkset0.
c          -------------------------------------------------------------
           if (nkset+tupdim.gt.nkset0) then
                do 30 i=1,nkset
                     kset(i)=invmap(kset(i))
30              continue
                do 40 k=1,tupdim
                     kset(nkset+k)=vtup(k)
40              continue
                nkset=nkset+tupdim
                nkset0=nkset
                btup=tpl
                biter=biter0
                call cpi4(maxn,nkset,kset,kset0)
c               --------------------------------------------------------
c               If the set found is of the size (look4) sought or
c               greater, return.
c               --------------------------------------------------------
                if (nkset0.ge.look4) goto 60
           endif
50    continue
c     ------------------------------------------------------------------
c     Sort the nkset0 indepedent vertices in kset0 to return.
c     ------------------------------------------------------------------
60    nq=0
      do 70 i=1,nkset0
           call insrtq(maxn2,iq,i,q,nq,kset0(i))
70    continue
      do 80 i=1,nkset0
           call removq(maxn2,iq,k,q,nq,kset0(i))
80      continue
c     ------------------------------------------------------------------
      return
      end





      subroutine mkds(adjn0,deg0,incdij,iniad0,max2a,maxa,maxn,maxn2,
     +                narcs0,nnode0,node1,node2,ptrn0)
c     ------------------------------------------------------------------
c
c     mkds: Builds adjacency data structure for graph, computes degree.
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c     
c          max2a  - array dimension = 2 * maxa 
c          maxa   - array dimension 
c          maxn   - array dimension 
c          maxn2  - array dimension 
c          narcs0 - number of arcs in graph
c          nnode0 - number of nodes in graph
c
c     ------------------------------------------------------------------
      integer   max2a,maxa,maxn,maxn2,narcs0,nnode0
c     ------------------------------------------------------------------
c
c     Passed input arrays:
c     
c          node1  - node 1 of arc
c          node2  - node 2 of arc
c
c     ------------------------------------------------------------------
      integer   node1(maxa),node2(maxa)
c     ------------------------------------------------------------------
c
c     Passed output arrays:
c     
c          adjn0  - node in adjacency list
c          deg0   - degree of node
c          incdij - incidence matrix (in vector form)
c          iniad0 - pointer to first vertex of adjacency list
c          ptrn0  - pointer to next vertex in adjacency list
c
c     ------------------------------------------------------------------
      integer   adjn0(max2a),deg0(maxn),incdij(maxn2),iniad0(maxn),
     +          ptrn0(max2a)
c     ------------------------------------------------------------------
c
c     Local scalars:
c     
c          i      - do loop index
c          index  - do loop index
c          n1     - node 1
c          n2     - node 2
c          ptradj - pointer to first available element in ad
c     
c     ------------------------------------------------------------------
      integer i,index,n1,n2,ptradj
c     ------------------------------------------------------------------
c     Initialize pointers to list of adjacent nodes and degrees.
c     ------------------------------------------------------------------
      ptradj=0
      do 10 i=1,nnode0
           iniad0(i)=0
           deg0(i)=0
10    continue
c     ------------------------------------------------------------------
c     Initialize adjacency matrix.
c     ------------------------------------------------------------------
      do 20 index=1,nnode0*nnode0
           incdij(index)=0
20    continue
c     ------------------------------------------------------------------
c     Build adjacency list and adjacency matrix and compute degrees
c     ------------------------------------------------------------------
      do 30 i=1,narcs0
           n1=node1(i)
           n2=node2(i)
           incdij((n1-1)*nnode0+n2)=1
           incdij((n2-1)*nnode0+n1)=1
           ptradj=ptradj+1
           adjn0(ptradj)=n2
           ptrn0(ptradj)=iniad0(n1)
           iniad0(n1)=ptradj
           deg0(n1)=deg0(n1)+1
           ptradj=ptradj+1
           adjn0(ptradj)=n1
           ptrn0(ptradj)=iniad0(n2)
           iniad0(n2)=ptradj
           deg0(n2)=deg0(n2)+1
30    continue
c     ------------------------------------------------------------------
      return
      end





      subroutine mktup(adj,adjn0,beta,deg0,errcnd,fdmiq,fdmq,incdij,
     +                 iniad0,iq,lowdeg,max2a,max2n2,maxn,maxn2,maxt,
     +                 maxtd,nlowdg,nnode0,ntuple,ptrn0,q,tupdim,tupind,
     +                 tuple,tuplst,vtup)
c     ------------------------------------------------------------------
c
c     mktup: Computes the list of tuples having freedom between
c            minfdm and minfdm + beta * (maxfdm-minfdm), where:
c
c               minfdm = min freedom over all pairs of nodes
c               maxfdm = max freedom over all pairs of nodes.
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          beta   - cut-off parameter
c          max2a  - array dimension
c          max2n2 - array dimension
c          maxn   - array dimension
c          maxn2  - array dimension
c          maxt   - maximum number of elements allowed in freedom heap
c          maxtd  - maximum tuple dimension
c          nlowdg - number of vertices of low degree to be used in
c                   formation of tuples
c          nnode0 - number of nodes in graph
c
c     ------------------------------------------------------------------
      real      beta
      integer   max2a,max2n2,maxn,maxn2,maxt,maxtd,nlowdg,nnode0
c     ------------------------------------------------------------------
c
c     Passed output scalars:
c
c          errcnd - error condition
c          ntuple - number of tuples
c          tupdim - dimension of tuple
c
c     ------------------------------------------------------------------
      integer   errcnd,ntuple,tupdim
c     ------------------------------------------------------------------
c
c     Local scalar:
c
c          sizeq  - size of vertex degree heap
c
c     ------------------------------------------------------------------
      integer   sizeq
c     ------------------------------------------------------------------
c 
c     Passed input arrays:
c
c          adjn0  - vertex in adjacency list (original)
c          deg0   - array of vertex degrees in original graph
c          incdij - graph incidence matrix
c          iniad0 - pointer to start of adjacency list (original)
c          ptrn0  - pointer to next element of adjacency list (original)
c
c     ------------------------------------------------------------------
      integer   adjn0(max2a),deg0(maxn),incdij(maxn2),iniad0(maxn),
     +          ptrn0(max2a)
c     ------------------------------------------------------------------
c 
c     Passed working arrays:
c
c          adj    - indicator array of adjacent nodes
c          fdmiq  - heap index array 
c          fdmq   - heap value array
c          iq     - heap index array
c          lowdeg - vertices of low degree
c          q      - heap value array
c          tupind - index of tuple 
c          tuple  - list of tuples
c          vtup   - tuple vertices 
c
c     ------------------------------------------------------------------
      integer   adj(maxn2),fdmiq(maxt),fdmq(maxt),iq(maxn2),
     +          lowdeg(maxn),q(maxn2),tupind(maxtd),tuple(max2n2),
     +          vtup(maxtd)
c     ------------------------------------------------------------------
c 
c     Passed output array:
c
c          tuplst - list of tuples
c
c     ------------------------------------------------------------------
      integer   tuplst(max2n2)
c     ------------------------------------------------------------------
c
c     Local scalars:
c
c          cutoff - freedom cutoff for tuple construction 
c          dv     - degree of vertex
c          freedm - freedom of tuple
c          i      - do loop index
c          index  - index
c          iv     - index of vertex
c          k      - do loop index
c          maxfdm - maximum freedom
c          minfdm - minimum freedom
c          nfdmq  - number of elements in freedom heap
c          tmkm1  - tupdim*(k-1)
c          tmtpm1 - tupdim*(tpl-1)
c          tpl    - tuple
c
c     ------------------------------------------------------------------
      integer   cutoff,dv,freedm,i,index,iv,k,maxfdm,minfdm,nfdmq,
     +          tmkm1,tmtpm1,tpl
c     ------------------------------------------------------------------
c     Sort nodes, using heap sort, in increasing order of degrees.
c     
c     Retrieve the nlowdg nodes with smallest degrees.  Place them in
c     array lowdeg.
c     ------------------------------------------------------------------
      nfdmq=0
      sizeq=0
      ntuple=0
      if (nlowdg .gt. nnode0) nlowdg=nnode0
      do 10 i=1,nnode0
           call insrtq(maxn2,iq,i,q,sizeq,deg0(i))
10    continue
      do 20 i=1,nlowdg
           call removq(maxn2,iq,iv,q,sizeq,dv)
           lowdeg(i)=iv
20    continue
      minfdm=nnode0
      if (tupdim.gt.nlowdg) tupdim=nlowdg
c     ------------------------------------------------------------------
c     Generate all tuples of lowdeg nodes of size tupdim and compute 
c     freedom for all those that are made up of nodes that are pairwise 
c     nonadjacent.  Sort freedoms, using heap sort, in decreasing order.
c
c     Generate first tuple = (1,2,...,tupdim).
c     ------------------------------------------------------------------
      sizeq=0
      do 30 i=1,tupdim
           tupind(i)=i
c          -------------------------------------------------------------
c          Push i into q.
c          -------------------------------------------------------------
           sizeq=sizeq+1
           q(sizeq)=i
30    continue
      call tupfdm(adj,adjn0,errcnd,fdmiq,fdmq,incdij,iniad0,lowdeg,
     +            max2a,max2n2,maxn2,maxn,maxt,maxtd,minfdm,nfdmq,
     +            nnode0,ntuple,ptrn0,tupdim,tupind,tuple,vtup)

c     ------------------------------------------------------------------
c     begin one level deeper
c     ------------------------------------------------------------------
40    if (sizeq .eq. 0) goto 80
      i=q(sizeq)
      tupind(i)=tupind(i)+1
      if (tupind(i).ne.nlowdg-tupdim+i) then
           do 50 k=i+1,tupdim
                tupind(k)=tupind(i)+k-i
c               --------------------------------------------------------
c               Push k into q.
c               --------------------------------------------------------
                sizeq=sizeq+1
                q(sizeq)=k
50         continue
      else
           do 60 k=i+1,tupdim
                tupind(k)=tupind(i)+k-i
60         continue
           call tupfdm(adj,adjn0,errcnd,fdmiq,fdmq,incdij,iniad0,
     +                 lowdeg,max2a,max2n2,maxn2,maxn,
     +                 maxt,maxtd,minfdm,nfdmq,nnode0,ntuple,ptrn0,
     +                 tupdim,tupind,tuple,vtup)
c          -------------------------------------------------------------
c          Pop i into q.
c          -------------------------------------------------------------
           if (sizeq.eq.0) then
                i=0
           else
                i=q(sizeq)
                sizeq=sizeq-1
           endif

           goto 40
      endif
70    continue 
      call tupfdm(adj,adjn0,errcnd,fdmiq,fdmq,incdij,iniad0,lowdeg,
     +            max2a,max2n2,maxn2,maxn,maxt,maxtd,minfdm,nfdmq,
     +            nnode0,ntuple,ptrn0,tupdim,tupind,tuple,vtup)
      if (tupind(i) .eq. nlowdg-tupdim+i) then
c          -------------------------------------------------------------
c          Pop i into q.
c          -------------------------------------------------------------
           if (sizeq.eq.0) then
                i=0
           else
                i=q(sizeq)
                sizeq=sizeq-1
           endif

           goto 40
      else
           i=q(sizeq)
           tupind(i)=tupind(i)+1
           goto 70
      endif

80    continue
c     ------------------------------------------------------------------
c     Setup the sorted tuple list.
c
c     Get first tuple (with largest freedom) from heap.
c     ------------------------------------------------------------------
      call removq(maxt,fdmiq,tpl,fdmq,nfdmq,freedm)
      index=tupdim*(tpl-1)
      do 90 i=1,tupdim
           tuplst(i)=tuple(index+i)
90    continue
c     ------------------------------------------------------------------
c     Compute cutoff for inclusion in tuple list.
c
c     All but one of the accepted tuples will have freedom between 
c     minfdm and minfdm + beta*(maxfdm-minfdm). One tuple will have
c     the largest freedom less than minfdm + beta*(maxfdm-minfdm).
c     ------------------------------------------------------------------
      maxfdm=-freedm
      cutoff=minfdm+ifix(beta*(maxfdm-minfdm))
c     ------------------------------------------------------------------
c     Get all tuples with freedom >= cutoff.
c     Get one tuple with largest freedom < cutoff.
c
c     Node pairs making up tuple list are placed in array tuplst.
c     ------------------------------------------------------------------
      do 110 k=2,ntuple
           call removq(maxt,fdmiq,tpl,fdmq,nfdmq,freedm)
           tmtpm1=tupdim*(tpl-1)
           tmkm1=tupdim*(k-1)
           do 100 i=1,tupdim
                tuplst(tmkm1+i)=tuple(tmtpm1+i)
100        continue
           if(-freedm.lt.cutoff) goto 120
110   continue
c     ------------------------------------------------------------------
c     Reset number of tuples.
c     ------------------------------------------------------------------
120   ntuple=k
c     ------------------------------------------------------------------
      return
      end






      subroutine tupfdm(adj,adjn0,errcnd,fdmiq,fdmq,incdij,iniad0,
     +                  lowdeg,max2a,max2n2,maxn2,maxn,maxt,maxtd,
     +                  minfdm,nfdmq,nnode0,ntuple,ptrn0,tupdim,tupind,
     +                  tuple,vtup)
c     ------------------------------------------------------------------
c
c     tupfdm: The freedom of a set of vertices is the number of vertices 
c             in the graph that are not adjacent to any member vertex 
c             of the set.
c
c             Compute the freedom of set define by lowdeg (tupind).
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          max2a  - array dimension
c          max2n2 - array dimension
c          maxn2  - array dimension
c          maxn   - array dimension
c          maxt   - array dimension
c          maxtd  - array dimension
c          nnode0 - number of nodes in original graph
c          tupdim - dimension of tuple
c
c     ------------------------------------------------------------------
      integer   max2a,max2n2,maxn2,maxn,maxt,maxtd,nnode0,tupdim
c     ------------------------------------------------------------------
c
c     Passed input/output scalars:
c
c          minfdm - minimum tuple freedom
c          nfdmq  - number of elements in freedopm heap 
c          ntuple - number of tuples generated so far
c
c     ------------------------------------------------------------------
      integer   minfdm,nfdmq,ntuple
c     ------------------------------------------------------------------
c
c     Passed output scalars:
c
c          errcnd - error condition
c
c     ------------------------------------------------------------------
      integer errcnd
c     ------------------------------------------------------------------
c
c     Local scalars:
c
c          freedm - tuple freedom
c          i      - do loop index
c          j      - do loop index
c          k      - do loop index
c          ntim1  - nnode0*(i-1)
c          nvim1  - nnode0*(vtup(i)-1)
c          ptr    - pointer
c          sumadj - sum of adjcent vertices to tuple
c          ttnm1  - tupdim*(ntuple-1)
c
c     ------------------------------------------------------------------
      integer  freedm,i,j,k,ntim1,nvim1,ptr,sumadj,ttnm1
c     ------------------------------------------------------------------
c
c     Passed input arrays:
c
c          adjn0  - adjacent nodes in original graph
c          incdij - incidence matrix
c          iniad0 - pointer to start of list of adjacent vertices
c          lowdeg - indices of vertices of lowest degree
c          ptrn0  - pointer to next element of adjacency list (original)
c          tupind - vertex indices of tuple
c
c     ------------------------------------------------------------------
      integer adjn0(max2a),incdij(maxn2),iniad0(maxn),lowdeg(maxn),
     +        ptrn0(max2a),tupind(maxtd)
c     ------------------------------------------------------------------
c
c     Passed input/output arrays:
c
c          fdmiq  - tuple index in freedom queue
c          fdmq   - tuple value in freedom queue
c
c     ------------------------------------------------------------------
      integer fdmiq(maxt),fdmq(maxt)
c     ------------------------------------------------------------------
c
c     Passed output arrays:
c
c          tuple  - array of tuples
c
c     ------------------------------------------------------------------
      integer tuple(max2n2)
c     ------------------------------------------------------------------
c
c     Passed work arrays:
c
c          adj    - indicator arrays of adjacent nodes
c          vtup   - vertex indices of tuple
c
c     ------------------------------------------------------------------
      integer adj(maxn2),vtup(maxtd)
c     ------------------------------------------------------------------
      do 10 i=1,tupdim
           vtup(i)=lowdeg(tupind(i))
10    continue
c     ------------------------------------------------------------------
c     Check if vertices in set vtup are pairwise nonadjacent.
c     Get next tuple if they are (return).
c     ------------------------------------------------------------------
      do 30 i=1,tupdim-1
           nvim1=nnode0*(vtup(i)-1)
           do 20 j=i+1,tupdim
                if (incdij(nvim1+vtup(j)).eq.1) return
20         continue
30    continue
c     ------------------------------------------------------------------
c     Vertices in set vtup are pairwise nonadjacent. 
c
c     Setup indicator arrays adj of adjacent nodes:
c
c             adj(nnode0*(vtup(i)-1)+u) = 1 if node u is adjacent 
c                                           to vtup(i)
c                                       = 0 otherwise
c
c
c     This array is used to compute freedom of the set vtup.
c     ------------------------------------------------------------------
      do 50 i=1,tupdim
           ntim1=nnode0*(i-1)
           do 40 k=1,nnode0
                adj(ntim1+k)=0
40         continue
50    continue
      do 70 i=1,tupdim
           ptr=iniad0(vtup(i))
           nvim1=nnode0*(vtup(i)-1)
60         if(ptr.gt.0) then
                adj(nvim1+adjn0(ptr))=1
                ptr=ptrn0(ptr)
                goto 60
           endif
70    continue
c     ------------------------------------------------------------------
c     Compute freedom of tuple (vtup) and insert into heap 
c     for sorting.
c
c     Find minimum value of freedom as freedoms are computed.
c     ------------------------------------------------------------------
      freedm=0
      do 100 k=1,nnode0
           do 80 i=1,tupdim
                if (k .eq. vtup(i)) goto 100 
80         continue
           sumadj=0
           do 90 i=1,tupdim
                sumadj=sumadj+adj(nnode0*(vtup(i)-1)+k)
90         continue
           if(sumadj.eq.0) freedm=freedm+1
100   continue
      ntuple=ntuple+1
      ttnm1=tupdim*(ntuple-1)
      do 110 i=1,tupdim
           tuple(ttnm1+i)=vtup(i)
110   continue
      if (ntuple .gt. maxt) then
           errcnd = 100
      else
           call insrtq(maxt,fdmiq,ntuple,fdmq,nfdmq,-freedm)
           if(freedm.lt.minfdm) minfdm=freedm
      endif
c     ------------------------------------------------------------------
      return
      end






      subroutine mkrgrf(adjn,adjn0,deg,fwdmap,iniad,iniad0,invmap,
     +                 max2a,max2n2,maxn,maxtd,narcs,nnode,nnode0,ptrn,
     +                 ptrn0,tupdim,vthere,vtup)
c     ------------------------------------------------------------------
c
c     mkrgrf:  Deletes nodes in vtup from the original graph generating 
c              the reduced graph for the GRASP iterations.
c              
c              Reduced graph is return in data structure
c              (nnode,narcs,iniad,ptrn,adjn,deg).
c
c              Arrays (fwdmap and invmap) that map corresponding 
c              vertices in original and reduced graphs are set up.
c
c                + Node i in reduced  graph is invmap(i) 
c                  in original graph. 
c
c                + Node i in original graph is fwdmap(i)
c                  in reduced  graph. 
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          max2a  - array dimension
c          max2n2 - array dimension
c          maxn   - array dimension
c          maxtd  - array dimension
c          nnode0 - number of nodes in original graph
c          tupdim - dimension of tuple
c
c     ------------------------------------------------------------------
      integer max2a,max2n2,maxn,maxtd,nnode0,tupdim
c     ------------------------------------------------------------------
c
c     Passed output scalars:
c
c          narcs  - number of arcs in reduced graph
c          nnode  - number of nodes in reduced graph
c
c     ------------------------------------------------------------------
      integer narcs,nnode
c     ------------------------------------------------------------------
c
c     Local scalars:
c
c          adjn0p - auxiliary variable
c          i      - do loop index
c          k      - do loop index
c          nz     - number of eliminated vertices in reduced graph
c          ptr    - pointer
c          ptrad  - pointer
c          ri     - auxiliary variable
c
c     ------------------------------------------------------------------
      integer adjn0p,i,k,nz,ptr,ptrad,ri
c     ------------------------------------------------------------------
c
c     Passed input arrays:
c
c          adjn0  - vertices of adjacency list of original graph
c          iniad0 - pointers to adjacency lists of original graph
c          ptrn0  - pointer to next element in adjacency list of 
c                   original graph
c          vtup   - vertices of tuple
c
c     ------------------------------------------------------------------
      integer adjn0(max2a),iniad0(maxn),ptrn0(max2a),
     +        vtup(maxtd)
c     ------------------------------------------------------------------
c
c     Passed output arrays:
c
c          adjn   - vertices of adjacency list of reduced graph
c          deg    - degrees of reduced graph
c          fwdmap - map of vertices to original graph
c          iniad  - pointers to adjacency lists of reduced graph
c          invmap - map of vertices to reduced graph
c          ptrn   - pointer to next element in adjacency list of 
c                   reduced graph
c
c     ------------------------------------------------------------------
      integer adjn(max2a),deg(maxn),fwdmap(maxn),iniad(maxn),
     +        invmap(maxn),ptrn(max2a)
c     ------------------------------------------------------------------
c
c     Passed work arrays:
c
c          vthere - indicates if vertex is in reduced graph
c
c     ------------------------------------------------------------------
      integer vthere(max2n2)
c     ------------------------------------------------------------------
c     Prepare forward (fwdmap) and inverse (invmap) maps.
c     ------------------------------------------------------------------
      do 10 i=1,nnode0
           vthere(i)=1
10    continue
      if (tupdim .gt. 0) then
           do 20 k=1,tupdim
                vthere(vtup(k))=0
20         continue
           do 40 k=1,tupdim
                ptr=iniad0(vtup(k))
30              if (ptr.gt.0) then
                     vthere(adjn0(ptr))=0
                     ptr=ptrn0(ptr)
                     goto 30
                endif
40         continue
           nz=0
           do 50 i=1,nnode0
                if (vthere(i).eq.0) then
                     nz=nz+1
                     fwdmap(i)=0
                else
                     fwdmap(i)=i-nz
                     invmap(i-nz)=i
                endif
50         continue
      else
           nz=0
           do 60 i=1,nnode0
                fwdmap(i)=i
                invmap(i)=i
60         continue
      endif


c     ------------------------------------------------------------------
c     Compress adjacency list.
c     ------------------------------------------------------------------
      nnode=nnode0-nz
      do 70 i=1,nnode
          iniad(i)=0
          deg(i)=0
70    continue
      narcs=0
      ptrad=0
      do 90 i=1,nnode0
           if (vthere(i).eq.1) then
                ptr=iniad0(i)
80              if (ptr.gt.0) then
                     ri=fwdmap(i)
                     adjn0p=adjn0(ptr)
                     if (vthere(adjn0p).eq.1) then
                          ptrad=ptrad+1
                          ptrn(ptrad)=iniad(ri)
                          iniad(ri)=ptrad
                          adjn(ptrad)=fwdmap(adjn0p)
                          deg(ri)=deg(ri)+1
                          narcs=narcs+1
                     endif
                     ptr=ptrn0(ptr)
                     goto 80
                endif
           endif
90    continue
      narcs=narcs/2
c     ------------------------------------------------------------------
      return
      end






      subroutine grspss(adjn,alpha,biter0,bkset,deg,degv,gitr,incdij,
     +                  iniad,invmap,ittype,kset,look4,lstype,lswhen,
     +                  max2a,max2n2,maxn,maxn2,maxitr,nkset,nkset0,
     +                  nnode,nnode0,out,ptrn,rcl,seed,tpl,tupdim,vset)
c     ------------------------------------------------------------------
c
c     grspss: Applies GRASP to construct an independent set in graph 
c             G=(iniad,ptrn,adjn,deg).  
c 
c             Set of size nkset is returned in array kset.
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          alpha  - GRASP RCL parameter
c          ittype - iteration summary type
c          look4  - size of independent set sought
c          lstype - type of local search
c          lswhen - indicates when local search is done
c          max2a  - array dimension = 2*maxa
c          max2n2 - array dimension = 2*maxn2
c          maxitr - maximum number of GRASP iterations
c          maxn   - array dimension
c          maxn2  - array dimension
c          nkset0 - number of elements in best independent set 
c                   of original graph
c          nnode  - number of nodes in reduced graph
c          nnode0 - number of nodes in original graph
c          out    - Fortran output device
c          tpl    - tuple under consideration
c          tupdim - dimension of tuple vector
c
c     ------------------------------------------------------------------
      real    alpha
      integer ittype,look4,lstype,lswhen,max2a,max2n2,
     +        maxitr,maxn,maxn2,nkset0,nnode,nnode0,out,tpl,tupdim
c     ------------------------------------------------------------------
c
c     Passed input/output scalars:
c
c          seed   - seed for pseudo random number generator
c
c     ------------------------------------------------------------------
      integer seed
c     ------------------------------------------------------------------
c
c     Passed output scalars:
c
c          biter0 - number of iterations to find best set
c          gitr   - do loop index (GRASP iterations)
c          nkset  - number of vertices in independent set of reduced
c                   graph
c
c     ------------------------------------------------------------------
      integer biter0,gitr,nkset
c     ------------------------------------------------------------------
c
c     Passed input arrays:
c
c          adjn   - vertex in adjacency list (copy)
c          deg    - degree of vertex (copy)
c          incdij - graph incidence matrix
c          iniad  - pointer to start of adjacency list (copy)
c          invmap - node i in reduced graph is invmap(i) in original 
c                   graph, passed to mkrgrf, grspss
c          ptrn   - pointer to next element of adjacency list (copy) 
c
c     ------------------------------------------------------------------
      integer adjn(max2a),deg(maxn),incdij(maxn2),iniad(maxn),
     +        invmap(maxn),ptrn(max2a)
c     ------------------------------------------------------------------
c
c     Passed work arrays:
c
c          bkset  - set of vertices in independent set
c          degv   - copy of original array deg, which is overwritten
c          rcl    - restricted candidate list 
c          vset   - work array passed to local search procedures 
c
c     ------------------------------------------------------------------
      integer bkset(maxn),degv(max2n2),rcl(maxn),vset(maxn2)
c     ------------------------------------------------------------------
c
c     Passed output arrays:
c
c          kset   - set of independent nodes
c
c     ------------------------------------------------------------------
      integer kset(maxn)
c     ------------------------------------------------------------------
c
c     Local scalars:
c
c          avgph1 - average size of set constructed in phase 1 of GRASP
c          nbkset - size of largest set found in phase 1
c
c     ------------------------------------------------------------------
      double precision  avgph1
      integer nbkset
c     ------------------------------------------------------------------
c     Initialize independent set size and average size of phase 1
c     independent sets.
c     ------------------------------------------------------------------
      nkset=nkset0-tupdim
      avgph1=0.0d0
c     ------------------------------------------------------------------
c     Save vertex degrees in degv.  Array deg is over-written and must
c     be recovered to start each GRASP iteration.
c     ------------------------------------------------------------------
      call cpi4(maxn,nnode,deg,degv)
c     ------------------------------------------------------------------
c     Do maxitr GRASP iterations.
c     ------------------------------------------------------------------
      do 20 gitr=1,maxitr
c          -------------------------------------------------------------
c          GRASP construction phase (phase 1).
c          -------------------------------------------------------------
           call build(adjn,alpha,bkset,deg,iniad,max2a,maxn,nbkset,
     +                nnode,ptrn,rcl,seed)
c          -------------------------------------------------------------
c          Store best solution, if found.
c          -------------------------------------------------------------
           if (nbkset .gt. nkset) then
                nkset=nbkset
                biter0=gitr
                call cpi4(maxn,nkset,bkset,kset)
                if (ittype.ge.1) then
                     write(out,10) tpl,gitr,nbkset+tupdim,nkset+tupdim
10                   format(' tuple: ',i5,'   itr: ',i7,'    size: ',i6,
     +                      '    best: ',i6)
                endif
                if (nkset+tupdim .ge. look4) return
           else
                if ((ittype.eq.2) .and. (gitr .eq. 1)) then
                     write(out,10) tpl,gitr,nbkset+tupdim,nkset+tupdim
                endif
           endif
c          -------------------------------------------------------------
c          Update average of phase 1 solutions.
c          -------------------------------------------------------------
           avgph1=(avgph1*(gitr-1)+nbkset)/gitr

c          -------------------------------------------------------------
c          GRASP local search phase.
c
c          If lstype = 0, there is no local search.
c                    = 1, local search is (1,2)-exchange.
c                    = 2, local search is (2,3)-exchange.
c
c          If lswhen = 0, local search is done only if phase 1 solution
c                         is better than the average phase 1 solution.
c                    = 1, local search is no matter what phase 1 
c                         solution was.
c          -------------------------------------------------------------
           if (lstype .ne. 0) then
                if (lswhen.eq.0) then
                     if (nbkset .ge. avgph1) then
                          if (lstype .eq. 1) then
c                              -----------------------------------------
c                              (1,2)-exchange local search
c                              -----------------------------------------
                               call gls12(adjn,bkset,incdij,iniad,
     +                                    invmap,max2a,maxn,maxn2,
     +                                    nbkset,nnode,nnode0,
     +                                    ptrn,vset)
                          else
c                              -----------------------------------------
c                              (2,3)-exchange local search
c                              -----------------------------------------
                               call gls23(adjn,bkset,incdij,iniad,
     +                                    invmap,max2a,maxn,maxn2,
     +                                    nbkset,nnode,nnode0,ptrn,
     +                                    vset)
                          endif
                     endif
                else
                     if (lstype .eq. 1) then
c                         ----------------------------------------------
c                         (1,2)-exchange local search
c                         ----------------------------------------------
                          call gls12(adjn,bkset,incdij,iniad,
     +                               invmap,max2a,maxn,maxn2,
     +                               nbkset,nnode,nnode0,ptrn,vset)
                     else
c                         ----------------------------------------------
c                         (2,3)-exchange local search
c                         ----------------------------------------------
                          call gls23(adjn,bkset,incdij,iniad,invmap,
     +                               max2a,maxn,maxn2,nbkset,nnode,
     +                               nnode0,ptrn,vset)
                     endif
                endif
           endif

c          -------------------------------------------------------------
c          Store best solution, if found.
c          -------------------------------------------------------------
           if (nbkset .gt. nkset) then
                nkset=nbkset
                call cpi4(maxn,nkset,bkset,kset)
                if (ittype.ge.1) then
                     write(out,10) tpl,gitr,nbkset+tupdim,nkset+tupdim
                endif
                if (nkset+tupdim .ge. look4) return
           else
                if ((ittype.eq.2) .and. (gitr .ge. 2)) then
                     write(out,10) tpl,gitr,nbkset+tupdim,nkset+tupdim
                endif
           endif
c          -------------------------------------------------------------
c          Restore array deg to its original state.
c          -------------------------------------------------------------
           call cpi4(maxn,nnode,degv,deg)
20    continue
c     gitr=maxitr
c     ------------------------------------------------------------------
      return
      end






      subroutine chkerr(alpha,beta,errcnd,look4,lstype,lswhen,
     +                  maxitr,maxtd,nlowdg,nnode0,seed,
     +                  tupdim)
c     ------------------------------------------------------------------
c
c     chkerr: Check if parameters have been set correctly.
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          alpha  - GRASP RCL parameter
c          beta   - tuple list size restriction parameter 
c          look4  - size of independent set sought
c          lstype - type of local search
c          lswhen - indicates when local search is done
c          maxitr - maximum number of GRASP iterations
c          maxtd  - array dimension
c          nlowdg - number of low degree vertices used to make tuples 
c          nnode0 - number of nodes in original graph
c          seed   - seed for pseudo random number generator
c          tupdim - dimension of tuple vector
c
c     ------------------------------------------------------------------
      real     alpha,beta
      integer  look4,lstype,lswhen,maxitr,
     +         maxtd,nlowdg,nnode0,seed,tupdim
c     ------------------------------------------------------------------
c
c     Passed input/output scalars:
c
c          errcnd - error condition
c
c     ------------------------------------------------------------------
      integer  errcnd
c     ------------------------------------------------------------------
      if (alpha. lt. 0.) then
           errcnd = 140
           return
      endif
      if (alpha. gt. 1.) then
           errcnd = 141
           return
      endif
      if (beta. lt. 0.) then
           errcnd = 150
           return
      endif
      if (beta. gt. 1.) then
           errcnd = 151
           return
      endif
      if (look4. gt. nnode0) then
           errcnd = 160
           return
      endif
      if ((lstype. ne. 0) .and. (lstype .ne. 1). and. 
     +    (lstype .ne. 2)) then
           errcnd = 170
           return
      endif
      if ((lswhen. ne. 0) .and. (lswhen .ne. 1)) then
           errcnd = 171
           return
      endif
      if (maxitr .le. 0) then
           errcnd = 180
           return
      endif
      if ((tupdim.ge.1) .and. (nlowdg .le. 0)) then
           errcnd = 190
           return
      endif
      if (seed.lt.1) then
           errcnd = 200
           return
      endif
      if (seed.gt.2147483647) then
           errcnd = 201
           return
      endif
      if (tupdim.lt.0) then
           errcnd = 210
           return
      endif
      if (tupdim.gt.maxtd) then
           errcnd = 211
           return
      endif
c     ------------------------------------------------------------------
      return
      end






      subroutine build(adjn,alpha,bkset,deg,iniad,max2a,maxn,nbkset,
     +                 nnode,ptrn,rcl,seed)
c     ------------------------------------------------------------------
c
c     build: GRASP construction phase (phase 1)
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          alpha  - GRASP RCL parameter
c          max2a  - array dimension = 2*maxa
c          maxn   - array dimension
c          nnode  - number of nodes in reduced graph
c
c     ------------------------------------------------------------------
      real    alpha
      integer max2a,maxn,nnode
c     ------------------------------------------------------------------
c
c     Passed input/output scalars:
c
c          nbkset - size of independent set found
c          seed   - seed for pseudo random number generator
c
c     ------------------------------------------------------------------
      integer nbkset,seed
c     ------------------------------------------------------------------
c
c     Passed input arrays:
c
c          adjn   - vertex in adjacency list (copy)
c          deg    - degree of vertex (copy)
c          iniad  - pointer to start of adjacency list (copy)
c          ptrn   - pointer to next element of adjacency list (copy) 
c
c     ------------------------------------------------------------------
      integer adjn(max2a),deg(maxn),iniad(maxn),ptrn(max2a)
c     ------------------------------------------------------------------
c
c     Passed output arrays:
c
c          bkset  - set of independent nodes
c
c     ------------------------------------------------------------------
      integer bkset(maxn)
c     ------------------------------------------------------------------
c
c     Passed working arrays:
c
c          rcl    - restricted candidate list
c
c     ------------------------------------------------------------------
      integer rcl(maxn)
c     ------------------------------------------------------------------
c
c     Local scalars and functions:
c
c          cutoff - cutoff for RCL
c          degi   - auxiliary scalar
c          i      - do loop index
c          maxdeg - maximum degree
c          mindeg - maximum degree
c          nrcl   - number of elements in RCL
c          nselct - vertex selected from RCL
c          ptr    - ptr
c          ptrv   - ptrv
c          u      - vertex
c          v      - vertex
c          vfound - indicator variable
c          vrtx   - vertex
c
c     ------------------------------------------------------------------
      integer cutoff,degi,i,maxdeg,mindeg,nrcl,nselct,ptr,ptrv
      integer u,v,vfound,vrtx
c     ------------------------------------------------------------------
c     Initialize size of indepedendent set
c     ------------------------------------------------------------------
      nbkset=0
c     ------------------------------------------------------------------
c     Independent set is constructed, one vertex at a time.
c     ------------------------------------------------------------------
10    continue
c          -------------------------------------------------------------
c          Set indicator that candidate vertex exists to false. 
c          -------------------------------------------------------------
           vfound=0
c          -------------------------------------------------------------
c          Find minimum and maximum degrees of all candidate nodes.
c          In process, determine if construction phase is done, i.e.
c          if there no more candidates to be included in the indep set.
c          -------------------------------------------------------------
           mindeg=nnode
           maxdeg=0
           do 20 i=1,nnode
                degi=deg(i)
                if (degi.gt.-1) then
c                    ---------------------------------------------------
c                    At least one vertex is a candidate to be in the
c                    independent set.
c                    ---------------------------------------------------
                     vfound=1
                     if (degi.lt.mindeg) then
                          mindeg=degi
                     endif
                     if (degi.gt.maxdeg) then
                          maxdeg=degi
                     endif
                endif
20         continue
c          -------------------------------------------------------------
c          Test if construction phase is done.
c          -------------------------------------------------------------
           if (vfound.eq.0) goto 60
c          -------------------------------------------------------------
c          Construction phase not done yet.
c
c          Compute cut-off for restricted candidate list (RCL).
c
c          A vertex will be put in the RCL if its degree with respect 
c          to vertices no in the independent set is between mindeg and 
c          mindeg + alpha*(maxdeg-mindeg).
c          -------------------------------------------------------------
           cutoff=mindeg + ifix(alpha*(maxdeg-mindeg))
c          -------------------------------------------------------------
c          Build the RCL.  Place the nrcl vertices in array rcl.
c          -------------------------------------------------------------
           nrcl=0
           do 30 i=1,nnode
                degi=deg(i)
                if (degi.gt.-1 .and. degi.le.cutoff) then
                     nrcl=nrcl+1
                     rcl(nrcl)=i
                endif
30         continue
c          -------------------------------------------------------------
c          Select a candidate at random from RCL to be put in the
c          independent set.
c
c          There are nrcl elements in array rcl.   The value nselct
c          is a pseudo random integer between 1 and nrcl, inclusive.
c          Subroutine randp returns a seed between 1 and 2147483647.
c          -------------------------------------------------------------
           call randp(seed)
           nselct=1+seed/(2147483647/nrcl)
c          -------------------------------------------------------------
c          Increment size of independent set and put vertex
c          vrtx = rcl(nselct) in independent set array bkset.
c          -------------------------------------------------------------
           nbkset=nbkset+1
           vrtx=rcl(nselct)
           bkset(nbkset)=vrtx
c          -------------------------------------------------------------
c          Adaptive component of GRASP. Adjust degrees of remaining 
c          nodes, thus changing the greedy function.
c
c          By convention, a vertex with degree < 0 is not a candidate.
c          Either it is in the independent set or is adjacent to a
c          vertex that is in the independent set.
c          -------------------------------------------------------------
           deg(vrtx)=-1

c          -------------------------------------------------------------
c          Scan adjacent vertices of independent vertex vrtx.
c          -------------------------------------------------------------
           ptr=iniad(vrtx)
40         if (ptr.gt.0) then
                v=adjn(ptr)
c               --------------------------------------------------------
c               If adjacent vertex (v) is not in the independent set,
c               scan its neighbors and decrease their degrees by one.
c               --------------------------------------------------------
                if (deg(v).ge.0) then
                    deg(v)=-1
                    ptrv=iniad(v)
50                  if(ptrv.gt.0) then
                         u=adjn(ptrv)
                         deg(u)=deg(u)-1
                         ptrv=ptrn(ptrv)
                         goto 50
                    endif
                endif
                ptr=ptrn(ptr)
                goto 40
           endif
c          -------------------------------------------------------------
c          Begin a new construction phase iteration.
c          -------------------------------------------------------------
           goto 10
60    continue
c     ------------------------------------------------------------------
      return
      end






      subroutine gls12(adjn,bkset,incdij,iniad,invmap,
     +                 max2a,maxn,maxn2,
     +                 nbkset,nnode,
     +                 nnode0,ptrn,vset)
c     ------------------------------------------------------------------
c
c     gls12: GRASP local search phase: (1,2)-exchange
c
c            For each vertex v in independent set S, examine all pairs
c            of nonadjacent vertices (u,w) not in the independent and
c            check if S - v + u + w is an independent set.  
c
c            If so, remove v from S and add u and w to S, thus 
c            increasing the size of the independent set by one.  
c
c            Repeat until local optimum (with respect to this 
c            neighborhood structure) is found.
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          max2a  - array dimension = 2*maxa
c          maxn   - array dimension
c          maxn2  - array dimension
c          nnode  - number of nodes in reduced graph
c          nnode0 - number of nodes in original graph
c
c     ------------------------------------------------------------------
      integer max2a,maxn,maxn2,nnode,nnode0
c     ------------------------------------------------------------------
c
c     Passed input/output scalars:
c
c          nbkset - size of independent set found
c
c     ------------------------------------------------------------------
      integer nbkset
c     ------------------------------------------------------------------
c
c     Passed input arrays:
c
c          adjn   - vertex in adjacency list (copy)
c          incdij - graph incidence matrix
c          iniad  - pointer to start of adjacency list (copy)
c          invmap - node i in reduced graph is invmap(i) in original 
c                   graph, passed to mkrgrf, grspss
c          ptrn   - pointer to next element of adjacency list (copy) 
c
c     ------------------------------------------------------------------
      integer adjn(max2a),incdij(maxn2),iniad(maxn),invmap(maxn),
     +        ptrn(max2a)
c     ------------------------------------------------------------------
c
c     Passed input/output array:
c
c          bkset  - set of vertices in independent set
c
c     ------------------------------------------------------------------
      integer bkset(maxn)
c     ------------------------------------------------------------------
c
c     Passed work array:
c
c          vset   - adjacency indicator
c
c     ------------------------------------------------------------------
      integer vset(maxn2)
c     ------------------------------------------------------------------
c
c     Local scalars:
c
c          bksetk - auxiliary variable = bkset(k) 
c          i      - do loop index
c          im1tn  - (invmap(u1)-1)*nnode0
c          index  - array index
c          invu2  - auxiliary variable = invmap(u2)
c          k      - do loop index
c          k1     - do loop index
c          k2     - do loop index
c          nvset  - number of nonadjacent nodes in vset
c          ptr    - pointer
c          u1     - vertex
c          u2     - vertex
c          v1     - vertex
c          w      - vertex
c
c     ------------------------------------------------------------------
      integer bksetk,i,im1tn,index,invu2,k,k1,k2,nvset,ptr,u1,
     +        u2,v1,w
c     ------------------------------------------------------------------
c     Local search of neighborhood N(S) of solution S starts.
c     ------------------------------------------------------------------
10    continue
c          -------------------------------------------------------------
c          For all vertices v1 in current indep set (bkset), check
c          if there exist nonadjacent vertices (u1,u2) such that u1
c          and u2 are not adjacent to bkset-v1.
c          -------------------------------------------------------------
           do 80 i=1,nbkset
                v1=bkset(i)
c               --------------------------------------------------------
c               Initialize array that indicates whether a node is
c               adjacent to set bkset-v1.
c               --------------------------------------------------------
                do 20 k=1,nnode
                    vset(k)=0
20              continue
c               --------------------------------------------------------
c               Mark node if it is adjacent to node in bkset-v1
c               --------------------------------------------------------
                do 40 k=1,nbkset
c                    ---------------------------------------------------
c                    If node is v1, get next node in bkset.
c                    ---------------------------------------------------
                     if (k.ne.i) then
                          bksetk=bkset(k)
c                         ----------------------------------------------
c                         Mark bkset node as adjacent and scan its 
c                         neighbors, marking them as adjacent.
c                         ----------------------------------------------
                          vset(bksetk)=1
                          ptr=iniad(bksetk)
30                        if(ptr.gt.0) then
                               w=adjn(ptr)
                               vset(w)=1
                               ptr=ptrn(ptr)
                               goto 30
                          endif
                     endif
40              continue
c               --------------------------------------------------------
c               Put all nvset nonadjacent nodes in vset.
c               --------------------------------------------------------
                nvset=0
                do 50 k=1,nnode
                     if (k.ne.v1 .and. vset(k).eq.0) then
                          nvset=nvset+1
                          vset(nvset)=k
                     endif
50              continue
c               --------------------------------------------------------
c               If there are enough vertices to improve independent set
c               (i.e. nvset .ge. 2), see if improvement occurs.
c               --------------------------------------------------------
                if (nvset.ge.2) then
c                    ---------------------------------------------------
c                    For all pairs (u1,u2) in vset, check if u1 and u2
c                    are nonadjacent.  
c
c                    If so, remove v1 from S and add
c                    u1 and u2 to S, increasing size of independent set.
c                    ---------------------------------------------------
                     do 70 k1=1,nvset-1
                          u1=vset(k1)
                          im1tn=(invmap(u1)-1)*nnode0
                          do 60 k2=k1+1,nvset
                                u2=vset(k2)
c                               ----------------------------------------
c                               Check if u1 and u2 are nonadjacent.
c                               ----------------------------------------
                                invu2=invmap(u2)
                                index=im1tn+invu2
                                if (incdij(index).ne.1) then
c                                    -----------------------------------
c                                    Nodes u1 and u2 are nonadjacent:
c                                    remove v1 from the independent set 
c                                    and add u1 and u2 to set, thus 
c                                    increasing size of independent set 
c                                    by one.
c                                    -----------------------------------
                                     bkset(i)=u1
                                     nbkset=nbkset+1
                                     bkset(nbkset)=u2
c                                    -----------------------------------
c                                    Restart local search in the
c                                    neighborhood of new solution.
c                                    -----------------------------------
                                     goto 10
                               endif
60                        continue
70                   continue
                endif
80         continue
      continue
c     ------------------------------------------------------------------
      return
      end






      subroutine gls23(adjn,bkset,incdij,iniad,invmap,
     +                 max2a,maxn,maxn2,
     +                 nbkset,nnode,nnode0,ptrn,vset)
c     ------------------------------------------------------------------
c
c     gls23: GRASP local search phase: (2,3)-exchange
c
c            For each pair of verteices (v,t) in independent set S, 
c            examine all triples of parwise nonadjacent vertices (u,w,y)
c            not in the independent and check if S - v - t + u + w + y 
c            is an independent set.  
c
c            If so, remove v and t from S and add u, w, adn y to S, thus 
c            increasing the size of the independent set by one.  
c
c            Repeat until local optimum (with respect to this 
c            neighborhood structure) is found.
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          max2a  - array dimension = 2*maxa
c          maxn   - array dimension
c          maxn2  - array dimension
c          nnode  - number of nodes in reduced graph
c          nnode0 - number of nodes in original graph
c
c     ------------------------------------------------------------------
      integer max2a,maxn,maxn2,nnode,nnode0
c     ------------------------------------------------------------------
c
c     Passed input/output scalars:
c
c          nbkset - size of independent set found
c
c     ------------------------------------------------------------------
      integer nbkset
c     ------------------------------------------------------------------
c
c     Passed input arrays:
c
c          adjn   - vertex in adjacency list (copy)
c          incdij - graph incidence matrix
c          iniad  - pointer to start of adjacency list (copy)
c          invmap - node i in reduced graph is invmap(i) in original 
c                   graph, passed to mkrgrf, grspss
c          ptrn   - pointer to next element of adjacency list (copy) 
c
c     ------------------------------------------------------------------
      integer adjn(max2a),incdij(maxn2),iniad(maxn),invmap(maxn),
     +        ptrn(max2a)
c     ------------------------------------------------------------------
c
c     Passed input/output array:
c
c          bkset  - set of vertices in independent set
c
c     ------------------------------------------------------------------
      integer bkset(maxn)
c     ------------------------------------------------------------------
c
c     Passed work array:
c
c          vset   - adjacency indicator
c
c     ------------------------------------------------------------------
      integer vset(maxn2)
c     ------------------------------------------------------------------
c
c     Local scalars:
c
c          bksetk - auxiliary variable = bkset(k) 
c          i      - do loop index
c          index  - array index
c          invu1  - auxiliary variable = invmap(u1)
c          invu2  - auxiliary variable = invmap(u2)
c          invu3  - auxiliary variable = invmap(u3)
c          j      - do loop index
c          k      - do loop index
c          k1     - do loop index
c          k2     - do loop index
c          k3     - do loop index
c          nvset  - number of nonadjacent nodes in vset
c          ptr    - pointer
c          u1     - vertex
c          u2     - vertex
c          u3     - vertex
c          w      - vertex
c
c     ------------------------------------------------------------------
      integer bksetk,i,index,invu1,invu2,invu3,j,k,k1,k2,k3,
     +        nvset,ptr,u1,u2,u3,w
c     ------------------------------------------------------------------
c     Local search of neighborhood N(S) of solution S starts.
c     ------------------------------------------------------------------
10    continue
c          -------------------------------------------------------------
c          For all pairs of vertices (v1,v2) in current indep set 
c          (v1 = bkset(i), v2 = bkset(j), check if there exist pairwise 
c          nonadjacent vertices (u1,u2,u3) such that u1, u2, and and u3
c          are not adjacent to bkset-v1-v2.
c          -------------------------------------------------------------
           do 100 i=1,nbkset-1
                do 90 j=i+1,nbkset
c                    ---------------------------------------------------
c                    Initialize array that indicates whether a node is
c                    adjacent to set bkset-v1-v2.
c                    ---------------------------------------------------
                     do 20 k=1,nnode
                          vset(k)=0
20                   continue
c                    ---------------------------------------------------
c                    Mark node if it is adjacent to node in bkset-v1-v2.
c                    ---------------------------------------------------
                     do 40 k=1,nbkset
c                         ----------------------------------------------
c                         Mark bkset node as adjacent.
c                         ----------------------------------------------
                          bksetk=bkset(k)
                          vset(bksetk)=1
c                         ----------------------------------------------
c                         If node is not v1 nor v2, scan neighbors of 
c                         bkset node, marking them as adjacent.
c                         ----------------------------------------------
                          if (k.ne.i .and. k.ne.j) then
                               ptr=iniad(bksetk)
30                             if(ptr.gt.0) then
                                    w=adjn(ptr)
                                    vset(w)=1
                                    ptr=ptrn(ptr)
                                    goto 30
                               endif
                          endif
40                   continue
c                    ---------------------------------------------------
c                    Put all nvset nonadjacent nodes in vset.
c                    ---------------------------------------------------
                     nvset=0
                     do 50 k=1,nnode
                          if (vset(k).eq.0) then
                               nvset=nvset+1
                               vset(nvset)=k
                          endif
50                   continue
c                    ---------------------------------------------------
c                    If there are enough vertices to improve independent
c                    set (i.e. nvset .ge. 3), see if improvement occurs.
c                    ---------------------------------------------------
                     if (nvset.ge.3) then
                          do 80 k1=1,nvset-2
c                            -------------------------------------------
c                            For all triples (u1,u2,u3) in vset, check
c                            if u1, u2, and u3 are pairwise 
c                            nonadjacent.  
c
c                            If so, remove (v1,v2) from S and add
c                            (u1,u2,u3) to S, thus increasing size of 
c                            the independent set by one.
c                            -------------------------------------------
                             u1=vset(k1)
                             invu1=invmap(u1)
                             do 70 k2=k1+1,nvset-1
                                u2=vset(k2)
                                invu2=invmap(u2)
                                index=(invu1-1)*nnode0+invu2
                                do 60 k3=k2+1,nvset
                                   u3=vset(k3)
                                   invu3=invmap(u3)
c                                  -------------------------------------
c                                  Check if u1, u2, and u3 are  pairwise
c                                  nonadjacent. 
c                                  -------------------------------------
                                   if (incdij(index).eq.1) goto 60
                                   index=(invu1-1)*nnode0+invu3
                                   if (incdij(index).eq.1) goto 60
                                   index=(invu2-1)*nnode0+invu3
                                   if (incdij(index).eq.1) goto 60
c                                  -------------------------------------
c                                  Nodes u1, u2, and u3 are pairwise
c                                  nonadjacent.
c
c                                  Remove v1, v2 from independent set 
c                                  and add u1, u2, and u3, thus 
c                                  increasing size of independent set 
c                                  by one.
c                                  -------------------------------------
                                   bkset(i)=u1
                                   bkset(j)=u2
                                   nbkset=nbkset+1
                                   bkset(nbkset)=u3
c                                  -------------------------------------
c                                  Restart local search in neighborhood 
c                                  of new solution.
c                                  -------------------------------------
                                   goto 10
60                              continue
70                           continue
80                        continue
                     endif
90              continue
100        continue
      continue
c     ------------------------------------------------------------------
      return
      end






      subroutine randp(ix)
c     ------------------------------------------------------------------
c
c     randp: Portable pseudo-random number generator.
c            Reference: L. Schrage, "A More Portable Fortran
c            Random Number Generator", ACM Transactions on
c            Mathematical Software, Vol. 2, No. 2, (June, 1979).
c
c     ------------------------------------------------------------------

      integer a,p,ix,b15,b16,xhi,xalo,leftlo,fhi,k
      data a/16807/,b15/32768/,b16/65536/,p/2147483647/

      xhi=ix/b16
      xalo=(ix-xhi*b16)*a
      leftlo=xalo/b16
      fhi=xhi*a+leftlo
      k=fhi/b15
      ix=(((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
      if (ix.lt.0) ix=ix+p
c     ------------------------------------------------------------------
      return
      end







      subroutine insrtq(dimq,iq,iv,q,sizeq,v)
c     ------------------------------------------------------------------
c
c     insrtq: Insert an element (v,iv) into a queue (q,iq).
c
c     ------------------------------------------------------------------
c     Passed input scalars:
c
c         dimq   - array dimension
c         iv     - heap element (index)
c         v      - heap element (value)
c
c     ------------------------------------------------------------------
      integer dimq,iv,v
c     ------------------------------------------------------------------
c     Passed input/output scalars:
c
c         sizeq  - size of heap
c
c     ------------------------------------------------------------------
      integer sizeq
c     ------------------------------------------------------------------
c     Passed input/output arrays:
c
c         iq     - heap (index)
c         q      - heap (value)
c
c     ------------------------------------------------------------------
      integer iq(dimq),q(dimq)
c     ------------------------------------------------------------------
c     Local scalars:
c
c         sq     - temporary size of heap
c         tsz    - temporary variable (sq/2)
c
c     ------------------------------------------------------------------
      integer sq,tsz

c     ------------------------------------------------------------------
c     Insert element into heap.
c     ------------------------------------------------------------------
      sizeq=sizeq+1
      q(sizeq)=v
      iq(sizeq)=iv
c     ------------------------------------------------------------------
c     Update heap to proper order.
c     ------------------------------------------------------------------
      sq=sizeq
10    tsz=sq/2
      if (tsz.ne.0) then
           if (q(tsz).ge.v) then
                q(sq)=q(tsz)
                iq(sq)=iq(tsz)
                sq=tsz
                goto 10
           endif
      endif
      q(sq)=v
      iq(sq)=iv
c     ------------------------------------------------------------------
      return
      end






      subroutine removq(dimq,iq,iv,q,sizeq,v)
c     ------------------------------------------------------------------
c
c     removq: Remove smallest element (v,iv) from a priority 
c             queue (q,iq).
c
c     ------------------------------------------------------------------
c     Passed input scalar:
c
c         dimq   - array dimension
c
c     ------------------------------------------------------------------
      integer  dimq
c     ------------------------------------------------------------------
c     Passed input/output scalar:
c
c         sizeq  - size of heap
c
c     ------------------------------------------------------------------
      integer sizeq
c     ------------------------------------------------------------------
c     Passed output scalars:
c
c         iv     - smallest element in heap (index)
c         v      - smallest element in heap (value)
c
c     ------------------------------------------------------------------
      integer iv,v
c     ------------------------------------------------------------------
c     Passed input/output arrays:
c
c         iq     - heap (index)
c         q      - heap (value)
c
c     ------------------------------------------------------------------
      integer iq(dimq),q(dimq)
c     ------------------------------------------------------------------
c     Local scalars:
c
c         ivtmp  - tmp smallest element in heap (index)
c         j      - heap counter (2*k)
c         k      - heap counter
c         szqd2  - sizeq/2
c         vtmp   - tmp smallest element in heap (value)
c
c     ------------------------------------------------------------------
      integer ivtmp,j,k,szqd2,vtmp

c     ------------------------------------------------------------------
c     Remove element from heap.
c     ------------------------------------------------------------------
      v=q(1)
      iv=iq(1)
      q(1)=q(sizeq)
      iq(1)=iq(sizeq)
      sizeq=sizeq-1
c     ------------------------------------------------------------------
c     Update heap to proper order.
c     ------------------------------------------------------------------
      k=1
      vtmp=q(k)
      ivtmp=iq(k)
      szqd2=sizeq/2
10    if (k .le. szqd2) then
           j=k+k
           if (j .lt. sizeq) then
                if (q(j) .gt. q(j+1)) j=j+1
           endif
           if (vtmp .gt. q(j)) then
                q(k)=q(j)
                iq(k)=iq(j)
                k=j
                goto 10
           endif
      endif
      q(k)=vtmp
      iq(k)=ivtmp
c     ------------------------------------------------------------------
      return
      end






      subroutine cpi4(dimn,n,x,y)
c     ------------------------------------------------------------------
c
c     cpi4: Copy integer array x of size n into y.  Dimension of
c           arrays is dimn.
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          dimn   - array dimension
c          n      - size of arrays x and y
c
c     ------------------------------------------------------------------
      integer dimn,n
c     ------------------------------------------------------------------
c
c     Passed input array:
c
c          x      - array
c
c     ------------------------------------------------------------------
      integer x(dimn)
c     ------------------------------------------------------------------
c
c     Passed output array:
c
c          y      - array
c
c     ------------------------------------------------------------------
      integer y(dimn)
c     ------------------------------------------------------------------
c
c     Local scalar:
c
c          i      - do loop index
c
c     ------------------------------------------------------------------
      integer i
c     ------------------------------------------------------------------
      do 10 i =1,n
           y(i)=x(i)
10    continue
c     ------------------------------------------------------------------
              return
      end

      program driver
c     ------------------------------------------------------------------
c
c     driver:  Driver for gmis, a GRASP for finding approximate
c              solutions to the maximum independent set problem.
c
c     authors: Mauricio G.C. Resende [mgcr@research.att.com]
c              Thomas A. Feo         [feo@emx.utexas.edu]
c              Stuart H. Smith       [stuart@alk.com]
c
c     ------------------------------------------------------------------
c
c     Dimensions of maximum graph size:
c
c          maxn  - max number of vertices
c          maxa  - max number of edges
c          maxt  - max number of tuples
c          maxtd - max tuple dimension
c
c     ------------------------------------------------------------------
      integer   maxa,maxn,maxt,maxtd
      parameter (maxn=1000,maxa=100000,maxtd=5,maxt=1000)
c     ------------------------------------------------------------------
c
c     Other array dimensions.
c
c          max2a  - max number of vertices
c          maxn2  - max number of edges
c          max2n2 - max number of edges
c
c     ------------------------------------------------------------------
      integer   max2a,maxn2,max2n2
      parameter (max2a=2*maxa,maxn2=maxn*maxn,max2n2=2*maxn2)
c     ------------------------------------------------------------------
c
c     Input/output device numbers passed to iniprt, gmis, outsol:
c
c          in     - input device
c          out    - output device
c
c     ------------------------------------------------------------------
      integer   in,out
      parameter (in=5,out=6)
c     ------------------------------------------------------------------
c
c     Scalars passed to gmis:
c
c          alpha  - GRASP RCL parameter
c          beta   - tuple list size restriction parameter 
c          bitr   - iteration best set was found
c          btup   - tuple on which best set was found
c          errcnd - error condition
c          ittype - iteration summary report
c                   = 0 - silent
c                   = 1 - show improvements only
c                   = 2 - show all iterations
c          look4  - size of independent set sought
c          lstype - type of local search
c          lswhen - indicates when local search is done
c          maxitr - maximum number of GRASP iterations
c          narcs0 - number of arcs in original graph
c          nkset0 - number of vertices in best independent set
c          nlowdg - number of low degree vertices used to make tuples 
c          nnode0 - number of nodes in original graph
c          seed   - seed for pseudo random number generator
c          tupdim - dimension of tuple vector
c          tuplmt - maximum number of tuples allowed
c
c     ------------------------------------------------------------------
      real      alpha,beta
      integer   biter,btup,errcnd,ittype,look4,lstype,lswhen,maxitr,
     +          narcs0,nkset0,nlowdg,nnode0,seed,tupdim,tuplmt
c     ------------------------------------------------------------------
c
c     Arrays passed to gmis:
c
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
c          kset0  - best set of independent nodes, in original graph
c          node1  - vertex 1 of egde
c          node2  - vertex 2 of edge
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
     +          invmap(maxn),iq(maxn2),kset(maxn),kset0(maxn),
     +          node1(maxa),node2(maxa),ptrn(max2a),ptrn0(max2a),
     +          q(maxn2),rcl(maxn),tupind(maxtd),tuple(max2n2),
     +          tuplst(max2n2),vtup(maxtd)


c     ------------------------------------------------------------------
c     Set default values to GRASP parameters.
c     ------------------------------------------------------------------
      call deflt(alpha,beta,errcnd,ittype,lstype,lswhen,maxitr,
     +           nlowdg,seed,tupdim)
c     ------------------------------------------------------------------
c     Read graph.
c     ------------------------------------------------------------------
      call readp(errcnd,in,maxa,maxn,narcs0,nnode0,node1,node2)
      if (errcnd .ge. 100) then
         call errprt(errcnd,out)
           stop
      endif
c     ------------------------------------------------------------------
c
c     Algorithm parameters:
c
c     look4  - size of mis we looking for
c
c     maxitr - maximum number of GRASP iterations
c
c     alpha  - restricted candidate list parameter
c
c     beta   - tuple restrition parameter
c
c     seed   - random number generator seed
c
c     lstype - local search type: = 0  .... no local search
c                                 = 1  .... (1,2)-exchange
c                                 = 2  .... (2,3)-exchange
c
c     lswhen - when is local search is done: = 0, only if phase 1
c                                                 sol'n is better than
c                                                 avg phase 1 sol'n
c                                            = 1, always (if lstype > 0)
c
c     tupdim - tuple dimension
c
c     nlowdg - number of vertices of lowest degree that make up tuples
c
c     ittype - type of iteration summary: = 0, silent
c                                         = 1, print improving solutions
c                                         = 2, print all iterations
c
c     ------------------------------------------------------------------
      look4=nnode0
      maxitr=1000
      alpha=0.1
      beta=0.1
      seed=270001
      lstype=1
      lswhen=0
      tuplmt=1000
      nlowdg=50
      ittype=0
c     ------------------------------------------------------------------
c     Print initial report.
c     ------------------------------------------------------------------
      call iniprt(alpha,beta,in,ittype,look4,lstype,lswhen,max2a,max2n2,
     +            maxa,maxitr,maxn,maxn2,maxt,maxtd,
     +            narcs0,nlowdg,nnode0,out,seed,tupdim,tuplmt)

c     ------------------------------------------------------------------
c     Execute GRASP.
c     ------------------------------------------------------------------
      call gmis(adj,adjn,adjn0,alpha,beta,biter,bkset,btup,deg,deg0,
     +          errcnd,fdmiq,fdmq,fwdmap,incdij,iniad,iniad0,invmap,
     +          iq,ittype,kset,kset0,look4,lstype,lswhen,max2a,max2n2,
     +          maxa,maxitr,maxn,maxn2,maxt,maxtd,narcs0,nkset0,nlowdg,
     +          nnode0,node1,node2,out,ptrn,ptrn0,q,rcl,seed,tupdim,
     +          tupind,tuple,tuplmt,tuplst,vtup)

c     ------------------------------------------------------------------
c     Print out best solution found.
c     ------------------------------------------------------------------
      if (errcnd .lt. 100) then
           call outsol(biter,btup,kset0,maxn,nkset0,out)
      endif
      call errprt(errcnd,out)
c     ------------------------------------------------------------------
      stop
      end





      subroutine deflt(alpha,beta,errcnd,ittype,lstype,lswhen,
     +                 maxitr,nlowdg,seed,tupdim)
c     ------------------------------------------------------------------
c
c     deflt: Sets default values for GRASP parameters.
c
c     ------------------------------------------------------------------
c
c     Passed output scalars:
c
c          alpha  - GRASP RCL parameter
c          beta   - tuple list size restriction parameter 
c          errcnd - error condition
c          ittype - type of iteration summary
c          lstype - type of local search
c          lswhen - indicates when local search is done
c          nlowdg - number of lowest degree nodes in tuples
c          maxitr - maximum number of GRASP iterations
c          seed   - seed for pseudo number generator
c          tupdim - dimension of tuple vector
c
c     ------------------------------------------------------------------
      real     alpha,beta
      integer  errcnd,ittype,lstype,lswhen,nlowdg,maxitr,
     +         seed,tupdim
c     ------------------------------------------------------------------
      alpha = 0.1
      beta = 0.1
      errcnd = 0
      ittype = 0
      lstype = 1
      lswhen = 0
      maxitr = 100000
      nlowdg = 50
      seed = 2700001
      tupdim = 0
c     ------------------------------------------------------------------
      return
      end






      subroutine readp(errcnd,in,maxa,maxn,m,n,node1,node2)
c     ------------------------------------------------------------------
c     readp: Input graph.
c     ------------------------------------------------------------------
c
c     Passed input parameters:
c
c          in     - input device number
c          maxa   - array dimension
c          maxn   - array dimension
c
c     ------------------------------------------------------------------
      integer in,maxa,maxn
c     ------------------------------------------------------------------
c     Passed input scalars:
c
c          n      - number of nodes
c          m      - number of arcs
c
c     ------------------------------------------------------------------
      integer n,m
c     ------------------------------------------------------------------
c     Passed output scalars:
c
c          errcnd - error condition
c
c     ------------------------------------------------------------------
      integer errcnd
c     ------------------------------------------------------------------
c     Passed output arrays:
c
c          node1  - node 1 of arc
c          node2  - node 2 of arc
c
c     ------------------------------------------------------------------
      integer node1(maxa),node2(maxa)
c     ------------------------------------------------------------------
c     Local parameters:
c
c          i     - do loop parameter
c
c     ------------------------------------------------------------------
      integer i
c     ------------------------------------------------------------------
c     Read graph dimensions.
c     ------------------------------------------------------------------
      read (in,*) n,m
c     ------------------------------------------------------------------
c     Check if arrays have been dimensioned correctly.
c     ------------------------------------------------------------------
      if (n.gt.maxn) then
           errcnd = 100
           return
      endif
      if (m.gt.maxa) then
           errcnd = 110
           return
      endif
c     ------------------------------------------------------------------
c     Read arcs.
c     ------------------------------------------------------------------
      do 10 i=1,m
          read (in,*) node1(i),node2(i)
10    continue
c     ------------------------------------------------------------------
      return
      end





      subroutine errprt(errcnd,out)
c     ------------------------------------------------------------------
c
c     errprt: Print out error message cooresponding to error condition
c             given by errcnd (returned from gmis).
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
c                 = 50   - terminated due to iteration limit without
c                          finding set of size look4.
c
c          unsuccessful termination
c
c                 = 100  - maxn too small
c                 = 110  - maxa too small
c                 = 120  - maxt too small
c                 = 130  - maxtd too small
c                 = 140  - alpha too small
c                 = 141  - alpha too big
c                 = 150  - beta too small
c                 = 151  - beta too big
c                 = 160  - look4 too big
c                 = 170  - lstype invalid
c                 = 171  - lswhen invalid
c                 = 180  - maxitr too small
c                 = 190  - nlowdg too small
c                 = 200  - seed too small
c                 = 201  - seed too big
c                 = 210  - tupdim too small
c                 = 211  - tupdim too big
c
c          out    - Fortran output device
c
c     ------------------------------------------------------------------
      integer   errcnd,out
c     ------------------------------------------------------------------
      if (errcnd .eq. 0) then
           write(out,10)
10         format(' gmis: terminated successfully')
           return
      endif
      if (errcnd .eq. 50) then
           write(out,20)
20         format(' gmis: iteration limit reached')
           return
      endif
      if (errcnd .eq. 100) then
           write(out,30)
30         format(' gmis: array too small - increase maxn')
           return
      endif
      if (errcnd .eq. 110) then
           write(out,40)
40         format(' gmis: array too small - increase maxa')
           return
      endif
      if (errcnd .eq. 120) then
           write(out,50)
50         format(' gmis: array too small - increase maxt')
           return
      endif
      if (errcnd .eq. 130) then
           write(out,60)
60         format(' gmis: array too small - increase maxtd')
           return
      endif
      if (errcnd .eq. 140) then
           write(out,70)
70         format(' gmis: error in value of alpha - too small')
           return
      endif
      if (errcnd .eq. 141) then
           write(out,80)
80         format(' gmis: error in value of alpha - too big')
           return
      endif
      if (errcnd .eq. 150) then
           write(out,90)
90         format(' gmis: error in value of beta - too small')
           return
      endif
      if (errcnd .eq. 151) then
           write(out,100)
100        format(' gmis: error in value of beta - too big')
           return
      endif
      if (errcnd .eq. 160) then
           write(out,110)
110        format(' gmis: error in value of look4 - too big')
           return
      endif
      if (errcnd .eq. 170) then
           write(out,120)
120        format(' gmis: error in value of lstype - invalid')
           return
      endif
      if (errcnd .eq. 171) then
           write(out,130)
130        format(' gmis: error in value of lswhen - invalid')
           return
      endif
      if (errcnd .eq. 180) then
           write(out,140)
140        format(' gmis: error in value of maxitr - too small')
           return
      endif
      if (errcnd .eq. 190) then
           write(out,150)
150        format(' gmis: error in value of nlowdg - too small')
           return
      endif
      if (errcnd .eq. 200) then
           write(out,160)
160        format(' gmis: error in value of seed - too small')
           return
      endif
      if (errcnd .eq. 201) then
           write(out,170)
170        format(' gmis: error in value of seed - too big')
           return
      endif
      if (errcnd .eq. 210) then
           write(out,180)
180        format(' gmis: error in value of tupdim - too small')
           return
      endif
      if (errcnd .eq. 211) then
           write(out,190)
190        format(' gmis: error in value of tupdim - too big')
           return
      endif
c     ------------------------------------------------------------------
      end





      subroutine iniprt(alpha,beta,in,ittype,look4,lstype,lswhen,
     +                  max2a,max2n2,maxa,maxitr,maxn,maxn2,maxt,maxtd,     
     +                  narcs0,nlowdg,nnode0,out,seed,tupdim,tuplmt)
c     ------------------------------------------------------------------
c
c     iniprt:  Prints initial report with code parameter settings and 
c              graph statistics.
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          alpha  - GRASP RCL parameter
c          beta   - tuple list size restriction parameter 
c          in     - input device
c          ittype - iteration summary type
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
c          out    - output device
c          seed   - seed for pseudo random number generator
c          tupdim - dimension of tuple vector
c          tuplmt - maximum number of tuples allowed
c
c     ------------------------------------------------------------------
      real     alpha,beta
      integer  in,ittype,look4,lstype,lswhen,max2a,max2n2,maxa,
     +         maxitr,maxn,maxn2,maxt,maxtd,narcs0,nlowdg,nnode0,
     +         out,seed,tupdim,tuplmt
c     ------------------------------------------------------------------
c
c     Local scalar:
c
c          d      - graph density
c          totspc - total array memory allocated
c
c     ------------------------------------------------------------------
      real d,totspc
c     ------------------------------------------------------------------
      write(out,10)
10    format(' --------gmis: GRASP for Max Independent Set------------')
      write(out,20)
20    format('  ')
      write(out,30)
30    format(' i/o device numbers::: ')
      write(out,40) in
40    format('      in             : ',i10)
      write(out,50) out
50    format('      out            : ',i10)
      write(out,20)
      write(out,60)
60    format(' memory allocation:::: ')
      write(out,70) maxn
70    format('      max nodes      : ',i10)
      write(out,80) maxa
80    format('      max arcs       : ',i10)
      write(out,90) maxt
90    format('      max tuples     : ',i10)
      write(out,100) maxtd
100   format('      max tuple dim  : ',i10)
      totspc=4*(2*(maxa+maxt+maxtd+max2n2)+3*maxn2+
     +       4*max2a+10*maxn)/1.0e6
      write(out,110) totspc
110   format('      array size (Mb): ',f10.1)
      write(out,20)
      write(out,120)
120   format(' algorithm control:::: ')
      write(out,130) seed
130   format('       seed          : ',i10)
      write(out,140) alpha
140   format('       alpha         : ',f10.2)
      write(out,150) beta
150   format('       beta          : ',f10.2)
      write(out,160) look4
160   format('       look4         : ',i10)
      write(out,170) maxitr
170   format('       maxitr        : ',i10)
      write(out,20)
      write(out,180)
180   format(' local search control: ')
      write(out,190) lstype
190   format('       lstype        : ',i10)
      write(out,200) lswhen
200   format('       lswhen        : ',i10)
      write(out,20)
      write(out,210)
210   format(' tuple control:::::::: ')
      write(out,220) nlowdg
220   format('       nlowdg        : ',i10)
      write(out,230) tupdim
230   format('       tupdim        : ',i10)
      write(out,240) tuplmt
240   format('       tuplmt        : ',i10)
      write(out,20)
      write(out,250)
250   format(' output control::::::: ')
      if (ittype.eq.0) then
           write(out,260)
260        format('       itr summary   :     silent')
      endif
      if (ittype.eq.1) then
           write(out,270)
270        format('       itr summary   :improve itr')
      endif
      if (ittype.eq.2) then
           write(out,280)
280        format('       itr summary   :  all iters')
      endif
      write(out,20)
      write(out,290)
290   format(' graph dimensions::::: ')
      d=float(narcs0)/float(nnode0*(nnode0-1)/2)
      write(out,300) nnode0
300   format('       nodes         : ',i10)
      write(out,310) narcs0
310   format('       arcs          : ',i10)
      write(out,320) d
320   format('       density       : ',f10.3)
      write(out,20)
c     ------------------------------------------------------------------
      return
      end





      subroutine outsol(biter,btup,kset,maxn,nkset,out)
c     ------------------------------------------------------------------
c
c     outsol: Print out size and vertices of independent set found. 
c
c     ------------------------------------------------------------------
c
c     Passed input scalars:
c
c          biter  - iteration best solution was found
c          btup   - tuple best solution was found
c          maxn   - dimension of kset array
c          nkset  - number of vertices in independent set
c          out    - output device number
c
c     ------------------------------------------------------------------
      integer  biter,btup,maxn,nkset,out
c     ------------------------------------------------------------------
c
c     Passed inpout array:
c
c          kset   - array if vertices of independent set
c
c     ------------------------------------------------------------------
      integer  kset(maxn)
c     ------------------------------------------------------------------
c
c     Local scalars:
c
c          i      - do loop index
c
c     ------------------------------------------------------------------
      integer  i
c     ------------------------------------------------------------------
c     Write size of largest independent set found.
c     ------------------------------------------------------------------
      write(out,10)
10    format(
     +' output--------------------------------------------------------')
      write(out,20) nkset
20    format(' Size of indep set : ',i10)
      write(out,30) btup
30    format(' Found on tuple    : ',i10)
      write(out,40) biter
40    format(' Iteration         : ',i10)
c     ------------------------------------------------------------------
c     Write vertices of largest independent set found.
c     ------------------------------------------------------------------
      write(out,50) (kset(i),i=1,nkset)
50    format(' Independent set   : ',
     +       4i10/(19x,': ',4i10))
      write(out,60)
60    format(
     +' --------------------------------------------------------------')

c     ------------------------------------------------------------------
      return
      end

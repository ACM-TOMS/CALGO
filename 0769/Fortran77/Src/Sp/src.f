c     ------------------------------------------------------------------
c     A Greedy Randomized Adaptive Search Procedure (GRASP) for the
c     Quadratic Assignment Problem (QAP).
c
c     Authors: M.G.C. Resende (AT&T Bell Laboratories)
c                [mgcr@research.att.com]
c              L.S. Pitsoulis (University of Florida)
c                [leonidas@deming.ise.ufl.edu]
c              P.M. Pardalos (University of Florida)
c                [pardalos@deming.ise.ufl.edu]
c
c     ver. 1.0 
c     ------------------------------------------------------------------
c     This file includes the following Fortran subroutines and
c     functions:
c
c          gqaps   - control subroutine for GRASP for QAP algorithm
c          srtcst  - sorts cost 
c          stage1  - stage 1 of GRASP construction phase
c          stage2  - stage 2 of GRASP construction phase
c          local   - 2-exchange local search for QAP 
c          mkbseq  - makes permutation vector b = (1,2,...,n)
c          mkassg  - changes the permutations for a new assignment
c          insrtq  - insert element into heap for sorting
c          removq  - remove element from heap
c          evalij  - evaluates the cost effect of swapping i and j
c          randp   - random number generator function
c
c     ------------------------------------------------------------------

      subroutine gqaps(n,n2,maxitr,alpha,beta,look4,seed,f,d,a,b,
     &                 optb,srtf,srtif,srtd,srtid,srtc,srtic,indexd,
     &                 indexf,cost,fdind,opta,bestv,iter,cp)
c     ------------------------------------------------------------------
c     gqaps: Subroutine for finding an approximate solution of a
c            sparse symmetric quadratic assignment problem.
c
c     ------------------------------------------------------------------
c     Parameters:
c
c          infty - a large integer
c
c     ------------------------------------------------------------------
      integer   infty
      parameter (infty = 2147483647)
c     ------------------------------------------------------------------
c     Passed input scalars:
c
c         n      - dimension of qap problem
c         n2     - n * n
c         maxitr - maximum number of GRASP iterations
c         look4  - if permutation of cost look4 or less is found gqapd 
c                  returns that permutation
c         bestv  - cost of best assignment found
c         iter   - number of GRASP iterations taken
c         alpha  - phase 1 parameter
c         beta   - phase 1 parameter
c
c     ------------------------------------------------------------------
      integer   n,n2,maxitr,look4
      real      alpha,beta
c     ------------------------------------------------------------------
c     Passed input/output scalar:
c
c          seed   - random number generator seed
c
c     ------------------------------------------------------------------
      integer   seed
c     ------------------------------------------------------------------
c     Passed output scalars: 
c
c     bestv       - cost of best assignment found
c     iter        - number of GRASP iterations taken
c
c     ------------------------------------------------------------------
      integer   bestv,iter
c     ------------------------------------------------------------------
c     Passed input arrays: 
c
c          f      - flow matrix stored as a 1-dimensional array,
c                   row by row
c          d      - distance matrix stored as a 1-dimensional array,
c                   row by row
c
c     ------------------------------------------------------------------
      integer f(n2),d(n2)
c     ------------------------------------------------------------------
c     Passed work arrays: 
c
c          a      - permutation vector
c          b      - permutation vector
c          cp     - sparse calculation matrix
c          srtf   - sorted F values
c          srtif  - sorted F values (indices)
c          srtd   - sorted D values
c          srtid  - sorted D values (indices)
c          srtc   - sorted cost values
c          srtic  - sorted cost values (indices)
c          indexf - indices of facilities in unsorted cost matrix
c          indexd - indices of locations in unsorted cost matrix
c          cost   - sorted cost matrix
c          fdind  - indices of sorted cost matrix
c
c     ------------------------------------------------------------------
      integer   a(n),b(n),srtf(n2),srtif(n2),srtd(n2),srtid(n2),
     &          srtc(n2),srtic(n2),indexf(n2),indexd(n2),cost(n2),
     &          fdind(n2),cp(n2)
c     ------------------------------------------------------------------
c     Passed output arrays: 
c
c          opta   - best permutation vector
c          optb   - best permutation vector
c
c     ------------------------------------------------------------------
      integer   opta(n),optb(n)
c     ------------------------------------------------------------------
c     Local scalars:
c
c         i      - facility index
c         j      - facility index
c         k      - location index
c         l      - location index
c         objv   - cost of permutation
c    
c     ------------------------------------------------------------------
      integer   i,j,k,l,objv
c     ------------------------------------------------------------------
c     Initialize cost of best assignment found to large number.
c     ------------------------------------------------------------------
      bestv = infty
c     ------------------------------------------------------------------
c     Sort the cost = f(i,j) * d(k,l) in increasing order to be used 
c     by stage1.
c     ------------------------------------------------------------------
      call srtcst(n,n2,f,d,beta,srtf,srtif,srtd,srtid,srtc,
     &            srtic,indexd,indexf,cost,fdind)
c     ------------------------------------------------------------------
c     Do GRASP iterations.
c     ------------------------------------------------------------------
      do 10 iter = 1,maxitr
c          -------------------------------------------------------------
c          Stage 1 of construction phase.
c          -------------------------------------------------------------
           call stage1(n,n2,i,j,k,l,seed,alpha,beta,objv,indexd,indexf,
     &                 fdind,cost,a,b)
c          -------------------------------------------------------------
c          Stage 2 of construction phase.
c          -------------------------------------------------------------
           call stage2(n,n2,f,d,alpha,srtc,srtic,a,b,i,j,k,l,seed,
     &                 objv,cp)
c          -------------------------------------------------------------
c          Local search phase.
c          -------------------------------------------------------------
           call local(n,n2,objv,f,d,a,b)
c          -------------------------------------------------------------
c          If cost of current permutation is best so far, save 
c          permutation and corresponding objective value.
c          -------------------------------------------------------------
           if (objv .lt. bestv) then
               call savsol(n,objv,bestv,a,opta)
c              ----------------------------------------------------------
c              If cost of assignment is as good as wished, return.
c              ----------------------------------------------------------
               if (bestv .le. look4) return
           endif
10    continue
c     -----------------------------------------------------------------
c     Adjust iteration counter.
c     -----------------------------------------------------------------
      iter = maxitr
c     ------------------------------------------------------------------
      return

c     end of grasps
      end



      subroutine srtcst(n,n2,f,d,beta,srtf,srtif,srtd,srtid,srtc,
     &                  srtic,indexd,indexf,cost,fdind)
c     ------------------------------------------------------------------
c     srtcst:  Sorts cost = f(i,j)*d(k,l) in increasing order.
c     ------------------------------------------------------------------
c     Passed input scalars:
c 
c         n      - qap dimension
c         n2     - n * n
c         beta   - construction phase parameter
c
c     ------------------------------------------------------------------
      integer n,n2
      real    beta
c     ------------------------------------------------------------------
c     Passed input arrays:
c
c         f      - flow matrix
c         d      - distance matrix
c
c     ------------------------------------------------------------------
      integer f(n2),d(n2)
c     ------------------------------------------------------------------
c     Passed work arrays: 
c
c         srtf   - sorted flow matrix (values)
c         srtif  - sorted flow matrix (indices)
c         srtd   - sorted distance matrix (values)
c         srtid  - sorted distance matrix (indices)
c         srtc   - sorted cost matrix (values)
c         srtic  - sorted cost matrix (indices)
c
c     ------------------------------------------------------------------
      integer srtf(n2),srtif(n2),srtd(n2),srtid(n2),srtc(n2),srtic(n2)
c     ------------------------------------------------------------------
c     Passed output arrays:
c
c         indexd - indices of locations in unsorted cost matrix
c         indexf - indices of facilities in unsorted cost matrix
c         cost   - sorted cost matrix
c         fdind  - indices of sorted cost matrix
c
c     ------------------------------------------------------------------
      integer indexd(n2),indexf(n2),cost(n2),fdind(n2)
c     ------------------------------------------------------------------
c     Local scalars: 
c
c         index  - index
c         sizec  - number of elements in cost heap
c         sized  - number of elements in distance heap
c         sizef  - number of elements in flow heap
c         dv     - distance value
c         fv     - flow value
c         dind   - distance index
c         find   - flow index
c         nbeta  - number of candidates
c         i      - do loop index
c         j      - do loop index
c     
c     ------------------------------------------------------------------
      integer index,sizec,sized,sizef,dv,fv,dind,find,nbeta,i,j
c     ------------------------------------------------------------------
c     Sort D in increasing order
c          F in decreasing order (-F in increasing order).
c
c     Keep only (n*n - n)*beta best in each sorting.
c
c     Initialize cardinalities of sorted sets of elements of D,F
c     and cost.
c     ------------------------------------------------------------------
      sized = 0
      sizef = 0
      sizec = 0
c     ------------------------------------------------------------------
c     Insert all non-diagonal elements of D into D-priority heap.
c     Insert all non-diagonal elements of -F into F-priority heap.
c     ------------------------------------------------------------------
      index = 0
      do 20 i = 1,n
           do 10 j = 1,n
              index = index + 1
              if (i .ne. j) then 
                 call insrtq(n2,d(index),index,sized,srtd,srtid)
                 call insrtq(n2,-f(index),index,sizef,srtf,srtif)
               endif
10          continue
20    continue
c     ------------------------------------------------------------------
c     Compute size of sorted sets.
c     ------------------------------------------------------------------
      nbeta = beta*(n*n - n)
c     ------------------------------------------------------------------
c     Remove the nbeta smallest D elements from D-priority heap, and 
c     the nbeta smallest -F elements from F-priority heap.
c     ------------------------------------------------------------------
      do 30 i = 1,nbeta
           call removq(n2,dv,dind,sized,srtd,srtid)
           call removq(n2,fv,find,sizef,srtf,srtif)
c          -------------------------------------------------------------
c          Cost is product of sorted flow and distance.
c          -------------------------------------------------------------
           cost(i)   = -dv*fv
           indexd(i) = dind
           indexf(i) = find
c          -------------------------------------------------------------
c          Insert cost into cost-priority heap.
c          -------------------------------------------------------------
           call insrtq(n2,cost(i),i,sizec,srtc,srtic)
30    continue
c     ------------------------------------------------------------------
c     Remove nbeta sorted cost elements from cost-priority heap.
c     ------------------------------------------------------------------
      do 40 i = 1,nbeta
           call removq(n2,cost(i),fdind(i),sizec,srtc,srtic)
40    continue
c     ------------------------------------------------------------------
      return

c     end of srtcst
      end



      subroutine stage1(n,n2,i,j,k,l,seed,alpha,beta,objv,indexd,indexf,
     &                  fdind,cost,a,b)
c     ------------------------------------------------------------------
c     stage1:  Builds the initial 2 assignments for the GRASP
c              construction phase (facility i to site k and
c              facility j to site l).
c     ------------------------------------------------------------------
c     Passed input scalars:
c
c         n      - qap dimension
c         n2     - n * n
c         alpha  - construction phase parameter
c         beta   - construction phase parameter
c
c     ------------------------------------------------------------------
      integer n,n2
      real    alpha,beta
c     ------------------------------------------------------------------
c     Passed input/output scalar:
c
c         seed   - random number generator seed      
c
c     ------------------------------------------------------------------
      integer seed
c     ------------------------------------------------------------------
c     Passed output scalars:
c
c         i      - facility index
c         j      - facility index
c         k      - location index
c         l      - location index
c         objv   - cost of initial 2 assignments (returned)    
c
c     ------------------------------------------------------------------
      integer i,j,k,l,objv
c     ------------------------------------------------------------------
c     Passed input arrays:
c
c         indexd - indices of locations in unsorted cost matrix
c         indexf - indices of facilities in unsorted cost matrix
c         fdind  - indices of sorted cost matrix
c         cost   - cost of assignment
c
c     ------------------------------------------------------------------
      integer indexd(n2),indexf(n2),fdind(n2),cost(n2)
c     ------------------------------------------------------------------
c     Passed output arrays:
c
c         a      - permutation array
c         b      - permutation array
c
c     ------------------------------------------------------------------
      integer a(n),b(n)
c     ------------------------------------------------------------------
c     Local scalars and functions:
c
c         nselct  - index of randomly selected element
c         dind    - distance index
c         find    - flow index
c         high    - upper bound of selection range
c         ii      - loop index
c         tmp     - temporary integer variable
c         randp   - random number generator function
c         xrand   - dummy probability
c
c     ------------------------------------------------------------------
      integer nselct,dind,find,high,ii,tmp
      real    randp,xrand
c     ------------------------------------------------------------------
c     Initialize permutations.
c     ------------------------------------------------------------------
      do 10 ii = 1,n
           a(ii) = ii
           b(ii) = ii
10    continue
c     ------------------------------------------------------------------
c     Select element, at random, from the best (n*n - n)*alpha cost 
c     elements.
c     ------------------------------------------------------------------
      xrand  = randp(seed)
      high   = alpha*beta*(n*n - n)
      nselct = 1 + seed/(2147483647/high)
c     ------------------------------------------------------------------
c     Initial assignment is facility i to location k
c                           facility j to location l.
c     
c     Cost of initial assignment is f(i,j) * d(k,l).
c     ------------------------------------------------------------------
      dind = indexd(fdind(nselct))
      find = indexf(fdind(nselct))
      i    = (find - 1)/n + 1
      j    = find - (i - 1)*n
      k    = (dind - 1)/n + 1
      l    = dind - (k - 1)*n
      objv = cost(nselct)
c     ------------------------------------------------------------------
c     Make initial assignments to permutation arrays:
c
c     Assign facility i to location k.
c     ------------------------------------------------------------------
      a(1) = i
      a(i) = 1
      b(1) = k
      b(k) = 1
c     ------------------------------------------------------------------
c     Assign facility j to location l.
c     ------------------------------------------------------------------
      do 20 ii = 1,n
           if (a(ii) .eq. j) then
                tmp   = a(2)
                a(2)  = j
                a(ii) = tmp
                goto 30
           endif
20    continue
30    do 40 ii = 1,n
           if (b(ii) .eq. l) then
                tmp   = b(2)
                b(2)  = l
                b(ii) = tmp
                goto 50
           endif
40    continue
c     ------------------------------------------------------------------
50    return

c     end of stage1
      end







      subroutine stage2(n,n2,f,d,alpha,srtc,srtic,a,b,i,j,k,l,seed,
     &                  objv,cp)
c     ------------------------------------------------------------------
c     build:  builds a randomized greedy permutation starting from
c             the assignments made in stage1.  Permutation is returned
c             in array a(*).
c     ------------------------------------------------------------------
c     Passed input scalars:
c
c         n      - problem dimension
c         n2     - n * n
c         i      - facility index
c         j      - location index
c         k      - facility index
c         l      - location index
c         alpha  - construction phase parameter
c
c     ------------------------------------------------------------------
      integer n,n2,i,j,k,l
      real    alpha
c     ------------------------------------------------------------------
c     Passed input/output scalars:
c
c         seed   - random number generator seed
c         objv   - cost of assignment
c     
c     ------------------------------------------------------------------
      integer seed,objv
c     ------------------------------------------------------------------
c     Passed input arrays:
c
c         f      - flow matrix
c         d      - distance matrix
c
c     ------------------------------------------------------------------
      integer f(n2),d(n2)
c     ------------------------------------------------------------------
c     Passed work arrays:
c
c         srtc   - sorted cost matrix (values)
c         srtic  - sorted cost matrix (indices)
c         cp     - zero cost assignment matrix
c
c     ------------------------------------------------------------------
      integer srtc(n2),srtic(n2),cp(n2)
c     ------------------------------------------------------------------
c     Passed input/output arrays:
c
c         a      - permutation array
c         b      - permutation array
c
c     ------------------------------------------------------------------
      integer a(n2),b(n2)
c     ------------------------------------------------------------------
c     Local scalars and functions
c
c         high   - upper bound of selection range
c         assign - do loop counter of assignments  
c         cost   - assignment cost
c         sizec  - number of cost elemnets in cost heap
c         nselct - selected index
c         fdind  - index of f d product 
c         scp    - zero cost counter variable
c         oldscp - previous iteration scp
c         assz   - zero cost assignments made
c         tl     - temporary index of l
c         tk     - temporary index of k
c         akm1tn - (a(k) - 1)*n
c         blm1tn - (b(l) - 1)*n
c         anm1tn - (a(n) - 1)*n
c         bnm1tn - (b(n) - 1)*n
c         randp  - random number generator function
c         xrand  - random number
c
c     ----------------------------------------------------------------
      integer high,assign,cost,sizec,nselct,fdind,
     &        akm1tn,blm1tn,anm1tn,bnm1tn,scp,oldscp,assz,tl,tk
      real    randp,xrand
c     ----------------------------------------------------------------
c     Initial phase of stage2 where the assignments resulting
c     in zero overall cost are inspected and kept. If the number 
c     of the available zero cost assignment pairs is 0 go to the 
c     second phase.
c     ----------------------------------------------------------------
c--->      objv = objv + objv
      assz = 2
      scp  = 0
c     ----------------------------------------------------------------
c     Construct the array cp of size scp that contains the indices 
c     of the zero cost assignment pairs, with respect to the initial 
c     two assignments of stage1.
c     ----------------------------------------------------------------
      do 30 k = assz + 1,n
         akm1tn = (a(k) - 1)*n
         do 20 l = assz + 1,n
            blm1tn = (b(l) - 1)*n
            cost = 0
            do 10 i = 1,2
               cost = cost + f(akm1tn + a(i))*d(blm1tn + b(i))
10          continue
            if (cost .eq. 0) then
               scp     = scp + 1
               cp(scp) = akm1tn + b(l)
            endif
20       continue
30    continue
      goto 60
c     ---------------------------------------------------------------
c     For the 4,5,...,assz zero cost assignment pairs, we only need 
c     to examine the previously constructed cp array. For every new 
c     zero cost assignment pair the size of the cp array reduces. 
c     ---------------------------------------------------------------
40    oldscp = scp
      scp    = 0
      do 50 i = 1,oldscp
         tk = (cp(i) - 1)/n + 1
         tl = cp(i) - (tk - 1)*n
         if (tk .eq. a(assz) .or. tl .eq. b(assz)) goto 50
         if (f((tk - 1)*n + a(assz)) .eq. 0 .or.
     &       d((tl - 1)*n + b(assz)) .eq. 0) then
            scp     = scp + 1
            cp(scp) = cp(i)
         endif
50    continue
c     --------------------------------------------------------------
c     From the available scp zero cost assignment pairs select one 
c     at random to be added in the partialy constructed solution.
c     --------------------------------------------------------------
60    high = scp
      if (scp .eq. 0) then 
         goto 70
      endif
      xrand  = randp(seed)
      nselct = 1 + seed/(2147483647/scp)
      fdind  = cp(nselct)
      call mkassg(n,assz + 1,fdind,a,b)
      assz   = assz + 1
      goto 40
c     --------------------------------------------------------------
c     Assignments assz + 1,4,..,n - 1 are made.
c     --------------------------------------------------------------
70    do 120 assign = assz + 1,n - 1
c        -----------------------------------------------------------
c        For all pairs not assigned yet, compute costs of all
c        possible assignments, w.r.t. to already-made assignments.
c        -----------------------------------------------------------
         sizec = 0
         do 100 k = assign,n
            akm1tn = (a(k) - 1)*n
            do 90 l = assign,n
               blm1tn = (b(l) - 1)*n
               cost = 0
               do 80 i = 1,assign - 1
c                 --------------------------------------------------
c                 Facility a(i) already assigned to location b(i):
c                 cost of assigning facility a(k) to location b(l)
c                 relative to assignment of  facility a(i) to 
c                 location b(i).
c                 --------------------------------------------------
                  cost = cost + f(akm1tn + a(i))*d(blm1tn + b(i))
80             continue
c              -----------------------------------------------------
c              Insert cost element into cost-priority heap for 
c              sorting.
c              -----------------------------------------------------
               call insrtq(n2,cost,akm1tn + b(l),sizec,srtc,srtic)
90           continue
100        continue
c        -----------------------------------------------------------
c        Select at random from the best alpha*sizec assignments.
c        -----------------------------------------------------------
         xrand  = randp(seed)
         high   = alpha*sizec
         nselct = 1 + seed/(2147483647/high)
         do 110 i = 1,nselct
            call removq(n2,cost,fdind,sizec,srtc,srtic)
110      continue
c        -----------------------------------------------------------
c        Make assignment.
c        -----------------------------------------------------------
c---->   objv = objv + 2*cost
         objv = objv + cost
         call mkassg(n,assign,fdind,a,b)
120   continue
      anm1tn = (a(n) - 1)*n
      bnm1tn = (b(n) - 1)*n
      do 130 i = 1,n - 1
c---->   objv = objv + 2*f(anm1tn + a(i))*d(bnm1tn + b(i))
         objv = objv + f(anm1tn + a(i))*d(bnm1tn + b(i))
130   continue
      objv = objv + objv
c     -----------------------------------------------------------------
      return

c     end of stage2
      end


      


      subroutine mkassg(n,assign,fdind,a,b)
c     -----------------------------------------------------------------
c     mkassg: given a facility to be assigned to a location as 
c             indicated by fdind, change the permutation vectors a and
c             b to reflect this assignment.
c     -----------------------------------------------------------------
c     Passed input scalars:
c
c               n       - problem dimension
c               assign  - assignment in consideration
c               fdind   - index of f d product
c
c     -----------------------------------------------------------------
      integer n,assign,fdind
c     -----------------------------------------------------------------
c     Passed input/output arrays: 
c
c               a       - permutation array
c               b       - permutation array
c
c     -----------------------------------------------------------------
      integer a(n), b(n)
c     -----------------------------------------------------------------      
c     Local scalars:
c         i      - facility index
c         j      - location index
c         k      - facility index
c         l      - location index
c         kinv   - index of k in inverted permutation 
c         linv   - index of l in inverted permutation
c         tmp    - temporary integer variable
c
c     -----------------------------------------------------------------     
      integer i,j,k,l,linv,kinv,tmp
c     -----------------------------------------------------------------
      kinv = (fdind - 1)/n + 1
      linv = fdind - (kinv - 1)*n
      do 10 i = assign,n
         if (a(i) .eq. kinv) then
            k = i
            goto 20
         endif  
10    continue
20    do 30 j = assign,n
         if (b(j) .eq. linv) then
            l = j
            goto 40
         endif  
30    continue
40    tmp       = a(assign)
      a(assign) = a(k)
      a(k)      = tmp
      tmp       = b(assign)
      b(assign) = b(l)
      b(l)      = tmp
c     -----------------------------------------------------------------
      return

c     end of mkassg
      end

      subroutine savsol(n,objv,bestv,a,opta)
c     ------------------------------------------------------------------
c     savsol:  Saves current best solution.
c     ------------------------------------------------------------------
c     Passed scalars:
c
c         n      - problem dimension
c         objv   - objective function value
c
c     ------------------------------------------------------------------
      integer n,objv
c     ------------------------------------------------------------------
c     Passed output scalars:
c
c         bestv  - best objective function value so far
c
c     ------------------------------------------------------------------
      integer bestv
c     ------------------------------------------------------------------
c     Passed input array:
c
c         a      - permutation array
c
c     ------------------------------------------------------------------
      integer a(n)
c     ------------------------------------------------------------------
c     Passed output array:
c     
c         opta   - array of best permutation so far
c
c     ------------------------------------------------------------------
      integer opta(n)
c     ------------------------------------------------------------------
c     Local scalar:
c
c         i      - loop index
c
c     ------------------------------------------------------------------
      integer i
c     ------------------------------------------------------------------
      do 10 i = 1,n
           opta(i) = a(i)
10    continue
      bestv = objv
c     ------------------------------------------------------------------
      return

c     end of savsol
      end



      subroutine local(n,n2,objv,f,d,a,b)
c     ------------------------------------------------------------------
c     local: Local 2-exchange on permutation array a. 
c            Return improved permutation array a and objv.
c     ------------------------------------------------------------------
c     Passed input scalars:
c
c         n      - problem dimension
c         n2     - n * n
c
c     ------------------------------------------------------------------
      integer n,n2
c     ------------------------------------------------------------------
c     Passed input/output scalars:
c
c         objv   - objective function value
c
c     ------------------------------------------------------------------
      integer objv
c     ------------------------------------------------------------------
c     Passed input arrays:
c
c         f      - flow matrix
c         d      - distance matrix
c
c     ------------------------------------------------------------------
      integer f(n2),d(n2)
c     ------------------------------------------------------------------
c     Passed input/output arrays:
c
c         a      - permutation array
c         b      - permutation array
c
c     ------------------------------------------------------------------
      integer a(n),b(n)
c     ------------------------------------------------------------------
c     Local scalars:
c
c         i      - loop index
c         j      - loop index
c         temp   - temp scalar used to swap a(i) and a(j)
c         xgain  - gain from switch
c         improv - objective function improvement
c
c     ------------------------------------------------------------------
      integer i,j,temp,xgain
      logical improv

c     ------------------------------------------------------------------
c     Make array b(*) = (1,2,3,...,n) for local search.
c     ------------------------------------------------------------------
      call mkbseq(n,a,b)
c     ------------------------------------------------------------------
c     Attempt to switch all pairs in permutation array a.
c     ------------------------------------------------------------------
10    improv = .false.
      do 30 i = 1,n - 1
           do 20 j = i + 1,n
c               --------------------------------------------------------
c               Evaluate cost difference by adopting switch of a(i) and 
c               a(j).
c               --------------------------------------------------------
                call evalij(n,n2,i,j,xgain,f,d,a)
c               --------------------------------------------------------
c               If switch improves cost, adopt it.
c               --------------------------------------------------------
                if (xgain .gt. 0) then
                     temp   = a(i)
                     a(i)   = a(j)
                     a(j)   = temp
                     objv   = objv - xgain
                     improv = .true.
                endif
20         continue
30    continue
c     ------------------------------------------------------------------
c     If no switch improves cost (improv = .false.) return; else repeat.
c     ------------------------------------------------------------------
      if (improv) goto 10
c     ------------------------------------------------------------------
      return

c     end of local
      end


      subroutine mkbseq(n,a,b)
c     ------------------------------------------------------------------
c     mkbseq: Change permutation arrays a and b to make b = (1,2,...,n).
c     ------------------------------------------------------------------
c     Passed input scalar:
c
c         n      - QAP dimension
c
c     ------------------------------------------------------------------
      integer n
c     ------------------------------------------------------------------
c     Passed input/output arrays:
c
c         a      - permutation array
c         b      - permutation array
c
c     ------------------------------------------------------------------
      integer a(n),b(n)
c     ------------------------------------------------------------------
c     Local scalars:
c
c         i      - loop index
c         j      - loop index
c         tmp    - temporary integer variable
c
c     ------------------------------------------------------------------
      integer i,j,tmp
c     ------------------------------------------------------------------
      do 20 i = 1,n - 1
           do 10 j = i + 1,n
                if (b(j) .eq. i) then
                     b(j) = b(i)
                     b(i) = i
                     tmp  = a(i)
                     a(i) = a(j)
                     a(j) = tmp
                     goto 20
                endif
10         continue
20    continue
c     ------------------------------------------------------------------
      return

c     end of mkbseq
      end


      subroutine insrtq(n2,v,iv,sizeq,q,iq)
c     ------------------------------------------------------------------
c     insrtq: Insert an element (v,iv) into a queue (q,iq).
c     ------------------------------------------------------------------
c     Passed input scalars:
c
c         n2     - n * n
c         v      - heap element (value)
c         iv     - heap element (index)
c
c     ------------------------------------------------------------------
      integer n2,v,iv
c     ------------------------------------------------------------------
c     Passed input/output scalar:
c
c         sizeq  - size of heap
c
c     ------------------------------------------------------------------
      integer sizeq
c     ------------------------------------------------------------------
c     Passed input/output arrays:
c
c         q      - heap (value)
c         iq     - heap (index)
c
c     ------------------------------------------------------------------
      integer q(n2),iq(n2)
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
      sizeq     = sizeq + 1
      q(sizeq)  = v
      iq(sizeq) = iv
c     ------------------------------------------------------------------
c     Update heap to proper order.
c     ------------------------------------------------------------------
      sq  = sizeq
      v   = q(sq)
      iv  = iq(sq)
10    tsz = sq/2
      if (tsz .ne. 0) then
           if (q(tsz) .gt. v) then
                q(sq)  = q(tsz)
                iq(sq) = iq(tsz)
                sq     = tsz
                goto 10
           endif
      endif
      q(sq)  = v
      iq(sq) = iv
c     ------------------------------------------------------------------
      return

c     end of insrtq
      end


      subroutine removq(n2,v,iv,sizeq,q,iq)
c     ------------------------------------------------------------------
c     removq: Remove smallest element (v,iv) from a priority 
c             queue (q,iq).
c     ------------------------------------------------------------------
c     Passed input scalar:
c
c         n2     - n * n
c
c     ------------------------------------------------------------------
      integer n2
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
c         v      - smallest element in heap (value)
c         iv     - smallest element in heap (index)
c
c     ------------------------------------------------------------------
      integer v,iv
c     ------------------------------------------------------------------
c     Passed input/output arrays:
c
c         q      - heap (value)
c         iq     - heap (index)
c
c     ------------------------------------------------------------------
      integer q(n2),iq(n2)
c     ------------------------------------------------------------------
c     Local scalars:
c
c         vtmp   - tmp smallest element in heap (value)
c         ivtmp  - tmp smallest element in heap (index)
c         k      - heap counter
c         j      - heap counter (2*k)
c         szqd2  - sizeq/2
c
c     ------------------------------------------------------------------
      integer vtmp,ivtmp,k,j,szqd2

c     ------------------------------------------------------------------
c     Remove element from heap.
c     ------------------------------------------------------------------
      v     = q(1)
      iv    = iq(1)
      q(1)  = q(sizeq)
      iq(1) = iq(sizeq)
      sizeq = sizeq - 1
c     ------------------------------------------------------------------
c     Update heap to proper order.
c     ------------------------------------------------------------------
      k     = 1
      vtmp  = q(k)
      ivtmp = iq(k)
      szqd2 = sizeq/2
10    if (k .le. szqd2) then
         j = k + k
         if (j .lt. sizeq) then
            if (q(j) .gt. q(j + 1)) j = j + 1
         endif
         if (vtmp .gt. q(j)) then
            q(k)  = q(j)
            iq(k) = iq(j)
            k     = j
            goto 10
         endif
      endif
      q(k)  = vtmp
      iq(k) = ivtmp
c     ------------------------------------------------------------------
      return

c     end of removq
      end




      real function randp(ix)
c     ------------------------------------------------------------------
c     randp: Portable pseudo-random number generator.
c            Reference: L. Schrage, "A More Portable Fortran
c            Random Number Generator", ACM Transactions on
c            Mathematical Software, Vol. 2, No. 2, (June, 1979).
c     ------------------------------------------------------------------
      integer a,p,ix,b15,b16,xhi,xalo,leftlo,fhi,k
      data a/16807/,b15/32768/,b16/65536/,p/2147483647/
c     ------------------------------------------------------------------
      xhi    = ix/b16
      xalo   = (ix - xhi*b16)*a
      leftlo = xalo/b16
      fhi    = xhi*a + leftlo
      k      = fhi/b15
      ix     = (((xalo - leftlo*b16) - p) + (fhi - k*b15)*b16) + k
      if (ix .lt. 0) ix = ix + p

      randp = float(ix)*4.656612875e-10
c     ------------------------------------------------------------------
      return

c     end of randp
      end



      subroutine evalij(n,n2,i,j,xgain,f,d,a)
c     ------------------------------------------------------------------
c     evalij: Computes the gain in objective function by switching
c             the locations of facilities i and j (i < j).
c     ------------------------------------------------------------------
c     Passed input scalars:
c
c         n      - QAP dimension
c         n2     - n * n
c         i      - permutation array index
c         j      - permutation array index
c
c     ------------------------------------------------------------------
      integer n,n2,i,j
c     ------------------------------------------------------------------
c     Passed output scalar:
c
c         xgain  - gain achieved by swapping i and j in 
c                  permutation
c
c     ------------------------------------------------------------------
      integer xgain
c     ------------------------------------------------------------------
c     Passed input arrays:
c
c         f      - flow matrix
c         d      - distance matrix
c         a      - permutation vector
c   
c     ------------------------------------------------------------------
      integer f(n2),d(n2),a(n)
c     ------------------------------------------------------------------
c     Local scalars:
c
c         k      - do loop index
c         aim1tn - (a(i)-1)*n
c         ajm1tn - (a(j)-1)*n
c         akm1tn - (a(k)-1)*n
c         ai     - a(i)
c         aj     - a(j)
c         ak     - a(k)
c         im1tn  - (i-1)*n
c         jm1tn  - (j-1)*n
c         km1tn  - (k-1)*n
c         dtmp1  - reusable distance computation
c         dtmp2  - reusable distance computation
c
c     ------------------------------------------------------------------
      integer k,aim1tn,ajm1tn,akm1tn,ai,aj,ak,im1tn,jm1tn,km1tn,dtmp1,
     &        dtmp2

c     ------------------------------------------------------------------
      xgain  = 0
      ai     = a(i)
      aj     = a(j)
      aim1tn = (ai - 1)*n
      ajm1tn = (aj - 1)*n
      im1tn  = (i - 1)*n
      jm1tn  = (j - 1)*n
      km1tn  = 0
      do 10 k = 1,n
           if (k .ne. i .and. k .ne. j) then
                ak     = a(k)
                akm1tn = (ak - 1)*n
                dtmp1  = d(km1tn + i) - d(km1tn + j)
                dtmp2  = d(im1tn + k) - d(jm1tn + k)
                xgain  = xgain + 
     &                   dtmp1*(f(akm1tn + ai) - f(akm1tn + aj)) +
     &                   dtmp2*(f(aim1tn + ak) - f(ajm1tn + ak))
           endif
           km1tn = km1tn + n
10    continue
      dtmp1 = d(im1tn + j) - d(jm1tn + i)
      xgain = xgain + dtmp1*(f(aim1tn + aj) - f(ajm1tn + ai))

c     ------------------------------------------------------------------
      return

c     end of evalij
      end






















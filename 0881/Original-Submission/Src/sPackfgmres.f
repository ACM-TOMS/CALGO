C*
C*  Copyright (C) CERFACS 1998
C*
C*  SOFTWARE LICENSE AGREEMENT NOTICE - THIS SOFTWARE IS BEING PROVIDED TO
C*  YOU BY CERFACS UNDER THE FOLLOWING LICENSE. BY DOWN-LOADING, INSTALLING
C*  AND/OR USING THE SOFTWARE YOU AGREE THAT YOU HAVE READ, UNDERSTOOD AND
C*  WILL COMPLY WITH THESE FOLLOWING TERMS AND CONDITIONS.
C*
C*  1 - This software program provided in source code format ("the " Source
C*  Code ") and any associated documentation (the " Documentation ") are
C*  licensed, not sold, to you.
C*
C*  2 - CERFACS grants you a personal, non-exclusive, non-transferable and
C*  royalty-free right to use, copy or modify the Source Code and
C*  Documentation, provided that you agree to comply with the terms and
C*  restrictions of this agreement. You may modify the Source Code and
C*  Documentation to make source code derivative works, object code
C*  derivative works and/or documentation derivative Works (called "
C*  Derivative Works "). The Source Code, Documentation and Derivative
C*  Works (called " Licensed Software ") may be used by you for personal
C*  and non-commercial use only. " non-commercial use " means uses that are
C*  not or will not result in the sale, lease or rental of the Licensed
C*  Software and/or the use of the Licensed Software in any commercial
C*  product or service. CERFACS reserves all rights not expressly granted
C*  to you. No other licenses are granted or implied.
C*
C*  3 - The Source Code and Documentation are and will remain the sole
C*  property of CERFACS. The Source Code and Documentation are copyrighted
C*  works. You agree to treat any modification or derivative work of the
C*  Licensed Software as if it were part of the Licensed Software itself.
C*  In return for this license, you grant CERFACS a non-exclusive perpetual
C*  paid-up royalty-free license to make, sell, have made, copy, distribute
C*  and make derivative works of any modification or derivative work you
C*  make of the Licensed Software.
C*
C*  4- The licensee shall acknowledge the contribution of the Source Code
C*  (using the reference [1]) in any publication of material dependent upon
C*  upon the use of the Source Code. The licensee shall use reasonable
C*  endeavours to notify the authors of the package of this publication.
C*
C*  [1] V. Frayssé, L. Giraud, S. Gratton, and J. Langou, A set of GMRES 
C*    routines for real and complex arithmetics on high performance
C*    computers, CERFACS Technical Report TR/PA/03/3, public domain software
C*    available on www.cerfacs/algor/Softs, 2003
C*
C*  5- CERFACS has no obligation to support the Licensed Software it is
C*  providing under this license.
C*
C*  THE LICENSED SOFTWARE IS PROVIDED " AS IS " AND CERFACS MAKE NO
C*  REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED. BY WAY OF EXAMPLE,
C*  BUT NOT LIMITATION, CERFACS MAKE NO REPRESENTATIONS OR WARRANTIES OF
C*  MERCHANTIBILY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF
C*  THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE ANY THIRD
C*  PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. CERFACS WILL NOT
C*  BE LIABLE FOR ANY CONSEQUENTIAL, INCIDENTAL, OR SPECIAL DAMAGES, OR ANY
C*  OTHER RELIEF, OR FOR ANY CLAIM BY ANY THIRD PARTY, ARISING FROM YOUR
C*  USE OF THE LICENSED SOFTWARE.
C*
C*  6- For information regarding a commercial license for the Source Code
C*  and Documentation, please contact Mrs Campassens (campasse@cerfacs.fr)
C*
C*  7- This license is effective until terminated. You may terminate this
C*  license at any time by destroying the Licensed Software.
C*
C*    I agree all the terms and conditions of the above license agreement
C*
*
        subroutine drive_sfgmres(n,nloc,m,lwork,work,
     &                         irc,icntl,cntl,info,rinfo)
*
*  Purpose
*  =======
*    drive_sfgmres is the driver routine for solving the linear system 
*  Ax = b using the Flexible  Generalized Minimal Residual iterative method 
*  with preconditioning.
*  This solver is implemented with a reverse communication scheme: control
*  is returned to the user for computing the 
*   - matrix-vector product
*   - preconditioning
*   - dot product
*  See the User's Guide for an example of use.
*
*
* Written : June 1996
* Authors : Luc Giraud, Serge Gratton, V. Fraysse
*             Parallel Algorithms - CERFACS
*
* Updated : April 1997
* Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton
*             Parallel Algorithms - CERFACS
*
* Updated : May 1998
* Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton
*             Parallel Algorithms - CERFACS
* Comments: Totally rewritten from the sfgmres implementation
*
* Updated : March 2005 - L. Giraud
* Purpose : 1- Add the capability to avoid explicit residual calculation at restart
*           2- Include room in the workspace to store the results of the dot products
*           Fix the bugs that appeared when M > Nloc
*
*  Arguments
*  =========
*
*  n      (input) INTEGER.
*          On entry, the dimension of the problem.
*          Unchanged on exit.
*
* nloc    (input) INTEGER.
*          On entry, the dimension of the local problem.
*          In a parallel distributed envirionment, this corresponds
*          to the size of the subset of entries of the right hand side
*          and solution allocated to the calling process.
*          Unchanged on exit.
*
*
*  m       (input) INTEGER
*          Restart parameter, <= N. This parameter controls the amount
*          of memory required for matrix H (see WORK and H).
*          Unchanged on exit.
*
*  lwork   (input) INTEGER
*          size of the workspace
*          lwork >= m*m + m*(2*nloc+5) + 5*nloc+1
*
*  work    (workspace) real/real array, length lwork
*          work contains the required vector and matrices stored in the 
*          following order :
*            x  (nloc,1)       : computed solution.
*            b  (nloc,1)       : right hand side.
*            r0 (nloc,1)       : vector workspace.
*            w  (nloc,1)       : vector workspace.
*            H  (m+1,m+1)      : Hessenberg matrix (full storage).
*            yCurrent (m,1)    : solution of the current LS
*            xCurrent (nloc,1) : current iterate
*            rotSin (m,1)      : Sine of the Givens rotation
*            rotCos (m,1)      : Cosine of the Givens rotation
*            V  (nloc,m)       : Krylov basis.
*            Z  (nloc,m)       : Right preconditioned Krylov basis.
*
*  irc     (input/output) INTEGER array. length 7
*            irc(1) : REVCOM   used for reverse communication
*                             (type of external operation)
*            irc(2) : COLX     used for reverse communication
*            irc(3) : COLY     used for reverse communication
*            irc(4) : COLZ     used for reverse communication
*            irc(5) : NBSCAL   used for reverse communication
*            irc(6) : pointer on the first free location in the workspace
*            irc(7) : size (expressed in # items) of the free space 
*                     in the workspace
*
* icntl    (input) INTEGER array. length 7
*            icntl(1) : stdout for error messages
*            icntl(2) : stdout for warnings
*            icntl(3) : stdout for convergence history
*            icntl(4) : 0 - modified Gram-Schmidt
*                       1 - iterative modified Gram-Schmidt
*                       2 - classical Gram-Schmidt
*                       3 - iterative classical Gram-Schmidt
*            icntl(5) : 0 - default initial guess x_0 = 0 (to be set)
*                       1 - user supplied initial guess
*            icntl(6) : maximum number of iterations
*            icntl(7) : 0 - use recurence formula at restart
*                       1 - default compute the true residual at each restart
*
* cntl     (input) real array, length 3
*            cntl(1) : tolerance for convergence
*            cntl(2) : scaling factor for normwise perturbation on A
*            cntl(3) : scaling factor for normwise perturbation on b
*
* info     (output) INTEGER array, length 3
*            info(1) :  0 - normal exit
*                      -1 - n < 1
*                      -2 - m < 1
*                      -3 - lwork too small
*                      -4 - convergence not achieved after icntl(7) iterations
*                      -5 - precondition type not set by user
*            info(2) : if info(1)=0 - number of iterations to converge
*                      if info(1)=-3 - minimum workspace size necessary
*            info(3) : optimal size for the workspace
*
* rinfo    (output) real 
*            if info(1)=0 
*              rinfo : backward error for the preconditioned system
*
* Input variables
* ---------------
       integer n, nloc, lwork, icntl(*)
       real   cntl(*)
       real   sA, sb
* Output variables
* ----------------
       integer  info(*)
       real    rinfo
* Input/Output variables
* ----------------------
       integer  m, irc(*)
       real work(*)
* Local variables
* ---------------
       integer xptr, bptr, wptr, r0ptr, Vptr, Zptr, Hptr,dotptr
       integer yCurrent,rotSin, rotCos, xCurrent
       integer sizeWrk, newRestart
       integer iwarn, ierr, ihist, compRsd
       real    rn, rx,rc
       real DZRO
       parameter (DZRO = 0.0e0)
*
       integer icheck
       save icheck
       DATA icheck /0/
*
       intrinsic ifix, float
*
*       Executable statements :
*
       ierr  = icntl(1)
       iwarn = icntl(2)
       ihist = icntl(3)
       compRsd = icntl(7)
*
       if (ierr.lt.0) ierr = 6
*
       if (compRsd.eq.1) then
          sizeWrk  = m*m + m*(2*nloc+5) + 5*nloc+1
       else
          sizeWrk  = m*m + m*(2*nloc+5) + 6*nloc+1
       endif
*
       if (icheck.eq.0) then
* Check the value of the arguments
         if ((n.lt.1).or.(nloc.lt.1)) then
            write(ierr,*)
            write(ierr,*)' ERROR FGMRES : '
            write(ierr,*)'     N < 1 '
            write(ierr,*)
            info(1) = -1
            irc(1)  = 0
            return
         endif
         if (m.lt.1) then
            write(ierr,*)
            write(ierr,*)' ERROR FGMRES : '
            write(ierr,*)'     N < 1 '
            write(ierr,*)
            info(1) = -2
            irc(1)  = 0
            return
         endif
*
         if ((icntl(4).eq.2).or.(icntl(4).eq.3)) then
* the workspace should be large enough to store the m dot-products
            sizeWrk  = sizeWrk  + m
         else
            sizeWrk  = sizeWrk  + 1
         endif
*
         if (iwarn.ne.0) then
           write(iwarn,*)
           write(iwarn,*) ' WARNING FGMRES : '
           write(iwarn,*) '       For M = ',m,' optimal value '     
           write(iwarn,*) '       for LWORK =  ', sizeWrk
           write(iwarn,*)
         endif
*
         if ((icntl(4).lt.0).or.(icntl(4).gt.3)) then
           icntl(4) = 0
           if (iwarn.ne.0) then
             write(iwarn,*)
             write(iwarn,*) ' WARNING FGMRES : '
             write(iwarn,*) '       Undefined orthogonalisation '  
             write(iwarn,*) '       Default MGS '
             write(iwarn,*)
           endif
         endif
         if ((icntl(5).ne.0).and.(icntl(5).ne.1)) then
           icntl(5) = 0
           if (iwarn.ne.0) then
             write(iwarn,*)
             write(iwarn,*) ' WARNING FGMRES : '
             write(iwarn,*) '       Undefined intial guess'
             write(iwarn,*) '       Default x0 = 0 '
             write(iwarn,*)
           endif
         endif
         if (icntl(6).le.0) then
           icntl(6) = n
           if (iwarn.ne.0) then
             write(iwarn,*)
             write(iwarn,*) ' WARNING FGMRES : '
             write(iwarn,*) '       Negative max number of iterations'  
             write(iwarn,*) '       Default N '
             write(iwarn,*)
           endif
         endif
         if ((icntl(7).ne.0).and.(icntl(7).ne.1)) then
           icntl(7) = 1
           write(iwarn,*)
           write(iwarn,*) ' WARNING FGMRES :'
           write(iwarn,*) '       Undefined strategy for the residual'
           write(iwarn,*) '       at restart'
           write(iwarn,*) '       Default 1 '
           write(iwarn,*)
         endif
* Check if the restart parameter is correct and if the size of the
*  workspace is big enough for the restart.
* If not try to fix correctly the parameters
*
         if ((m .gt. n).or.(lwork.lt.sizeWrk)) then
           if (m .gt. n) then
             m = n
             if (iwarn.ne.0) then
               write(iwarn,*)
               write(iwarn,*) ' WARNING FGMRES : '
               write(iwarn,*) '       Parameter M bigger than N'  
               write(iwarn,*) '       New value for M ',m
               write(iwarn,*)
             endif
             if (compRsd.eq.1) then
                sizeWrk  = m*m + m*(2*nloc+5) + 5*nloc+1
             else
                sizeWrk  = m*m + m*(2*nloc+5) + 6*nloc+1
             endif
             if ((icntl(4).eq.2).or.(icntl(4).eq.3)) then
* the workspace should be large enough to store the m dot-products
                sizeWrk  = sizeWrk  + m
             else
                sizeWrk  = sizeWrk  + 1
             endif
           endif 
           if ((lwork.lt.sizeWrk).and.(n.eq.nloc)) then
* Compute the maximum size of the restart according to the memory space
             rn  = float(n)
             rx  = 2.0*rn + 5.0
             rc  = 5.0*rn + 1 - float(lwork)
*
* Update the linear part of the second order equation to be solved
             if ((icntl(4).eq.2).or.(icntl(4).eq.3)) then
               rx = rx + 1
             endif
* Update the constant part of the second order equation to be solved
*             
             if (icntl(7).eq.0) then
               rc = rc + 2.0*rn
             endif
             newRestart = ifix((-rx+sqrt(rx**2-4.0*rc))/2.0)
             if (newRestart.gt.0) then
               m = newRestart
               if (iwarn.ne.0) then
                 write(iwarn,*)
                 write(iwarn,*)' WARNING FGMRES : '
                 write(iwarn,*) '      Workspace too small for M'  
                 write(iwarn,*)'       New value for M ',m
                 write(iwarn,*)
               endif
             else
               write(ierr,*)
               write(ierr,*)' ERROR FGMRES : '
               write(ierr,*)'     Not enough space for the problem'
               write(ierr,*)'     the space does not permit any m'
               write(ierr,*)
               info(1) = -3
               irc(1)  = 0
               return
             endif
           endif
           if ((lwork.lt.sizeWrk).and.(n.ne.nloc)) then
              write(ierr,*)
              write(ierr,*)' ERROR FGMRES : '
              write(ierr,*)'     Not enough space for the problem'
              write(ierr,*)
              info(1) = -3
              irc(1)  = 0
              return
           endif
         endif
*
         info(3) = sizeWrk
         icheck = 1
*
* save the parameters in the history file
*
         if (ihist.ne.0) then
           write(ihist,'(10x,A39)') 'CONVERGENCE HISTORY FOR FGMRES'
           write(ihist,*)
           write(ihist,'(A30,I2)') 'Errors are displayed in unit: ',ierr 
           if (iwarn.eq.0) then
             write(ihist,'(A27)') 'Warnings are not displayed:'
           else
             write(ihist,'(A32,I2)') 'Warnings are displayed in unit: ',
     &                               iwarn
           endif 
           write(ihist,'(A13,I7)') 'Matrix size: ',n
           write(ihist,'(A19,I7)') 'Local matrix size: ',nloc
           write(ihist,'(A9,I7)') 'Restart: ',m
           if (icntl(4).eq.0) then
             write(ihist,'(A21)') 'Modified Gram-Schmidt'
           elseif (icntl(4).eq.1) then
             write(ihist,'(A31)') 'Iterative modified Gram-Schmidt'
           elseif (icntl(4).eq.2) then
             write(ihist,'(A22)') 'Classical Gram-Schmidt'
           else
             write(ihist,'(A32)') 'Iterative classical Gram-Schmidt'
           endif
           if (icntl(5).eq.0) then
             write(ihist,'(A29)') 'Default initial guess x_0 = 0'
           else
             write(ihist,'(A27)') 'User supplied initial guess'
           endif
           write(ihist,'(A30,I5)') 'Maximum number of iterations: ',
     &                              icntl(6)
           if (icntl(7).eq.1) then
             write(ihist,'(A33)') 'True residual computed at restart'
           else
             write(ihist,'(A30)') 'Recurrence residual at restart'
           endif
           write(ihist,'(A27,E8.2)') 'Tolerance for convergence: ', 
     &                                cntl(1) 
* 
           write(ihist,'(A53)') 
     &       'Backward error on the unpreconditioned system Ax = b:'
           sA       = cntl(2)
           sb       = cntl(3)
           if ((sA.eq.DZRO).and.(sb.eq.DZRO)) then
             write(ihist,'(A39)') 
     &       '    the residual is normalised by ||b||'
           else
             write(ihist,'(A34)') '    the residual is normalised by '
             write(ihist,'(A8,E8.2,$)') '        ', sA
             write(ihist,'(A11,E8.2)') ' * ||x|| + ', sb
           endif
*
           write(ihist,'(A31,I7)') 'Optimal size for the workspace:',
     &                              info(3)
           write(ihist,*) 
           write(ihist,'(A32,$)') 'Convergence history: b.e. on the'
           write(ihist,'(A22)') ' preconditioned system'
           write(ihist,'(A11,$)') ' Iteration '
           write(ihist,'(A27)') '  Arnoldi b.e.    True b.e.'
         endif
*
       endif
* setup some pointers on the workspace
       xptr     = 1
       bptr     = xptr + nloc
       dotptr   = bptr + nloc
       if ((icntl(4).eq.2).or.(icntl(4).eq.3)) then
         r0ptr = dotptr + m
       else
         r0ptr = dotptr + 1
       endif
       wptr     = r0ptr + nloc
       Hptr     = wptr + nloc
       yCurrent = Hptr + (m+1)*(m+1)
       xCurrent = yCurrent + m
       rotSin   = xCurrent + nloc
       rotCos   = rotSin + m
       Vptr     = rotCos + m
       Zptr     = lwork - (m*nloc)+1
*
       call sfgmres(nloc,m,lwork,work(bptr),work(xptr),work(dotptr),
     &            work(Hptr),work(wptr),work(r0ptr),work(Vptr),
     &            work(Zptr),work(yCurrent),work(xCurrent),
     &            work(rotSin),work(rotCos),irc,icntl,cntl,info,rinfo)
*
       if (irc(1).eq.0) then
         icheck = 0
       endif
*
       return
       end
*
        subroutine sfgmres(n,m,lwork,b,x,dot,H,w,r0,V,Z,yCurrent,
     &                xCurrent,rotSin,rotCos,irc,icntl,cntl,info,rinfo)
*
*
*  Purpose
*  =======
*  sfgmres solves the linear system Ax = b using the
*  Flexible Generalized Minimal Residual iterative method
*
* When preconditioning is used we solve :
*     A M_2^{-1} y = b
*     x = M_2^{-1} y
*
*   Convergence test based on the normwise backward error for
*  the preconditioned system
*
* Written : June 1996
* Authors : Luc Giraud, Serge Gratton, V. Fraysse
*             Parallel Algorithms - CERFACS
*
* Updated : April 1997
* Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton
*             Parallel Algorithms - CERFACS
*
* Updated : March 1998
*           Pb with F90 on DEC ws
*           cure : remove "ZDSCAL" when used to initialize vectors to zero
*
* Updated : May 98
* Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton
*             Parallel Algorithms - CERFACS
*
* Updated : February 2001 - Luc Giraud
* Purpose : In complex version, initialization to zero performed in complex
*           arithmetic to avoid implicit conversion by the compiler.
*
* Updated : July 2001 - L. Giraud, J. Langou
* Purpose : Avoid to compute the approximate solution at each step of
*           the Krylov space construction when sA is zero.
*
* Updated : November 2002 - S. Gratton
* Purpose : Use Givens rotations conform to the classical definition.
*           No impact one the convergence history.
*
* Updated : November 2002 - L. Giraud
* Purpose : Properly handle the situation when the convergence is obtained
*           exactly at the "IterMax" iteration
*
* Updated : March 2005 - L. Giraud
* Purpose :  Use Givens rotations from BLAS
*
* Updated : March 2005 - L. Giraud
* Purpose : Include room in the workspace to store the results of the dot products
*           Fix the bugs that appeared when M > Nloc
*
*  Arguments
*  =========
*
*  n       (input) INTEGER.
*           On entry, the dimension of the problem.
*           Unchanged on exit.
*
*
*  m        (input) INTEGER
*           Restart parameter, <= N. This parameter controls the amount
*           of memory required for matrix H (see WORK and H).
*           Unchanged on exit.
*
*  lwork    (input) INTEGER
*           the size of the workspace provided by the user.
*           Unchanged on exit.
*
*  b        (input) real/real
*           Right hand side of the linear system.
*
*  x        (output) real/real
*           Computed solution of the linear system.
*
*  dot      (workspace) real/real
*           Store the results of the dot product calculation
*
*  H        (workspace)  real/real
*           Hessenberg matrix built within dgmres
*
*  w        (workspace)  real/real
*           Vector used as temporary storage
*
*  r0       (workspace)  real/real
*           Vector used as temporary storage
*
*  V        (workspace)  real/real
*            Basis computed by the Arnoldi's procedure.
*
* Z         (workspace)  real/real
*            Set of preconditioned V vectors
*  
*  yCurrent (workspace) real/real
*           solution of the current LS
*
* xCurrent  (workspace) real/real
*           current iterate
*
* rotSin    (workspace) real/real
*           Sine of the Givens rotation
*
* rotCos    (workspace) real/real
*           Cosine of the Givens rotation
*
*  irc     (input/output) INTEGER array. length 7
*            irc(1) : REVCOM   used for reverse communication
*                             (type of external operation)
*            irc(2) : COLX     used for reverse communication
*            irc(3) : COLY     used for reverse communication
*            irc(4) : COLZ     used for reverse communication
*            irc(5) : NBSCAL   used for reverse communication
*            irc(6) : pointer on the first free location in the workspace
*            irc(7) : size (expressed in # items) of the free space
*                     in the workspace
*
* icntl    (input) INTEGER array. length 6
*            icntl(1) : stdout for error messages
*            icntl(2) : stdout for warnings
*            icntl(3) : stdout for convergence history
*            icntl(4) : 0 - modified Gram-Schmidt
*                       1 - iterative modified Gram-Schmidt
*                       2 - classical Gram-Schmidt
*                       3 - iterative classical Gram-Schmidt
*            icntl(5) : 0 - default initial guess x_0 = 0 (to be set)
*                       1 - user supplied initial guess
*            icntl(6) : maximum number of iterations
*            icntl(7) : 1 - default compute the true residual at each restart
*                       0 - use recurence formula at restart
*
* cntl     (input) real array, length 3
*            cntl(1) : tolerance for convergence
*            cntl(2) : scaling factor for normwise perturbation on A
*            cntl(3) : scaling factor for normwise perturbation on b
*
* info     (output) INTEGER array, length 2
*            info(1) :  0 - normal exit
*                      -1 - n < 1
*                      -2 - m < 1
*                      -3 - lwork too small
*                      -4 - convergence not achieved after icntl(7) iterations
*            info(2) : if info(1)=0 - number of iterations to converge
*                      if info(1)=-3 - minimum workspace size necessary
*            info(3) : optimal size for the workspace
*
* rinfo    (output) real 
*            if info(1)=0 
*              rinfo : backward error for the linear system
*
* Input variables
* ---------------
        integer  n, m, lwork, icntl(*)
        real b(n)
        real    cntl(3)
*
* Output variables
* ----------------
       integer  info(3)
       real    rinfo
*
* Input/Output variables
* ----------------------
       integer  irc(7)
       real x(n), dot(m), H(m+1,m+1), w(n)
       real r0(n), V(n,m+1), Z(n,m)
       real xCurrent(n), rotSin(m), yCurrent(m)   
       real rotCos(m)
*
* Local variables
* ---------------
       integer  j, jH, iterOut, nOrtho, iterMax, initGuess, iOrthog
       integer  xptr, bptr, wptr, r0ptr, Vptr, Hptr, yptr, xcuptr
       integer  Zptr, dotptr
       integer  iwarn, ihist, compRsd
       integer  MGS, IMGS, CGS, ICGS
       parameter (MGS = 0, IMGS = 1, CGS = 2, ICGS = 3)
       real    beta, bn, sA, sb, bea, be, temp
       real    dloo, dnormw, dnormx, dnormres, trueNormRes
       real dVi, aux
       real auxHjj, auxHjp1j
*
       real ZERO, ONE
       parameter (ZERO = 0.0e0, ONE = 1.0e0)
       real DZRO,DONE
       parameter (DZRO = 0.0e0, DONE = 1.0e0)
*
* External functions
* ------------------
       real    snrm2
       external snrm2
*
* Saved variables
* ---------------
      save iterOut, jH, beta, bn, dnormres, retlbl, j
      save sA, sb, dnormx, trueNormRes, bea, be
      save dloo, nOrtho, compRsd
*
* Intrinsic function
* ------------------
        intrinsic abs, sqrt 
*
* Reverse communication variables
* -------------------------------
       integer matvec, precondRight, prosca
       parameter(matvec=1, precondRight=3, prosca=4)
       integer retlbl
       DATA retlbl /0/
*
*       Executable statements
*
* setup some pointers on the workspace
       xptr     = 1
       bptr     = xptr + n
       dotptr   = bptr + n
       if ((icntl(4).eq.2).or.(icntl(4).eq.3)) then
         r0ptr = dotptr + m
       else
         r0ptr = dotptr + 1
       endif
       wptr     = r0ptr + n
       Hptr     = wptr + n
       yptr     = Hptr + (m+1)*(m+1)
       xcuptr   = yptr + m
       Vptr     = xcuptr + n +  2*m
       Zptr     = lwork - (m*n)+1
*
       iwarn      = icntl(2)
       ihist      = icntl(3)
       iOrthog    = icntl(4)
       initGuess  = icntl(5)
       iterMax    = icntl(6)
*
       if (retlbl.eq.0) then
         compRsd    = icntl(7)
       endif
*
       if (retlbl.ne.0) then
          if (retlbl.eq.5) then
            goto 5
          else if (retlbl.eq.11) then
            goto 11
          else if (retlbl.eq.18) then
            goto 18
          else if (retlbl.eq.21) then
            goto 21
          else if (retlbl.eq.26) then
            goto 26
          else if (retlbl.eq.32) then
            goto 32
          else if (retlbl.eq.33) then
            goto 33
          else if (retlbl.eq.34) then
            goto 34 
          else if (retlbl.eq.38) then
            goto 38
          else if (retlbl.eq.41) then
            goto 41
          else if (retlbl.eq.43) then
            goto 43
          else if (retlbl.eq.61) then
            goto 61
          else if (retlbl.eq.68) then
            goto 68
          endif
        endif
*
*
* intialization of various variables
*
        iterOut  = 0
        beta     = DZRO
*
        if (initGuess.eq.0) then
          do j=1,n
            x(j) = ZERO
          enddo
        endif
*
*        bn = snrm2(n,b,1)
*
        irc(1) = prosca
        irc(2) = bptr
        irc(3) = bptr
        irc(4) = dotptr
        irc(5) = 1
        retlbl = 5
        return
 5      continue
        bn = sqrt((dot(1)))
*
        if (bn.eq.DZRO) then
          do j=1,n
            x(j) = ZERO
          enddo  
          if (iwarn.ne.0) then
            write(iwarn,*)
            write(iwarn,*) ' WARNING FGMRES : '
            write(iwarn,*) '       Null right hand side'
            write(iwarn,*) '       Solution set to zero'
            write(iwarn,*)
          endif
          jH = 0
          bea = DZRO
          be  = DZRO
          write(ihist,'(I5,11x,E8.2,$)') jH,bea
          write(ihist,'(7x,E8.2)') be
          info(1)  = 0
          info(2)  = 0
          rinfo    = DZRO
          irc(1)   = 0
          retlbl = 0
          return
        endif
*
* Compute the scaling factor for the backward error on the 
*  unpreconditioned sytem
*
       sA       = cntl(2)
       sb       = cntl(3)
       if ((sA.eq.DZRO).and.(sb.eq.DZRO)) then
         sb = bn
       endif
*
*
* Compute the first residual
*           Y = AX : r0 <-- A x
*
* The residual is computed only if the initial guess is not zero
*
       if (initGuess.ne.0) then
         irc(1) = matvec
         irc(2) = xptr
         irc(4) = r0ptr
         retlbl = 11
         return
       endif
 11    continue
       if (initGuess.ne.0) then
         do j=1,n
           r0(j) = b(j)-r0(j)
         enddo
       else
         call scopy(n,b,1,r0,1)
       endif 
*
* Normalize the first Krylov vector
*
       call scopy(n,r0,1,w,1)
*
*       beta = snrm2(n,w,1)
*
*
        irc(1) = prosca
        irc(2) = wptr
        irc(3) = wptr
        irc(4) = dotptr
        irc(5) = 1
        retlbl = 18
        return
 18     continue
        beta = sqrt((dot(1)))
*
        if (beta .eq. DZRO) then
*  The residual is exactly zero : x is the exact solution
          info(1) = 0
          info(2) = 0
          rinfo    = DZRO
          irc(1)   = 0
          retlbl = 0 
          jH = 0
          bea = DZRO
          be  = DZRO
          write(ihist,'(I5,11x,E8.2,$)') jH,bea
          write(ihist,'(7x,E8.2)') be
          if (iwarn.ne.0) then
            write(iwarn,*)
            write(iwarn,*) ' WARNING GMRES : '
            write(iwarn,*) '       Intial residual is zero'
            write(iwarn,*) '       initial guess is solution'
            write(iwarn,*)
          endif
          return
        endif
*
* Update the free space descriptor
*
        irc(6) = Vptr
        irc(7) = lwork - irc(6)
*
        aux = ONE/beta
        do j=1,n
          V(j,1) = ZERO
        enddo
        call saxpy(n,aux,w,1,V(1,1),1)
*
*
        irc(6) = irc(6) + n 
        irc(7) = irc(7) - n 
*
*       Most outer loop : sfgmres iteration
*
*       REPEAT
 7      continue
*
*
        H(1,m+1)=beta
        do j=1,m
          H(j+1,m+1) = ZERO
        enddo
*
*        Construction of the hessenberg matrix WORK and of the orthogonal
*        basis V such that AV=VH 
*
        jH = 1
 10     continue
* Remark : this  do loop has been written with a while do
*          because the
*               " do jH=1,restart "
*         fails with the reverse communication.
*      do  jH=1,restart
*
*
* Compute the preconditioned residual
*
*           Z(1,jH) <-- M_2^{-1} V(1,jH)
*
           irc(7) = irc(7) - n 
*
           irc(1) = precondRight
           irc(2) = vptr + (jH-1)*n
           irc(4) = Zptr + (m-jH)*n
           retlbl = 21
           return
 21      continue
*
*           Y = AX : w <-- A Z(1,jH)
*
         irc(1) = matvec
         irc(2) = Zptr + (m-jH)*n
         irc(4) = wptr
         retlbl = 26
         return
 26      continue
*
*
* Orthogonalization using either MGS or IMGS
*  
* initialize the Hessenberg matrix to zero in order to be able to use
*     IMGS as orthogonalization procedure.
       do j=1,jH
         H(j,jH) = ZERO
       enddo
       nOrtho = 0
 19    continue
         nOrtho = nOrtho +1
         dloo   = DZRO
*
         if ((iOrthog.eq.MGS).or.(iOrthog.eq.IMGS)) then
* MGS
*
*           do j=1,jH
*
            j = 1
*           REPEAT
         endif
 23      continue
         if ((iOrthog.eq.MGS).or.(iOrthog.eq.IMGS)) then
*
*             dVi     = sdot(n,V(1,j),1,w,1)
*
              irc(1) = prosca
              irc(2) = vptr + (j-1)*n
              irc(3) = wptr
              irc(4) = dotptr
              irc(5) = 1
              retlbl = 32
              return
          endif
 32       continue
          if ((iOrthog.eq.MGS).or.(iOrthog.eq.IMGS)) then
              dVi     = dot(1)
              H(j,jH) = H(j,jH) + dVi
              dloo    = dloo + abs(dVi)**2
              aux = -ONE*dVi
              call saxpy(n,aux,V(1,j),1,w,1)
              j = j + 1
              if (j.le.jH) goto 23
*          enddo_j
         else
* CGS
* produit scalaire groupe
*
*           call sgemv('C',n,jH,ONE,V(1,1),n,w,1,ZERO,r0,1)
*
            irc(1) = prosca
            irc(2) = vptr
            irc(3) = wptr
            irc(4) = dotptr
            irc(5) = jH
            retlbl = 34
            return
          endif
 34       continue
          if ((iOrthog.eq.CGS).or.(iOrthog.eq.ICGS)) then
*
           call saxpy(jH,ONE,dot,1,H(1,jH),1)
           call sgemv('N',n,jH,-ONE,V(1,1),n,dot,1,ONE,w,1)
           dloo = snrm2(jH,dot,1)**2
         endif
*
*         dnormw = snrm2(n,w,1)
*
         irc(1) = prosca
         irc(2) = wptr
         irc(3) = wptr
         irc(4) = dotptr
         irc(5) = 1
         retlbl = 33
         return
 33      continue
         dnormw = sqrt((dot(1)))
*
         if ((iOrthog.eq.CGS).or.(iOrthog.eq.ICGS)) then
* IMGS / CGS orthogonalisation
           dloo = sqrt(dloo)
* check the orthogonalization quality
           if ((dnormw.le.dloo).and.(nOrtho.lt.3)) then
             goto 19
           endif
         endif
*
       H(jH+1,jH) = dnormw
       if ((jH.lt.m).or.(icntl(7).eq.0)) then
         aux = ONE/dnormw
         do j=1,n
           V(j,jH+1) = ZERO
         enddo
         call saxpy(n,aux,w,1,V(1,jH+1),1)
       endif
*
       irc(6) = irc(6) + n
       irc(7) = irc(7) - n
* Apply previous Givens rotations to the new column of H
       do j=1,jH-1
         call srot(1, H(j,jH), 1, H(j+1,jH), 1, (rotCos(j)),
     &             rotSin(j))
       enddo
       auxHjj = H(jH,jH)
       auxHjp1j= H(jH+1,jH)
       call srotg(auxHjj, auxHjp1j, temp,rotSin(jH))
       rotCos(jH) = temp
* Apply current rotation to the rhs of the least squares problem
       call srot(1, H(jH,m+1), 1, H(jH+1,m+1), 1, (rotCos(jH)),
     &            rotSin(jH))
*
* zabs(H(jH+1,m+1)) is the residual computed using the least squares
*          solver
* Complete the QR factorisation of the Hessenberg matrix by apply the current
* rotation to the last entry of the collumn
       call srot(1, H(jH,jH), 1, H(jH+1,jH), 1, (rotCos(jH)), 
     &           rotSin(jH))
       H(jH+1,jH) = ZERO
*
* Get the Least square residual
*
       dnormres = abs(H(jH+1,m+1))
       if (sA.ne.DZRO) then
*
* Compute the solution of the current linear least squares problem
*
          call scopy(jH,H(1,m+1),1,yCurrent,1)
          call strsv('U','N','N',jH,H,m+1,yCurrent,1)
*
* Permute yCurrent according to the reverse ordering storage of Z before to
* compute the current solution
*
          do j=1,jH/2
            aux              = yCurrent(j)
            yCurrent(j)      = yCurrent(jH-j+1)
            yCurrent(jH-j+1) = aux
          enddo
*
* Compute the value of the new iterate 
*
          call sgemv('N',n,jH,ONE,Z(1,m-jH+1),n,
     &            yCurrent,1,ZERO,xCurrent,1)
*
          call scopy(n,xCurrent,1,r0,1)
          call scopy(n,x,1,xCurrent,1)
          call saxpy(n,ONE,r0,1,xCurrent,1)
*
*         dnormx = snrm2(n,xCurrent,1)
*
          irc(1) = prosca
          irc(2) = xcuptr
          irc(3) = xcuptr
          irc(4) = dotptr
          irc(5) = 1
          retlbl = 38
          return
       else
          dnormx    = DONE
       endif
 38    continue
       if (sA.ne.DZRO) then
         dnormx = sqrt((dot(1)))
       endif
*
       bea = dnormres/(sA*dnormx+sb)
*
* Check the convergence based on the Arnoldi Backward error for the
* preconditioned system
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then  
* 
* The Arnoldi Backward error indicates that sgmres might have converge
* enforce the calculation of the true residual at next restart
         compRsd = 1
*
*  If the update of X has not yet been performed
         if (sA.eq.DZRO) then
*
* Compute the solution of the current linear least squares problem
*
            call scopy(jH,H(1,m+1),1,yCurrent,1)
            call strsv('U','N','N',jH,H,m+1,yCurrent,1)
*
* Permute yCurrent according to the reverse ordering storage of Z before to
* compute the current solution
*
            do j=1,jH/2
              aux              = yCurrent(j)
              yCurrent(j)      = yCurrent(jH-j+1)
              yCurrent(jH-j+1) = aux
            enddo
*
* Compute the value of the new iterate 
*
            call sgemv('N',n,jH,ONE,Z(1,m-jH+1),n,
     &            yCurrent,1,ZERO,xCurrent,1)
*
            call scopy(n,xCurrent,1,r0,1)
* Update the current solution
            call scopy(n,x,1,xCurrent,1)
            call saxpy(n,ONE,r0,1,xCurrent,1)
         endif
*
         call scopy(n,xCurrent,1,r0,1)
* Compute the true residual, the Arnoldi one may be unaccurate
*
*           Y = AX : w  <-- A r0
*
         irc(1) = matvec
         irc(2) = r0ptr
         irc(4) = wptr
         retlbl = 41
         return
       endif
 41    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
*
         do j=1,n
           w(j) = b(j) - w(j)
         enddo
* Compute the norm of the unpreconditioned residual
*
*        trueNormRes = snrm2(n,w,1)
*
         irc(1) = prosca
         irc(2) = wptr
         irc(3) = wptr
         irc(4) = dotptr
         irc(5) = 1
         retlbl = 43
         return
       endif
 43    continue
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then
         trueNormRes = sqrt((dot(1)))
         dnormres    = trueNormRes
*
         be = dnormres/(sA*dnormx+sb)
* Save the backward error on a file if convergence history requested
         if (ihist.ne.0) then
           write(ihist,'(I5,11x,E8.2,$)') iterOut*m+jH,bea
           write(ihist,'(7x,E8.2)') be
         endif
*
       endif
*
*
* Check again the convergence
       if ((bea.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then   
         if ((be.le.cntl(1)).or.(iterOut*m+jH.ge.iterMax)) then   
            call scopy(n,xCurrent,1,x,1)
*
* Return the backward errors
             rinfo = be
             if (be.le.cntl(1)) then
               info(1) = 0
               if (ihist.ne.0) then
                 write(ihist,*)
                 write(ihist,'(A20)') 'Convergence achieved'
               endif
             else
               if (iwarn.ne.0) then
                 write(iwarn,*)
                 write(iwarn,*) ' WARNING FGMRES : '
                 write(iwarn,*) '       No convergence after '
                 write(iwarn,*) iterOut*m+jH,' outer iterations '
                 write(iwarn,*)
               endif
               if (ihist.ne.0) then
                 write(ihist,*)
                 write(ihist,*) ' WARNING FGMRES : '
                 write(ihist,*) '       No convergence after '
                 write(ihist,*) iterOut*m+jH,' outer iterations '
                 write(ihist,*)
               endif
               info(1) = -4
             endif
             if (ihist.ne.0) then
               write(ihist,'(A27,$)') 'B.E. on the '
               write(ihist,'(A10,E8.2)') 'system:   ', rinfo
             endif
             info(2) = iterOut*m+jH
             if (ihist.ne.0) then
               write(ihist,'(A10,I2)') 'info(1) = ',info(1)
               write(ihist,'(A32,I5)') 
     &                'Number of iterations (info(2)): ',info(2)  
             endif
             irc(1)  = 0
             retlbl  = 0
             return
           endif
       else
* Save the backward error on a file if convergence history requested
         if (ihist.ne.0) then
           write(ihist,'(I5,11x,E8.2,$)') iterOut*m+jH,bea
           write(ihist,'(9x,A2)') '--'
         endif
*
       endif  
*
       jH = jH + 1
       if (jH.le.m) then
         goto 10
       endif
*
       iterOut = iterOut + 1
*
* we have completed the Krylov space construction, we restart if
* we have not yet exceeded the maximum number of iterations allowed.
*
       if ((sA.eq.DZRO).and.(bea.gt.cntl(1))) then
*
* Compute the solution of the current linear least squares problem
*
         jH = jH - 1
         call scopy(jH,H(1,m+1),1,yCurrent,1)
         call strsv('U','N','N',jH,H,m+1,yCurrent,1)
*
* Permute yCurrent according to the reverse ordering storage of Z before to
* compute the current solution
*
         do j=1,jH/2
           aux              = yCurrent(j)
           yCurrent(j)      = yCurrent(jH-j+1)
           yCurrent(jH-j+1) = aux
         enddo
*
* Compute the value of the new iterate 
*
         call sgemv('N',n,jH,ONE,Z(1,m-jH+1),n,
     &            yCurrent,1,ZERO,xCurrent,1)
*
         call scopy(n,xCurrent,1,r0,1)
         call scopy(n,x,1,xCurrent,1)
         call saxpy(n,ONE,r0,1,xCurrent,1)
       endif
*
         call scopy(n,xCurrent,1,x,1)
*
* Compute the residual
*
         call scopy(n,x,1,w,1)
*
         if (compRsd.eq.1) then
           irc(1) = matvec
           irc(2) = wptr
           irc(4) = r0ptr
           retlbl = 61
           return
         endif
 61      continue
         if (compRsd.eq.1) then
           do j=1,n
             r0(j) = b(j) - r0(j)
           enddo
*
           call scopy(n,r0,1,w,1)
*
*          beta = snrm2(n,w,1)
*
           irc(1) = prosca
           irc(2) = wptr
           irc(3) = wptr
           irc(4) = dotptr
           irc(5) = 1
           retlbl = 68
           return
         endif
 68      continue
*
         if (compRsd.eq.1) then
           beta = sqrt((dot(1)))
         else
* Use recurrence to approximate the residual at restart
           beta = abs(H(m+1,m+1))
* Apply the Givens rotation is the reverse order
           do j=m,1,-1
             H(j,m+1)   = ZERO
             call srot(1, H(j,m+1), 1, H(j+1,m+1), 1,
     &               (rotCos(j)), -rotSin(j))
           enddo
*
* On applique les vecteurs V
*
           call sgemv('N',n,m+1,ONE,v,n,H(1,m+1),1,ZERO,w,1)
*
         endif
*
         do j=1,n
            V(j,1) = ZERO
         enddo
         aux = ONE/beta
         call saxpy(n,aux,w,1,V(1,1),1)
*
         irc(6) = Vptr + n
         irc(7) = lwork - irc(6)
*
         goto 7
*
        end
*
*
        subroutine init_sfgmres(icntl,cntl)
*
*  Purpose
*  =======
*    Set default values for the parameters defining the characteristics
* of the Gmres algorithm.
*  See the User's Guide for an example of use.
*
*
* Written : April 1997
* Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton
*             Parallel Algorithms - CERFACS
*
*  Arguments
*  =========
*
* icntl    (input) INTEGER array. length 6
*            icntl(1) : stdout for error messages
*            icntl(2) : stdout for warnings
*            icntl(3) : stdout for convergence history
*            icntl(4) : 0 - modified Gram-Schmidt
*                       1 - iterative modified Gram-Schmidt
*                       2 - classical Gram-Schmidt
*                       3 - iterative classical Gram-Schmidt
*            icntl(5) : 0 - default initial guess x_0 = 0 (to be set)
*                       1 - user supplied initial guess
*            icntl(6) : maximum number of iterations
*            icntl(7) : 1 - default compute the true residual at each restart
*                       0 - use recurence formaula at restart
*
* cntl     (input) real array, length 5
*            cntl(1) : tolerance for convergence
*            cntl(2) : scaling factor for normwise perturbation on A
*            cntl(3) : scaling factor for normwise perturbation on b
*
* Output variables
* ----------------
       integer icntl(*)
       real   cntl(*)
*
       icntl(1) = 6
       icntl(2) = 6
       icntl(3) = 0
       icntl(4) = 0
       icntl(5) = 0
       icntl(6) = 100
       icntl(7) = 1
* 
       cntl(1) = 1.0 e -7
       cntl(2) = 0.0 e 0
       cntl(3) = 0.0 e 0
*
       return
       end

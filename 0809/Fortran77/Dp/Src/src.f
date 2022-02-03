c
      subroutine preqn ( n, m, iop, iprob, jcg, s, y, r, z, w, lw, iw, 
     *                   liw, build, info, mssg )
       
      integer  n, m, iop, iprob, jcg, lw, iw(*), liw, info, mssg
      double   precision s(n), y(n), r(n), z(n), w(*)
      logical  build
c
c***********************************************************************
c
c     Algorithm PREQN: Fortran Subroutines for Preconditioning the 
c     Conjugate Gradient Method
c
c
c        Jose Luis Morales                 Jorge Nocedal                      
c        Departamento de Matematicas       Department of Electrical
c        Instituto Tecnologico Autonomo    and Computer Engineering,
c        de Mexico. Rio Hondo 1,           Northwestern University
c        Mexico DF, CP 04530               Evanston, Il 60208,
c        MEXICO.                           USA.
c        jmorales@gauss.rhon.itam.mx       nocedal@ece.nwu.edu
c                                          www.ece.nwu.edu/~nocedal
c
c***********************************************************************
c
c     1. INTRODUCTION
c     ---------------
c
c     This collection of subroutines computes preconditioners for the
c     conjugate gradient method (CG), when applied to a sequence of 
c     linear systems
c
c                  A_i x = b_i,  i=1, 2,..., t,                    (1)
c
c
c     where the matrices A_i are symmetric, positive definite and of 
c     size n, and  x and  b_i are n-dimensional vectors. It is assumed 
c     that the matrices A_i are related (or constant), but the right 
c     hand side vectors b_i are arbitrary. 
c
c     The routines compute quasi-Newton preconditioners H_i. Thus
c     we solve (1) for all i>1 by means of the CG method applied to the 
c     preconditioned systems
c
c
c                  H_i A_i x = H_i b_i     i=2,3,...,t            (2)
c
c
c     The preconditioners are constructed using limited memory BFGS 
c     updating. Briefly, the preconditioner H_{i+1} for the system  
c     A_{i+1} x = b_{i+1} is computed using information gathered by the 
c     CG method while solving the previous system  A_i x = b_i. No
c     preconditioner is provided for the first system A_1 x = b_1, and 
c     it is therefore assumed that this system is solved without 
c     preconditioning. 
c
c     A detailed description of the quasi-Newton preconditioners and 
c     numerical results illustrating their performance is given in:
c
c         Morales, JL and Nocedal, J. (1997). Automatic Preconditioning
c         by Limited Memory Quasi-Newton Updating. Technical Report OTC 
c         97/08, Optimization Technology Center, Argonne National 
c         Laboratory and Northwestern University. (See also
c         www.ece.nwu.edu/~nocedal/PSfiles)
c
c
c     2  OVERVIEW OF THE PACKAGE
c     --------------------------
c
c     The package is easily called from any CG routine. To illustrate 
c     how this is done, let us consider the following simplified 
c     description of the CG method applied to the i-th problem in the 
c     sequence (1):
c
c     
c     CG METHOD 
c
c     Compute  r^0 = A_i x^0 - b,  for some initial approximation x^0;
c     define H_0 = I.
c
c     for j=1,2, ...
c
c       1.  compute the preconditioned residual z^{j-1} = H_i r^{j-1}
c
c       2.  check CONVERGENCE
c
c       3.  compute scalars alpha, beta
c
c       4.  compute direction p^j = - z^{j-1} + beta p^{j-1}
c
c       5.  compute x^j = x{j-1} + alpha p^j
c
c       6.  compute r^j = r{j-1} + alpha A_ip^j 
c
c     end for
c     
c     The routines in the package provide the implementation of step 1. 
c     A salient feature of the package is that a single call to the main 
c     routine in the package, at step 1 of the CG method, is used to 
c     perform three different tasks:
c
c       (a)  provide the most recently generated correction pair (s,y) 
c            to PREQN, so that this pair is stored, if appropriate;
c       (b)  compute the preconditioned residual;
c       (c)  inform the package if a new problem in the sequence (1) is
c            beginning to be solved; in this case the old preconditioner 
c            is removed from storage and a new preconditioner based on 
c            the newly collected correction pairs is formed.
c            
c     We illustrate the design and use of the package in the next 
c     paragraphs.   
c
c
c     2.1 Solving the first system.
c     -----------------------------
c
c       We solve the first problem A_1 x = b_1 without any precondi-
c       tioning. At every iteration of this unpreconditioned CG run
c       we provide PRECOND with the new pair of vectors
c
c            s =    p^j  = [ x^{j+1} - x^j ]/alpha
c            y = A_1p^j  = [ r^{j+1} - r^j ]/alpha
c
c       which provide information for constructing the precondi-
c       tioner for the next problem in the sequence. The pair is trans-
c       mitted by means of the call
c
c              call PREQN ( ...  i, j, s, y, r, z  ... ), 
c
c
c        where   i = 1  specifies that the first system is solved
c                j      specifies the iteration of the CG method
c                s      is a vector holding p^{j-1}, the scaled 
c                       difference x^j - x^{j-1}, 
c                y      is a vector holding A_ip^{j-1}, the scaled
c                       difference r^j - r^{j-1}
c                r      is a vector containing the residual r^{j-1}
c                z      on OUTPUT holds the residual  r^{j-1}
c
c       
c       The pair (s, y) will be saved or removed by PRECOND depending
c       on the choice of strategy for collecting correction pairs 
c       specified by the parameter IOP, and depending on the amount of 
c       storage available. (For the first problem the preconditioned
c       residual z is identical to the residual r.)
c
c
c
c     2.2 Building a new preconditioner and computing z^0 = H_i r^0
c     -------------------------------------------------------------
c
c       When the solution of a new problem in the sequence (1) commences 
c       (that is when i > 1 and j = 1) the call to PRECOND instructs it 
c       to remove the previous preconditioner and store in its place a 
c       new preconditioner based on the pairs (s,y) generated while 
c       solving the previous problem i-1. After this is done, the 
c       package computes the initial preconditioned residual z^0 for the 
c       CG iteration. These tasks are performed by a single call to 
c       PRECOND: 
c
c              call preqn ( ...  i, j, s, y, r, z  ... ), 
c
c
c        where   i > 1  specifies that the CG is solving the i-th 
c                       subproblem
c                j = 1  specifies that this is the first iteration of 
c                       the CG method
c                s      is a dummy argument
c                y      is a dummy argument
c                r      is a vector holding the residual r^0 
c                z      on OUTPUT contains the preconditioned initial 
c                       residual H_i r^0 
c          
c
c       We note that s and y are dummy variables in this case, because 
c       the last CG iteration terminated in step 2, and therefore no 
c       new correction pairs were generated since the latest call to 
c       PRECOND which took place the line above, in step 1. 
c
c
c     2.3 Storing pairs and computing z^{j-1} = H_i r^{j-1}.
c     ------------------------------------------------------
c
c       Most of the time, PREQN will be called with i > 1 and j > 1. In
c       this case, PREQN is provided a new correction pair (s,y), and
c       computes the preconditioned residual. These tasks are performed 
c       with the call:
c
c
c              call preqn ( ...  i, j, s, y, r, z  ... ), 
c
c
c        where   i > 1  specifies that the j-th system is solved
c                j > 1  specifies the CG iteration number
c                s      on INPUT is a vector holding p^{j-1}
c                y      is a vector holding the product A_i p^{j-1}
c                r      contains the (j-1)-st residual 
c                z      contains, on output, the (j-1)-st preconditioned 
c                       residual
c
c
c     2.4 Reusing preconditioners.
c     ----------------------------
c
c       In some cases it is appropriate to use the same preconditioner 
c       for several problems. We have added some flexibility to the 
c       package by allowing this option. A preconditioner can be 
c       reused for several problems by setting the logical variable 
c       BUILD to .false. 
c
c       Suppose that MAIN is a user program that calls PREQN repeatedly
c       to solve a sequence of linear systems (1), and suppose that the 
c       user wants to compute the preconditioner only once and use it to
c       solve all the problems in the sequence with i > 2. Thus MAIN 
c       would typically contain a loop, as follows: 
c
c
c       PROGRAM MAIN
c
c       for i=1,..., t
c
c           if i=2 then 
c              BUILD = .TRUE.
c           else
c              BUILD = .FALSE.
c           end if
c
c           call preqn ( ...  i, j, s, y, r, z,  ... , BUILD, ... ), 
c
c       end for 
c
c       This option may be useful, for example, when the matrices A_i
c       are all the same, or when a CG iteration performed very few 
c       iterations.
c
c       In general, we recommend that a new preconditioner be computed 
c       after each linear system is solved. This is done by setting
c                   BUILD = .TRUE.
c       before every call to PREQN. 
c
c
c     3  THE CALLING SEQUENCE AND DESCRIPTION OF THE PARAMETERS
c     ---------------------------------------------------------
c
c     In all cases the package is invoked as:
c
c        call  preqn ( n, m, iop, iprob, jcg, s, y, r, z, w, lw, iw, 
c    *                 liw, build, info, mssg ),
c
c     where the parameters have the following meaning:
c
c
c     N    is an integer variable specifying the dimension of the 
c          problem. An error message is  printed if  N  is zero or 
c          negative. A warning message is printed if  N = 1.
c
c     M    is an integer variable specifying the maximum number of 
c          correction pairs used to build the quasi-Newton precondi-
c          tioner. Values of M in the range [4,20] are recommended. 
c          If PREQN is called with IOP = 1 then  M  must be an even 
c          number. An error message is printed if  M  is less than 2. 
c
c     IOP  is an integer variable specifying the scheme for the 
c          selection of correction pairs.
c
c          IOP = 1  the pairs are selected as a uniform sample  
c                   throughout the CG cycle,
c 
c          IOP = 2  the last  M computed pairs are selected.
c
c          If the uniform sampling scheme is chosen (IOP = 1) an 
c          extra pair may be stored and used in the construction of 
c          the preconditioner, as it is explained now. The first  M 
c          pairs are selected by the uniform sampling scheme. A routine 
c          in the package checks if the last pair produced by the CG 
c          iteration was selected, and if not, saves it. In this case 
c          the number of correction pairs becomes  M + 1. Numerical 
c          experience shows a slight improvement in performance if the 
c          last pair is always included, and since by incorporating it
c          the storage requirements and computational effort increase 
c          only modestly, this strategy has been implemented in the 
c          uniform sampling scheme. The use of  IOP = 1  requires  M  
c          to be an even integer. A warning message is displayed if 
c          this is not the case, then  M  is set to  M - 1.
c
c     IPROB  is an integer variable that specifies the current 
c          problem being solved. 
c
c     JCG  is an integer variable specifying the current CG iteration 
c          number.
c     
c     S    is a real array of dimension N containing the direction p. 
c
c     Y    is a real array of dimension N containing the product Ap.
c          
c     R    is a real array of dimension N containing the residual at
c          iteration JCG of the CG method.
c 
c     Z    is a real array of dimension N. On OUTPUT Z contains the 
c          preconditioned residual.
c
c     W    is a real array of dimension no less than 
c          4(M+1)N + 2(M+1) + 1.  It is used as real workspace, as
c          follows:
c         
c          1) The first 2(M+1)N positions contain the correction pairs
c             generated by the CG iteration on the current problem. See
c             the documentation of subroutine INIPR1 for more details;
c          2) 2(M+1)N + 2(M+1) + 1 positions are reserved to hold the
c             current preconditioner and to assist in the computation
c             of the preconditioned residual. See the comments in
c             subroutine INIPR2. 
c
c          More specifically, the storage demands are:
c
c             pairs (s,y),                  4(M+1)N positions
c             constants rho = 1/y^Ts         (M+1)  positions
c             space for the computation      
c             of products H_i * v            (M+1)  positions
c             scale                            1    position
c             -----------------------------------------------
c                 T O T A L             4(M+1)N + 2(M+1) + 1 positions 
c
c
c     LW   is an integer variable specifying the length of the array W 
c          as declared in the calling subprogram. An error message is 
c          printed if the value of  LW  is less than  
c          4(M+1)N + 2(M+1) + 1.
c
c     IW   is an integer array of dimension no less than 2M + 30. It 
c          is used as integer workspace. A more detailed description
c          of its contents is given an the documentation of subroutines
c          INIPR1 and INIPR2. 
c
c     LIW  is an integer specifying the length of the array  IW  as 
c          declared in the calling subprogram. An error message is 
c          printed if  IW  is less than  2M + 30.
c
c     BUILD is a logical variable that indicates whether the pre-
c          conditioner for the problem  IPROB - 1  should be built. 
c          It should be set to  .FALSE.  before the first call to 
c          PREQN. This variable allows the reuse of a preconditioner 
c          for several problems.  
c
c          BUILD = .TRUE.  a preconditioner is computed using the 
c                   information collected during the CG cycle for 
c                   problem  IPROB - 1.
c
c          BUILD = .FALSE.  the information collected during the 
c                   previous CG cycle is ignored and superseeded 
c                   by the information collected during the current
c                   CG cycle.
c
c     INFO is an integer variable. On output it gives information
c          on the status of the routine.
c
c          INFO =  0   normal termination
c          INFO = -1   error in input parameters
c          INFO = -2   insufficient storage in the workspace arrays
c
c     MSSG is an integer variable specifying the output unit used to 
c          display error messages.
c
c         
c**********************************************************************
c
c     5  OTHER SUBROUTINES CALLED
c     ---------------------------
c
c
c     dcopy, ddot, daxpy from  BLAS.
c
cc
c         
c**********************************************************************
c
c     6  TESTING PREQN
c     ----------------
c
c     Users are encouraged to run the sample calling program 
c     driver1.f  provided with PREQN.
c
c
c***********************************************************************
c
      integer is1, iy1, liperm, llistp, npr, irho, ialph, is2, iy2, 
     *        iscl
c
c     Check parameters iprob and jcg
c
      if ( iprob.lt.1 .or. jcg.lt.1 ) then 
         info = -1
         write(mssg,900)
         return
      end if
c
      if ( iprob.eq.1 ) then
c--------
c        This is the first linear system.
c--------
         if ( jcg.eq.1 ) then
c        ---
c           First CG iteration. Check parameters and initialize 
c        ---
            call chkpar ( n, m, iop, lw, liw, info, mssg )
            if ( info.ne.0 ) then
               return
            else
               call inipr1 ( n, m, iop, iw )
               is1    = iw(5)
               iy1    = iw(6)
               liperm = iw(8)
               llistp = iw(9)
            end if
         else
c        ---
c           Subsequent CG iterations. Store correction pairs.
c        ---
            call storep ( n, m, iop, jcg-1, s, y, w(is1), w(iy1), iw,
     *                    iw(liperm), iw(llistp), mssg )
c
         end if
c
c        No preconditioning ( r --> z ) for the first system
c
         call dcopy ( n, r, 1, z, 1 )
c
      else
c--------
c        Subsequent linear systems
c--------
         if ( jcg.eq.1 ) then
c        ---          
c           First CG iteration. 
c        ---
c           Build the preconditioner with the information collected 
c           during the solution of the previous system iprob-1, or 
c           reuse old preconditioner. If preconditioner is reused, 
c           then just do the cleaning
c           
            if ( build ) then 
               call inipr2 ( n, m, iop, iw )
c
               npr   = iw(11)
               irho  = iw(12)
               ialph = iw(13)
               is2   = iw(14)
               iy2   = iw(15)
               iscl  = iw(16)
c
c              Check that the number of collected pairs is non zero. 
c              Print a warning message if that is the case            
c
               if ( npr.eq.0 ) then
                  write(mssg,910)
                  return
               end if
c
               call buildp ( n, npr, w(is1), w(iy1), w(is2), w(iy2),
     *                       w(irho), w(iscl), iw(liperm) )
            end if
c
c           Clean integer space for the new CG cycle
c
            call cleanp ( iop, iw )

         else
c        ---
c           Subsequent CG iterations. Store the pair (s,y)
c        ---
            call storep ( n, m, iop, jcg-1, s, y, w(is1), w(iy1), iw,
     *                    iw(liperm), iw(llistp), mssg )
         end if
c
c        Compute the matrix vector product z = H*r
c
         call prodct ( n, npr, w(irho), w(ialph), w(is2), w(iy2),
     *                 w(iscl), r, z )
c
      end if
c
 900  format(//,2x,' ** Error in PREQN: incorrect input: iprob or jcg',
     *            /,' are less than 1',/)
 910  format(//,2x,' ** Warning in PREQN: CG collected no pairs ',/,
     *          2x,'    Check the CG routine ',/)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine chkpar ( n, m, iop, lw, liw, info, mssg )
       
      integer  n, m, iop, lw, liw, info, mssg
c
c***********************************************************************
c
c     This routine checks the input parameters and the size of working 
c     space. If an error is detected then INFO is assigned a nonzero 
c     value. The routine also checks the value of M for IOP=1. For this 
c     option M must be an even number; if it is not, the routine changes 
c     it to M-1. Messages are printed in all cases.
c 
c***********************************************************************
c
c     PARAMETERS
c
c     N    is an integer variable specifying the dimension of the 
c          problem
c
c     M    is an integer variable specifying the maximum number of 
c          correction pairs used to build the preconditioner. The 
c          routine allows a minimum of two pairs (M=2)
c
c     IOP  is an integer variable specifying the storing scheme for 
c          the correction pairs {s,y}.
c
c     LW   is an integer variable specifying the length of the array 
c          W  as declared in the calling subprogram.
c
c     LIW  is an integer specifying the length of the array IW as 
c          declared in the calling subprogram.
c
c     INFO is an integer variable. On output it gives information
c          on the status of the routine.
c
c          INFO =  0   normal termination
c          INFO = -1   error in input parameters
c          INFO = -2   insufficient storage in arrays W, IW
c
c     MSSG is a positive integer variable specifying the unit at 
c          which  error messages will be displayed
c
c**********************************************************************
c
      integer maux
c
      info = 0
      if ( n.lt.1 .or. m.lt.2 .or. iop.lt.1 .or. iop.gt.2 ) then 
         info = -1 
         write(mssg,900) 
         return
      else      
         if ( lw.lt. (m+1)*(4*n+2) + 1 .or. liw.lt. 2*m + 30 ) then 
            info = -2
            write(mssg,920)
            return
         end if
      end if
c
c     Check that the maximum number of pairs is an even number
c     in the sampling option (IOP=1)
c
      if ( iop.eq.1 ) then
         maux = m
         maux = 2*(maux/2)
         if ( maux.lt.m ) write(mssg,940) 
         m = maux
      end if
c
c     Check if n = 1. Print a warning message
c
      if ( n.eq.1 )  write(mssg,950)
c

 900  format(//,2x,' ** Error in PREQN: incorrect input parameters',/)
 920  format(//,2x,' ** Error in PREQN: insufficient storage',/)
 940  format(//,2x,' ** Warning in PREQN: m is not an even number',/
     *         ,2x,'    m-1 pairs will be stored',/)
 950  format(//,2x,' ** Warning in PREQN: n is 1',/) 
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine inipr1 ( n, m, iop, iw )
       
      integer  n, m, iop, iw(*)
c
c***********************************************************************
c
c     This routine is called only once, at the first CG iteration for 
c     the first linear system in the sequence. The routine initializes
c     the portions of the arrays IW and W that store correction pairs 
c     generated by the CG iteration, as well as pointers. The new 
c     correction pairs are stored in arrays called S1 and Y1 in other
c     routines
c
c     DESCRIPTION
c
c     iw(1)  : NPR     CG iteration counter used in several routines
c     iw(2)  : COUNT   a variable needed by routine STOREP when IOP=1
c     iw(3)  : ICYCLE    "        "        "     "     "    "
c     iw(4)  : M + 1   the maximum number of pairs when IOP=1.
c     iw(5)  : pointer to S1 in array W. 
c     iw(6)  : pointer to Y1 in array W. 
c     iw(7)  : it is an offset set to 28
c     iw(8)  : pointer to IPERM in IW. IPERM is a vector of pointers
c              keeping the chronological order of the arrays. For
c              example, if IPERM = [1 4 3 2], then the oldest pair is
c              in the first position, and the newest pair in the second
c              position.
c              
c     iw(9)  : pointer to LISTP in  IW. LISTP is the list of indices of 
c              pairs selected within the CG cycle.  For example, if  
c              LISTP = [2 9 17 25 ] then the first pair was generated 
c              during the second CG iteration, and so on.
c     iw(10) : empty 
c
c***********************************************************************
c
c     PARAMETERS
c
c     N    is an integer variable specifying the dimension of the 
c          problem
c
c     M    is an integer variable specifying the maximum number of 
c          pairs used to build the quasi-Newton preconditioner.
c
c     IOP  is an integer variable specifying the storing scheme of 
c          the  pairs {s,y}.
c
c     W    is a real array of dimension no less than 
c          4(M+1)N + 2(M+1) + 1.  It is used as real workspace.
c
c     IW   is an integer array of dimension no less than 2M + 30. 
c          It is used as integer workspace.
c
c**********************************************************************
c
      integer mp1, mn, ioff, i
c
      mp1  = m + 1
c
c     Set counters and pointers 
c
      iw(1) = 0
      if ( iop.eq.1 ) then
         iw(2) = 1
         iw(3) = 1
      end if
      iw(4)  = mp1
c
c     Set pointers for s1, y1, (the arrays that contain the new 
c     correction pairs) and for iperm, listp
c
      mn    = mp1*n
      iw(5) = 1
      iw(6) = iw(5) + mn     
c
      ioff  = 28
      iw(7) = ioff
      iw(8) = iw(7) + 1
      iw(9) = iw(8) + mp1
c
c     Initialize space in iw from i=10 to 30 + 2M.
c
      do 10 i=10, 30 + 2*m
         iw(i) = 0
 10   continue
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine inipr2 ( n, m, iop, iw )
       
      integer  n, m, iop, iw(*)
c
c***********************************************************************
c
c     This routine is called at the beginning of every new linear
c     system (except for the first one). It reorganizes the information 
c     gathered by the CG method during the solution of the previous
c     problem. O output, the array IW holds the following information:
c
c     DESCRIPTION
c
c     iw(11) : the effective number (NPR) of pairs at the completion of
c              the CG iteration. This is the number of pairs in the new 
c              preconditioner.
c
c     iw(12) : pointer to RHO in array W              
c     iw(13) : pointer to ALPHA in array W.
c     iw(14) : pointer to S2 in arrays W.
c     iw(15) : pointer to Y2 in array W. 
c     iw(16) : pointer to GAMMA in array W
c
c     S2, Y2, RHO, ALPHA, and GAMMA are required to build the new 
c     preconditioner
c
c***********************************************************************
c
c     PARAMETERS
c
c
c     N    is an integer variable specifying dimension of the problem
c
c     M    is an integer variable specifying maximum number of pairs
c          that must be saved to build the quasi-Newton preconditioner.
c
c     IOP  is an integer variable specifying the storing scheme of the 
c          pairs {s,y}.
c
c     IW   is an integer array of dimension no less than 2M + 30. It 
c          is used as integer workspace.
c
c**********************************************************************
c
      integer npr, mp1, nprn, listp
c
c     Determine the number of pairs npr in the new preconditioner. 
c     (Recall that the last correction pair generated by the CG 
c     iteration is added, if necessary, when IOP = 1.) 
c
      npr   = iw(1)
      mp1   = iw(4)
c     ioff  = iw(7)
      listp = iw(9) 
c
c     For IOP = 1 test if the last pair is already stored. If the
c     number of stored pairs is not greater than m, then we know that
c     this last correction pair has already been stored.
   
      if ( iop.eq.1 ) then
         if ( npr.gt.m ) then
            if ( iw(listp + mp1 - 1) .ne. iw(listp + mp1 - 2) ) then
               npr = mp1
            else 
               npr = m
            end if
         end if
      else
         if ( npr.ge.m ) npr = m
      end if
c
      nprn = npr*n
c
c     Space for rho, alpha, s2, y2 and gamma
c  
      iw(11) = npr
      iw(12) = 2*mp1*n + 1
      iw(13) = iw(12)  + npr
      iw(14) = iw(13)  + npr
      iw(15) = iw(14)  + nprn
      iw(16) = iw(15)  + nprn
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine cleanp ( iop, iw ) 
       
      integer  iop, iw(*)
c
c***********************************************************************
c
c     This routine resets the working space. It is called after building 
c     the preconditioner by means of routine BUILDP
c
c***********************************************************************
c
c     PARAMETERS
c
c     IOP  is an integer variable specifying the storing scheme of the 
c          pairs {s,y}.
c
c     IW   is an integer array of dimension no less than 2M + 30. It 
c          is used as integer workspace.
c
c***********************************************************************
c
      integer i, mp1, ioff
c
      mp1  = iw(4)
      ioff = iw(7)
c
c     Reset the space
c
      iw(1) = 0
c      
      if ( iop.eq.1 ) then
         iw(2) = 1
         iw(3) = 1
      end if
c
      do 10 i=1, 2*mp1
         iw(ioff+i) = 0
 10   continue
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine storep ( n, m, iop, jcg, s, y, s1, y1, iw, iperm,
     *                    listp, mssg )
       
      integer  n, m, iop, jcg, iw(*), iperm(*), listp(*), mssg
      double   precision s(n), y(n), s1(*), y1(*)
c
c***********************************************************************
c
c     This routine is called at every CG iteration, except for the 
c     first iteration for each problem. The routine receives a new
c     correction pair (s,y) and determines if it is to be stored,
c     according to the amount of storage used so far and following
c     the storage scheme determined by the parameter IOP. The arrays
c     S1 and Y1 are used to store the new correction pairs, and hold
c     at most M of them. If the new correction pair is to be stored and 
c     the arrays S1 and 11 are full, then the routine determines which 
c     pair has to be deleted from S1 and Y1 to make room for the 
c     accepted pair. The routine keeps a record of the pairs that are 
c     being held. 
c
c     PARAMETERS
c
c     N      is an integer specifying the dimension of the problem
c
c     M      is an integer specifying  the maximum number of pairs
c            in the preconditioner
c
c     IOP    is an integer variable, chosen by the user, indicating the
c            scheme used to select correction pairs
c
c     JCG    is an integer variable that holds the CG iteration
c            at which pair (s,y) was collected
c
c     S      double precision array of dimension N storing the new 
c            CG  direction  p_jcg
c
c     Y      double precision array of dimension N storing the new
c            product  A*p_jcg
c
c     S1     is a double precision array of dimension MN that holds
c            all the collected directions p
c
c     Y1     is a double precision array of dimension MN that holds
c            all the collected products A*p
c
c     LISTP  is an integer array of dimension  M+1  that holds the 
c            iteration of the CG method at which the pair (s,y) was 
c            collected
c
c     IPERM  is an array of pointers of dimension M. IPERM(i) is the 
c            position in the arrays S1, Y1 in which the i-th pair is 
c            stored. 
c
c     MSSG   is an integer specifying the output unit in which errors 
c            will be displayed
c
c***********************************************************************
c
      integer mp1, m2, count, cycle, istrt, jl, i, ii, iold
c
c     Set counters
c
      mp1  = iw(4)
c
c     For IOP=1 store the new pair also in position m+1. 
c    
      if ( iop.eq.1 ) then
         iperm(mp1) = mp1
         listp(mp1) = jcg
c
         call dcopy ( n, s, 1, s1(m*n+1), 1 )
         call dcopy ( n, y, 1, y1(m*n+1), 1 )
      end if
c
c     When the storage is not full we can keep the chronological
c     order of the pairs and there is no need to remove pairs.
c
      if ( jcg.le.m ) then
c
c        Save the pair in position j = jcg.
c
         listp(jcg) = jcg
         iperm(jcg) = jcg
         istrt    = (jcg-1)*n + 1
         call dcopy ( n, s, 1, s1(istrt), 1 )
         call dcopy ( n, y, 1, y1(istrt), 1 )
c
c        Update the variable that counts how many pairs have been stored 
c
         iw(1) = iw(1) + 1
         return
c
      else
c
c        Check the entering pair. If the pair is accepted then
c        compute index of leaving pair
c
         if ( iop.eq.1 ) then
            m2    = m/2 
            count = iw(2)
            cycle = iw(3)
            jl = 1 + (m2 + count - 1)*2**cycle 
c
            if ( jcg.eq.jl ) then
               jl = 1 + (2*count - 1)*2**(cycle - 1)
            else
               return
            end if
         else
            jl = jcg - m 
         end if
c
c        Store/delete information. First find the leaving pair in listp.
c        Then move information backwards to make room in the last
c        position.
c
         i = 1
 10      if ( listp(i).eq.jl ) then
            iold = iperm(i)
c
c           Make room in the last position
c
            do 20 ii=i, m-1
               iperm(ii) = iperm(ii+1)
               listp(ii) = listp(ii+1)
 20         continue
c
c           Add the index of the new pair to LISTP. Create the 
c           pointer for the new pair. Insert the new pair. 
c
            listp(m) = jcg
            iperm(m) = iold
c
            istrt = (iold-1)*n + 1
            call dcopy ( n, s, 1, s1(istrt), 1 ) 
            call dcopy ( n, y, 1, y1(istrt), 1 )
         else
            i = i + 1
c
c        This condition may happen just when the data structures
c        are corrupted.
c
            if ( i.gt.m ) then 
               write(mssg,*) ' We could not find the leaving pair'
               return
            end if
            go to 10
         end if
c
      end if
c
c     Update counters and leave
c
      iw(1) = iw(1) + 1
c
      if ( iop.eq.1 ) then
         count = count + 1
         if ( count.gt.m2 ) then
            count = 1
            cycle = cycle + 1
         end if
         iw(2) = count
         iw(3) = cycle
      end if
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine prodct ( n, npr, rho, alpha, s2, y2, scale, r, z )
       
      integer  n, npr
      double   precision  rho(*), alpha(*), s2(*), y2(*), scale, r(n),
     *                    z(n)
c
c***********************************************************************
c
c     This routine computes the product H times an arbitrary vector v,
c     where H is a limited memory preconditioner stored in the form
c     of m pairs (s,y). The product is computed in 4MN floating point
c     operations.
c
c***********************************************************************
c
c     PARAMETERS
c
c     N      is an integer variable specifying dimension of 
c            the problem
c
c     NPR    is an integer specifying number of pairs in the 
c            preconditioner. If NPR = 0 the routine returns
c            in Z the input array R. (See the description of Z).
c
c     RHO    is a double precision array that holds the NPR 
c            quantities 1/y^Ts
c
c     ALPHA  is a double precision array that is used to hold 
c            intermediate computations
c  
c     S2     is a double precision array that stores the S's
c
c     Y2     is a double precision array that stores the Y's
c
c     SCALE  is a double precision variable specifying scale
c
c     R      is a double precision array holding the residual
c
c     Z      is a double precision array holding the preconditioned 
c            residual  H*R
c
c*********************************************************************
c
      integer  j, istrt 
      double   precision beta, almbet, alphj, ddot
c
c     Solve the unpreconditioned case (NPR = 0) and leave
c
      call dcopy ( n, r, 1, z, 1 )
      if ( npr.eq.0 )                 return
c
c     First sweep 
c
      do 10 j = npr, 1, -1
         istrt    = (j-1)*n + 1
         alpha(j) = rho(j)*ddot ( n, s2(istrt), 1, z, 1 )     
c
c        Update vector z
c
         alphj = - alpha(j)
         call daxpy ( n, alphj, y2(istrt), 1, z, 1 )
 10   continue
c
c     Now get the initial Hessian approximation H(0). 
c     It is the identity matrix scaled by the last pair. 
c
      call dscal ( n, scale, z, 1 )
c
c     Second sweep. 
c 
      do 30 j=1, npr, 1
         istrt  = (j-1)*n + 1
         beta   = rho(j)*ddot ( n, y2(istrt), 1, z, 1 )
         almbet = alpha(j) - beta
         call daxpy ( n, almbet, s2(istrt), 1, z, 1 )
 30   continue
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine buildp ( n, npr, s1, y1, s2, y2, rho, scale, iperm )
       
      integer  n, npr, iperm(*)
      double   precision s1(*), y1(*), s2(*), y2(*), rho(*), scale
c
c***********************************************************************
c     
c     This routine computes a quasi-Newton preconditioner using  
c     NPR pairs collected by the CG method. First the pairs stored
c     in arrays the S1, Y1 are moved to the arrays S2, Y2. This 
c     involves reordering the pairs according to their chronological 
c     order. The order is recorded in the vector IPERM. Then the 
c     scalars RHO(i) and GAMMA are computed and stored.
c
c***********************************************************************
c
c     PARAMETERS
c
c     N      is an integer variable specifying the dimension of the 
c            problem
c
c     NPR    is an integer specifying the number of pairs in the 
c            preconditioner
c
c     RHO    is a double precision array that holds the NPR quantities
c            1/y^Ts
c            
c
c     ALPHA  is a double precision array that is used to hold 
c            intermediate computations
c  
c     S1     is a double precision array that stores the directions p
c            collected by the CG method. 
c
c     Y1     is a double precision array that stores the products Ap
c            collected by the CG method. 
c
c     S2     is a double precision array that stores the vectors p;
c            they are now stored in chronological order
c
c     Y2     is a double precision array that stores the products Ap
c            in chronological order. 
c
c     RHO    is a double precision array that holds the NPR quantities
c            1/y^Ts
c
c     SCALE  is a double precision variable holding the scale
c
c     IPERM  is an integer array with pointers to the correction pairs. 
c            IPERM(i) gives the position of pair (s_i,y_i) in the arrays 
c            S1, Y1
c
c**************************************************************************
c
      integer istr1, istr2, istrt, j, ipj
      double  precision one, skyk, yknrm, ddot, dnrm2
      data    one/1.0d0/
c
c     Reallocate the pairs
c
      do 10 j=1, npr
         ipj   =  iperm(j)
         istr1 =  n*(ipj-1) + 1
         istr2 =  n*(  j-1) + 1 
         call dcopy ( n, s1(istr1), 1, s2(istr2), 1 )
         call dcopy ( n, y1(istr1), 1, y2(istr2), 1 )
 10   continue
c
c     Compute  rho(j)
c
      do 20 j=1, npr
         istrt  = n*(j-1) + 1
         skyk   = ddot ( n, s2(istrt), 1, y2(istrt), 1 )
         rho(j) = one/skyk
 20   continue
c
c     Compute the scale using the last pair
c
      istrt = n*(npr-1) + 1
      yknrm = dnrm2 ( n, y2(istrt), 1 )
      yknrm = yknrm*yknrm
c
      scale = skyk/yknrm
c
      return
      end














c
c     The routines pcg.f and pcg1.f implement the Preconditioned 
c     Conjugate Gradient Method (PCG) to solve a symmetric and positive
c     definite system
c   
c                         A*x = b.
c
c     They show how the PREQN package is called from within a
c     Conjugate Gradient code. A user will normally replace these two 
c     routines by his/her own implementation of the PCG method. 
c
c-----------------------------------------------------------------------
c        
      subroutine pcg ( n, nz, a, adiag, jptra, indra, x, b, anorm, 
     *                 w,
     *                 m, iop, iprob, wp, lwp, iwp, liwp, 
     *                 build, info, 
     *                 maxit, tolcg, iout, mssg ) 

      integer  n, nz, jptra(*), indra(*), m, iop, iprob, lwp, iwp(*), 
     *         liwp, info, maxit, iout, mssg
      double   precision a(nz), adiag(n), x(n), b(n), anorm,
     *         w(4*n), wp(*), tolcg
      logical  build
c
c***********************************************************************
c  
c     This routine partitions the workspace and calls pcg1.f, which
c     implements that PCG method. For a description of the parameters,
c     see routine pcg1.f below.
c
c
c***********************************************************************
c        
      integer  l1, l2, l3, l4
c
c     Memory partition
c
      l1 = 1
      l2 = l1 + n
      l3 = l2 + n
      l4 = l3 + n
c
      call pcg1 ( n, nz, a, adiag, jptra, indra, x, b, anorm, 
     *            w(l1), w(l2), w(l3), w(l4),
     *            m, iop, iprob, wp, lwp, iwp, liwp, 
     *            build, info, 
     *            maxit, tolcg, iout, mssg )
c      
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine pcg1   ( n, nz, a, adiag, jptra, indra, x, b, anorm,
     *                    r, v, gv, pr, 
     *                    m, iop, iprob, wp, lwp, iwp, liwp, 
     *                    build, info, 
     *                    maxit, tolcg, iout, mssg ) 
c
      integer  n, nz, jptra(*), indra(*), m, iop, iprob, lwp, iwp(*), 
     *         liwp, info, maxit, iout, mssg
      double   precision a(nz), adiag(n), x(n), b(n), anorm,
     *         r(n), v(n), gv(n), pr(n), wp(*), tolcg
      logical  build

c
c***********************************************************************
c 
c    This routine implements the PCG method. The preconditioner is
c    computed by the PREQN package
c
c***********************************************************************
c 
c
c     DESCRIPTION OF THE PARAMETERS
c     -----------------------------
c
c     N    is an integer variable specifying the dimension of A
c
c     NZ   is an integer variable specifying the number of nonzeros
c          in A
c     A    is a real array holding the nonzero elements of A in sparse
c          column oriented format. Given that A is symmetric, only its 
c          lower triangle is stored.
c
c     ADIAG  is a real array of dimension NZ containing the diagonal 
c          elements of matrix A.
c
c     JPTRA  integer array of dimension N+1 holding pointers to the 
c          columns of A in the column oriented format.
c
c     INDRA  integer array of dimension NZ containing row indices of 
c          the elements in A.
c
c     X    real array of dimension N. On input it contains an initial 
c          approximation. On output it contains the last approximation 
c          found by the CG method.
c
c     B    real array of dimension N holding the elements of the 
c          right hand side
c
c     ANORM  real variable holding the infinity norm of the matrix A. 
c          This quantity is used in the termination test.
c
c     W    double precision array of length 4N. It is used as 
c          real work space to hold intermediate computations.
c
c     M    is an integer variable specifying the maximum number of pairs
c          used in the quasi-Newton preconditioner.
c
c     IOP  is an integer variable specifying the storing scheme of the 
c          pairs {sk,yk}.
c
c          IOP = 1  the pairs are sampled in order to have
c                   an almost uniform distribution.
c
c          IOP = 2  the last M consecutive pairs are stored.
c
c     IPROB  is an integer labeling the linear system being solved 
c
c     WP   is a real array of dimension no less than 
c          4(M+1)N + 2(M+1) + 1.  It is used as real workspace to
c          compute the quasi-Newton preconditioner.
c
c     LWP  is an integer variable specifying length of the array W as 
c          declared in the calling program.
c
c     IWP  is an integer array of dimension no less than 2M + 30. It 
c          is used as integer workspace to compute the quasi-Newton 
c          preconditioner.
c
c     LIWP is an integer specifying the length of the array IW as 
c          declared in the calling subprogram.
c
c     BUILD is a logical variable that determines if the preconditioner
c          stored while solving IPROB-1 should be built. This variable
c          allows the reuse of a preconditioner.
c
c     INFO  is an integer variable. On output it gives information
c          about the status of the routine.
c
c          INFO =  0   normal termination
c          INFO = -1   error in input parameters
c          INFO = -2   insufficient storage in the workspace arrays
c          INFO = -4   initial residual is sufficiently small
c         
c    MAXIT integer variable specifying the maximum number of CG 
c          iterations
c
c    TOLCG Tolerance to stop the CG iteration.
c
c    IOUT  is an integer variable that controls the amount and type of
c          output
c
c          IOUT = 0   no output printed
c          IOUT = >0  output at every IOUT iterations of the CG method, 
c                    plus output at first and last iteration.
c
c    MSSG  is an integer variable specifying the output unit used 
c          to display error messages
c
c
c***********************************************************************
c

      integer  jcg, i
      double   precision dnrmif, bnorm, xnorm, rnorm, resrel, 
     *                   rho, ddot, beta, rhold, vgv, alpha
c
c     Check input 
c
      info =  0

      if ( n.le.0 .or. m.lt.0 .or. maxit .le.0 ) then 
         write(mssg,900)
         info = -1
         return
      end if
c
c     Initial printing
c
      if ( iout.ne.0 ) write(6,920) m, iop, iprob
c     
c     Compute the initial residual. First compute the sparse product
c     Ax = v using the routine dsymvs from the MINPACK-2 package.
c
      call dsymvs ( n, a, adiag, jptra, indra, x, v )
      do 100 i = 1, n
         r(i)  = b(i) - v(i)
 100  continue
c
c     if the initial residual is sufficiently small, terminate
c
      bnorm = dnrmif ( n, b )
      rnorm = dnrmif ( n, r )
      if ( rnorm.le.bnorm*tolcg ) then
         info = -4
         return
      end if
c
c*************************
c     Begin CG iteration
c*************************
c
      do 300 jcg = 1, maxit
c
c        Compute the preconditioned residual using the PREQN package
c
         call preqn ( n, m, iop, iprob, jcg, v, gv, r, pr, wp, lwp,
     *                iwp, liwp, build, info, mssg )
         if ( info.ne.0 ) then
            write(mssg,940)
c
            return
         end if 
c--------
c        Check convergence. Print the relative residual at the initial
c        and final iterations, as well as at selected iterations.
c        The stopping test is  
c
c               ||r|| < tolcg*( ||A|| ||x|| + ||b|| )
c
c        where ||*|| denotes the infinity norm.         
c
         xnorm  = dnrmif ( n, x )
         rnorm  = dnrmif ( n, r )
         resrel = rnorm /( anorm*xnorm + bnorm )
c
         if ( iout.ne.0 ) then
            if ( jcg.eq.1 ) then 
               write(6,960) jcg-1, resrel
            else
               if ( mod(jcg-1,iout).eq.0 ) write(6,960) jcg-1, resrel
            end if
         end if
c
         if ( resrel.le.tolcg ) then 
            if ( iout.ne.0 ) write(6,960) jcg-1, resrel
            return
         end if
c--------
c
c        Compute a new direction v       
c
         rho = ddot ( n, r, 1, pr, 1 )
c
         if ( jcg.eq.1 ) then
            call dcopy ( n, pr, 1, v, 1 )
         else
            beta = rho/rhold
c
            do 200 i=1, n
               v(i) = pr(i) + beta*v(i)
 200        continue
         end if
c
c     Compute the matrix-vector product Av = gv
c
         call dsymvs ( n, a, adiag, jptra, indra, v, gv )
c
         vgv = ddot ( n, v, 1, gv, 1 )
c
c     Compute step length, new iterate and new residual
c
         alpha = rho / vgv
         call daxpy ( n,  alpha,  v, 1, x, 1 )
         call daxpy ( n, -alpha, gv, 1, r, 1 )        
c
         rhold = rho
 300  continue
c
c**************************
c     end of CG iteration
c**************************
c
 900     format(//,2x,' * Error in PCG: incorrect input parameters',/)
 920     format(////,15x,'*** New Problem ***',/, 2x,
     *           'Preconditioned Conjugate Gradient method using ',/,
     *       2x, 'a quasi-Newton preconditioner ', //,
     *       5x, 'number of pairs in the preconditioner = ', i3,/,
     *       5x, 'sampling strategy                     = ', i3,/,
     *       5x, 'number of problem                     = ', i3,//,
     *       5x,'||rel|| is the norm of the relative residual',/)
 940     format(//,2x,' * Abnormal termination in PCG.', 
     *                ' Error in PREQN',/)
 960     format( 5x,' iter = ',i6, 2x,' ||rel|| = ',
     *           e14.8)

      return
      end


      
     






















c                                                                      
c TENSOLVE:  A Software Package for Solving Systems of Nonlinear
c            Equations and Nonlinear Least Squares Problems Using
c            Tensor Methods.
c
c AUTHORS:   Ali Bouaricha                                         
c            Argonne National Laboratory   
c            MCS Division                                               
c            e-mail: bouarich@mcs.anl.gov                               
c AND                                                                   
c            Robert B. Schnabel                                         
c            University of colorado at Boulder                          
c            Department of computer Science                             
c            e-mail: bobby@cs.colorado.edu                              
c                                                                        
c DATE:      Version of January, 1997                                         
c                                                                        
c Purpose of Tensolve:
c                                                                        
c       TENSOLVE finds roots of systems of n nonlinear equations in n
c       unknowns, or minimizers of the sum of squares of m > n nonlinear
c       equations in n unknowns.  It allows the user to choose between
c       a tensor method based on a quadratic model and a standard method
c       based on a linear model. Both models calculate an analytic or
c       finite-difference Jacobian matrix at each iteration.  The tensor
c       method augments the linear model with a low-rank, second-order
c       term that is chosen so that the model is hardly more expensive
c       to form, store, or solve than the standard linear model.  Either
c       a line search or trust-region step-selection strategy may be 
c       selected.  The tensor method is significantly more efficient 
c       than the standard method in iterations, function evaluations, and
c       time.  It is especially useful on problems where the Jacobian at
c       the solution is singular.
c       The software can be called with an interface where the user 
c       supplies only the function, number of nonlinear equations, number 
c       of variables, and starting point; default choices are made for 
c       all the other input parameters.  An alternative interface allows 
c       the user to specify any input parameters that are different from 
c       the defaults.
c       
c List of subroutine and function names called by TENSOLVE:
c
c TS2DTR,TSBSLV,TSCHKI,TSCHKJ,TSCPMU,TSCPSS,TSDLOD,TSD1SV,TSDFCN,TSDFLT,
c TSDUMJ,TSFAFA,TSFDFJ,TSFRMT,TSFSCL,TSFSLV,TSJMUV,TSJQTP,TSLMIN,TSLMSP,
c TSLSCH,TSMAFA,TSMDLS,TSMFDA,TSMFDV,TSMGSA,TSMSDA,TSMSDV,TSMSLV,TSNECI,
c TSNESI,TSNESV,TSNSTP,TSPRMV,TSRSLT,TSQ1P1,TSQFCN,TSQLFC,TSQMLV,TSQMTS,
c TSQMUV,TSQRFC,TSRSID,TSSCLF,TSSCLJ,TSSCLX,TSSLCT,TSSMIN,TSSMRD,TSSQP1,
c TSSSTP,TSSTMX,TSTRUD,TSUDQV,TSUNSF,TSUNSX,TSUPSF,TSUTMD.
c
c Packages called by TENSOLVE:
c
c UNCMIN (R. B. Schnabel, J. E. Koontz, and B. E. Weiss, 
c "A Modular System of Algorithms of Unconstrained Minimization", 
c Acm Trans. Math. Softw., 11 (1985), 419-440).
c
c BLAS called by TENSOLVE:
c
c LEVEL 1 BLAS: DASUM,DAXPY,DCOPY,DDOT,DNRM2,DSCAL,DSWAP,IDAMAX
c LEVEL 2 BLAS: DGEMV
c
c Parameters and Default Values for the interfaces TSNECI and TSNESI:
c ==================================================================
c
c Following each variable name in the list below appears a one- or
c two-headed arrow symbol of the form ->, <-, and <-->.
c These symbols signify that the variable is for input, output, and 
c input-output, respectively.
c The symbol EPSM in some parts of this section designates the machine 
c precision.

c MAXM->: A positive integer specifying the row dimension of the work 
c array WRKNEM in the user's calling program. It must satisfy 
c MAXM >= M+N+2. The provision of MAXM, MAXN, and MAXP allows 
c the user the flexibility of solving several problems with different 
c values of M and N one after the other, with the same work arrays.

c MAXN->: A positive integer specifying the row dimension of the work 
c array WRKNEN in the user's calling program. It must satisfy 
c MAXN >= N+2.  

c MAXP->: A positive integer specifying the row dimension of the work 
c array WRKUNC in the user's calling program.  It must satisfy 
c MAXP >= NINT(sqrt(N)), where NINT is a function that rounds to the 
c nearest integer.

c X0->: An array of length N that contains an initial estimate 
c of the solution x*. 

c M->: A positive integer specifying the number of nonlinear equations.

c N->: A positive integer specifying the number of variables in the 
c problem.  

c TYPX->:  An array of length N in which the typical size of the 
C components of X is specified. The typical component sizes should be 
c positive real scalars. If a negative value is specified, its absolute 
c value will be used. If 0.0 is specified, 1.0 will be used. This 
c vector is used to determine the scaling matrix, Dx. Although the 
c package may work reasonably well in a large number of instances without 
c scaling, it may fail when the components of x* are of radically 
c different magnitude and scaling is not invoked. If the sizes 
c of the parameters are known to differ by many orders of magnitude, then 
c the scale vector TYPX should definitely be used. For example, if 
c it is anticipated that the range of values for the iterates xk would be

c                   x1 in [-10e+10,10e+10]
c                   x2 in [-10e+2,10e+4]
c                   x3 in [-6*10e-6,9*10e-6]

c then an appropriate choice would be TYPX = (1.0e+10,1.0e+3,7.0e-6). 
c Module TSDFLT returns TYPX = (1.0,...,1.0).  

c TYPF->: An array of length M in which the typical size of the components 
c of F is specified. The typical component sizes should be positive real 
c scalars.  If a negative value is specified, its absolute value will be 
c used. If 0.0 is specified, 1.0 will be used. This vector is used to 
c determine the scaling matrix DF. TYPF should be chosen so that all 
c the components of DF(x) have similar typical magnitudes at points not 
c too near a root, and should be chosen in conjunction with FTOL.  It is 
c important to supply values of TYPF when the magnitudes of the components 
c of F(x) are expected to be very different.  If the magnitudes of the 
c components of F(x) are similar, the choice DF = I suffices.  Module 
c TSDFLT returns TYPF = (1.0,...,1.0).  

c ITNLIM->:  Positive integer specifying the maximum number of iterations 
c to be performed before the program is terminated. Module TSDFLT returns 
c ITNLIM = 150. If the user specifies ITNLIM <= 0, the module TSCHKI will 
c supply the value 150.

c JACFLG->: Integer designating whether or not an analytic Jacobian has 
c been supplied by the user.
c JACFLG = 0 : No analytic Jacobian supplied.  The Jacobian is obtained 
c by finite differences. 
c JACFLG = 1 : Analytic Jacobian supplied. 
c The module TSDFLT returns the value 0.  If the user specifies an illegal 
c value, the module TSCHKI will supply the value 0.  

c GRADTL->: Positive scalar giving the tolerance at which the scaled 
c gradient of f(x) = 0.5*F(x)-trans*F(x) is considered close enough to 
c zero to terminate the algorithm. The scaled gradient is a measure of 
c the relative change in F in each direction xj divided by the relative 
c change in xj. The module TSDFLT returns the value EPSM**(1/3).  If the 
c user specifies a negative value, the module TSCHKI will supply 
c the value EPSM**(1/3).

c STEPTL->: A positive scalar providing the minimum allowable relative 
c step length. STEPTL should be at least as small as 10**(-d), where d 
c is the number of accurate digits the user desires in the solution x*. 
c The program may terminate prematurely if STEPTL is too large.  Module 
c TSDFLT returns the value EPSM**(2/3).  If the user specifies a negative 
c value, the module TSCHKI will supply the value EPSM**(2/3).

c FTOL->: A positive scalar giving the tolerance at which the scaled 
c function DF*F(x) is considered close enough to zero to terminate the 
c algorithm. The program is halted if ||DF*F(x)|| (in the infinity norm) 
c is <= FTOL. This is the primary stopping condition for nonlinear 
c equations; the values of TYPF and FTOL should be chosen so that this 
c test reflects the user's idea of what constitutes a solution to the 
c problem. The module TSDFLT returns the value EPSM**(2/3). If the 
c user specifies a negative value, the module TSCHKI will supply the 
c value EPSM**(2/3).

c METHOD->: An integer designating which method to use. 
c METHOD = 0 : Newton or Gauss-Newton algorithm is used. 
c METHOD = 1 : Tensor algorithm is used.  
c Module TSDFLT returns value 1. If the user specifies an illegal value, 
c module TSCHKI will reset METHOD to 1.

c GLOBAL->: An integer designating which global strategy to use.
c GLOBAL = 0 : Line search is used.
c GLOBAL = 1 : Two-dimensional trust region is used. 
c Module TSDFLT returns value of 0. If the user specifies an illegal 
c value, module TSCHKI will reset GLOBAL to 0.

c STEPMX->: A positive scalar providing the maximum allowable scaled step 
c length ||Dx*(x+ - xc)||2, where Dx = diag(1/TYPX_j). STEPMX is used to 
c prevent steps that would cause the nonlinear equations problem to 
c overflow, and to prevent the algorithm from leaving the area of 
c interest in parameter space.  STEPMX should be chosen small enough 
c to prevent these occurrences but should be larger than any anticipated 
c "reasonable" step. Module TSDFLT returns the value STEPMX = 10e+3.
c If the user specifies a nonpositive value, module TSCHKI sets STEPMX 
c to 10e+3.

c DLT->: A positive scalar giving the initial trust region radius. When
c the line search strategy is used, this parameter is ignored. For the 
c trust region algorithm, if DLT is supplied, its value should reflect 
c what the user considers a maximum reasonable scaled step length at 
c the first iteration. If DLT = -1.0, the routine uses the length of 
c the Cauchy step at the initial iterate instead. The module TSDFLT 
c returns the value -1.0. If the user specifies a nonpositive value, 
c module TSCHKI sets DLT = -1.0.

c IPR->: The unit on which the package outputs information. TSDFLT returns 
c the value 6. 

c WRKUNC->: Workspace used by UNCMIN. The user must declare this
c array to have dimensions MAXP*LUNC in the calling routine.

c LUNC->: A positive integer specifying the column dimension of the work 
c array WRKUNC in the user's calling program. It must satisfy 
c LUNC >= 2*NINT(sqrt(N))+4.

c WRKNEM->: Workspace used to store the Jacobian matrix, the function
c values matrix FV, the tensor matrix ANLS, and working vectors. The 
c user must declare this array to have dimensions MAXM*LNEM in the 
c calling routine.  

c LNEM->: A positive integer specifying the column dimension of the work 
c array WRKNEM in the user's calling program. It must satisfy 
c LNEM >= N+2*NINT(sqrt(N))+11.

c WRKNEN->: Workspace used to store the matrix S of previous 
c directions, the matrix SHAT of linearly independent directions, and 
c working vectors. The user must declare this array to have dimensions
c MAXN*LNEN in the calling routine.

c LNEN->: A positive integer specifying the column dimension of the work 
c array WRKNEN in the user's calling program. It must satisfy 
c LNEN >= 2*NINT(sqrt(N))+9.

c IWRKN->: Workspace used to store the integer working vectors. The user 
c must declare this array to have dimensions at least MAXN*LIN in the
c calling routine.

c LIN->: A positive integer specifying the column dimension of the work 
c array IWRKN in the user's calling program. It must satisfy 
c LIN >= 3.

c FVEC->: The name of a user-supplied subroutine that evaluates the 
c function F at an arbitrary vector X.  The subroutine must 
c be declared EXTERNAL in the user's calling program and must conform 
c to the usage
c                      CALL FVEC(X, F, M, N),
c where X is a vector of length N and F is a vector of length M. The 
c subroutine must not alter the values of X.

c JAC->: The name of a user-supplied subroutine that evaluates the first 
c derivative (Jacobian) of the function F(x). The subroutine must be 
c declared EXTERNAL in the user's program and must conform to the usage 
c                      CALL JAC(X, AJA, MAXM, M, N)
c where X is a vector of length N and the 2-dimensional array AJA of row 
c dimension MAXM and column dimension N is the analytic Jacobian of F at 
c X. When using the interface TSNECI, if no analytic Jacobian is supplied 
c (JACFLG = 0), the user must use the dummy name TSDUMJ as the value of 
c this parameter.

c MSG<-->: An integer variable that the user may set on input to inhibit 
c certain automatic checks or to override certain default characteristics 
c of the package. (For the short call it should be set to 0.) There are 
c four "message" features that can be used individually or in combination 
c as discussed below. 
c MSG = 0 : Values of input parameters, final results, and termination code 
c are printed. 
c MSG = 2 : Do not check user's analytic Jacobian routine against its 
c finite difference estimate.  This may be necessary if the user knows the 
c Jacobian is properly coded, but the program aborts because the comparative 
c tolerance is too tight.  Do not use MSG = 2 if the analytic acobian is 
c not supplied. 
c MSG = 8 : Suppress printing of the input state, the final results, and 
c the stopping condition.
c MSG = 16 : Print the intermediate results; that is, the input state, 
c each iteration including the current iterate xk, 0.5*||DF*F(xk)||2**2, 
c and grad(f(x)) = J(x)-trans*DF**2 F(x), and the final results including 
c the stopping conditions. 
c The user may specify a combination of features by setting MSG to 
c the sum of the individual components. The module TSDFLT returns a value 
c of 0. On exit, if the program has terminated because of erroneous 
c input, MSG contains an error code indicating the reason. 
c MSG = 0   : No error. 
c MSG = -1  : Illegal dimension, M <= 0. 
c MSG = -2  : Illegal dimension, N <= 0. 
c MSG = -3  : Illegal dimension, MAXM < M+N+2.  
c MSG = -4  : Illegal dimension, MAXN < N+2.  
c MSG = -5  : Illegal dimension, MAXP < NINT(sqrt(N)).
c MSG = -6  : Illegal dimension, LUNC < 2*NINT(sqrt(N))+4. 
c MSG = -7  : Illegal dimension, LNEM < N+2*NINT(sqrt(N))+11. 
c MSG = -8  : Illegal dimension, LNEN < 2*NINT(sqrt(N))+9.
c MSG = -9  : Illegal dimension, LIN  < 3. 
c MSG = -10 : Program asked to override check of analytic Jacobian 
c against finite difference estimate, but routine JAC not 
c supplied (incompatible input).
c MSG = -11  : Probable coding error in the user's analytic Jacobian 
c routine JAC. Analytic and finite difference Jacobian do not agree 
c within the assigned tolerance. 

c XP<-: An array of length N containing the best approximation 
c to the solution x*. (If the algorithm has not converged, the final 
c iterate is returned).

c FP<-: An array of length M containing the function value F(XP).

c GP<-: An array of length N containing the gradient of the 
c function 0.5*||F(x)||2**2 at XP. 

c TERMCD<-:  An integer specifying the reason for termination.
c TERMCD = 0 : No termination criterion satisfied (occurs if package 
c terminates because of illegal input).
c TERMCD = 1 : function tolerance reached.  The current iteration is 
c probably a solution. 
c TERMCD = 2 : gradient tolerance reached.  For nonlinear least 
c squares, the current iteration is probably a solution; for nonlinear 
c equations, it could be a solution or a local minimizer.
c TERMCD = 3 : Successive iterates within step tolerance.  The 
c current iterate may be a solution, or the algorithm is making very slow
c progress and is not near a solution.
c TERMCD = 4 : Last global step failed to locate a point lower 
c than XP. It is likely that either XP is an approximate solution 
c of the problem or STEPTL is too large.
c TERMCD = 5 : Iteration limit exceeded. 
c TERMCD = 6 : Five consecutive steps of length STEPMX have been taken. 
c
c ===========================================================================
c Begin TENSOLVE
c ===========================================================================

        SUBROUTINE TS2DTR(AJA,SHAT,ANLS,DT,G,GBAR,XC,METHOD,NWTAKE,
     +                    STEPMX,STEPTL,EPSM,MXTAKE,DLT,FQ,MAXM,MAXN,
     +                    M,N,P,CURPOS,PIVOT,PBAR,ITN,IERR,FLAG,
     +                    DXN,DFN,FVEC,D,XPLSP,ADT,AG,TEMP,VN,VNP,VNS,
     +                    WRK1,CONST1,CONST2,FNORM,XPLS,FP,FPLS,RETCD)

        INTEGER MAXM,MAXN,M,N,P,ITN,METHOD,IERR,FLAG
        DOUBLE PRECISION STEPMX,STEPTL,EPSM,DLT,FPLS
        INTEGER CURPOS(N),PIVOT(N),PBAR(N),RETCD
        DOUBLE PRECISION DT(N),G(N),GBAR(N),XC(N)
        DOUBLE PRECISION XPLS(N),FP(M),XPLSP(N),AJA(MAXM,N),D(M)
        DOUBLE PRECISION TEMP(M),SHAT(MAXN,P),ANLS(MAXM,P),VNS(M)
        DOUBLE PRECISION VN(M),VNP(M),FQ(M),DXN(N),DFN(M)
        DOUBLE PRECISION ADT(N),AG(N),WRK1(M),CONST1(P),CONST2(P)
        LOGICAL NWTAKE,MXTAKE

C**********************************************************************
C THIS ROUTINE FINDS A NEXT ITERATE BY A 2-DIMENSIONAL TRUST REGION.
C**********************************************************************
C
C
C       INPUT PARAMETERS :
C       -----------------
C
C       AJA    : JACOBIAN AT THE CURRENT ITERATE
C       SHAT   : MATRIX OF PAST LINEARLY INDEPENDENT DIRECTIONS
C                AFTER A QL FACTORIZATION
C       ANLS   : TENSOR TERM MATRIX
C       DT     : CURRENT STEP
C       G      : GRADIENT AT CURRENT ITERATE
C       GBAR   : STEEPEST DESCENT DIRECTION (= -G)
C       XC     : CURRENT ITERATE
C       METHOD : METHOD TO USE
C                  =  0  : STANDARD METHOD USED
C                  =  1  : TENSOR METHOD USED
C       NWTAKE : LOGICAL VARIABLE WITH THE FOLLOWING MEANINGS:
C                NWTAKE  =  .TRUE.  : STANDARD STEP TAKEN
C                NWTAKE  =  .FALSE. : TENSOR STEP TAKEN
C       STEPMX : MAXIMUM STEP ALLOWED
C       STEPTL : STEP TOLERANCE
C       EPSM   : MACHINE PRECISION
C       MXTAKE : BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED
C       FQ     : FUNCTION VALUE AT CURRENT ITERATE MULTIPLIED BY
C                ORTHOGONAL MATRICES
C       MAXM   : LEADING DIMENSION OF AJA AND ANLS
C       MAXN   : LEADING DIMENSION OF SHAT
C       M,N    : DIMENSIONS OF PROBLEM
C       P      : COLUMN DIMENSION OF THE MATRICES SHAT AND ANLS
C       CURPOS : PIVOT VECTOR (USED DURING THE FACTORIZATION OF THE
C                JACOBIAN FROM COLUMN 1 TO N-P)
C       PIVOT  : PIVOT VECTOR ( USED DURING THE FACTORIZATION OF THE
C                JACOBIAN FROM COLUMN N-P+1 TO N)
C       PBAR   : PIVOT VECTOR (USED DURING THE FACTORIZATION OF THE
C                JACOBIAN IF IT IS SINGULAR
C       FNORM  :  0.5 * || FC ||**2
C       ITN    : ITERATION NUMBER
C       IERR   : RETURN CODE FROM THE QRP FACTORIZATION ROUTINE:
C                IERR = 0 : NO SINGULARITY OF JACOBIAN DETECTED
C                IERR = 1 : SINGULARITY OF JACOBIAN DETECTED
C       FLAG   : RETURN CODE WITH THE FOLLOWING MEANINGS :
C                FLAG  =  0 : NO SINGULARITY DETECTED DURING
C                             FACTORIZATION OF THE JACOBIAN FROM
C                             COLUMN 1 TO N
C                FLAG  =  1 : SINGULARITY DETECTED DURING FACTORIZATION
C                           OF THE JACOBIAN FROM COLUMN 1 TO N-P
C                FLAG  =  2 : SINGULARITY DETECTED DURING FACTORIZATION
C                           OF THE JACOBIAN FROM COLUMN N-P+1 TO N
C        DXN   : DIAGONAL SCALING MATRIX FOR X
C        DFN   : DIAGONAL SCALING MATRIX FOR F
C        FVEC  : SUBROUTINE TO EVALUATE THE USER'S FUNCTION
C        D,XPLSP,ADT,AG,TEMP,VN,VNP,VNS : WORKING VECTORS
C        WRK1,CONST1,CONST2,X: WORKING VECTORS
C
C       INPUT-OUTPUT PARAMETERS :
C       ------------------------
C
C       DLT    : INITIAL TRUST RADIUS (= -1.0D0) IF IT IS NOT SUPPLIED
C                BY THE USER ON ENTRY AND CURRENT TRUST RADIUS ON EXIT
C
C       OUTPUT PARAMETERS :
C       -------------------
C
C       XPLS   : NEXT ITERATE
C       FP     : FUNCTION VALUE AT NEXT ITERATE
C       FPLS   : 0.5 * || FP ||**2
C       RETCD  :  RETURN CODE FROM SUBROUTINE (SEE SUBROUTINE TSTRUD
C                 FOR MEANING )
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DAXPY,DCOPY,DDOT,DNRM2,DSCAL
C       TENSOLVE      ...  TSPRMV,TSUTMD,TSJMUV,TSUDQV,TSSMIN,TSRSID,
C       TENSOLVE      ...  TSTRUD
C
C***********************************************************************

        INTEGER I
        DOUBLE PRECISION FNORM,RES,ALPH,SUM,RESG,OPTIM
        DOUBLE PRECISION SCRES,FPLSP,RRES,RRESG
        DOUBLE PRECISION DNRM2,DDOT
        DOUBLE PRECISION ZERO,ONE
        LOGICAL DTAKEN
        INTRINSIC SQRT
        EXTERNAL FVEC
        DATA ZERO,ONE/0.0D0,1.0D0/

        DTAKEN = .FALSE.
        RETCD = 4
        IF(DLT.EQ.-ONE) THEN

c set DLT to length of Cauchy step

          ALPH = DNRM2(N,G,1)
          ALPH = ALPH**2
          CALL TSPRMV(VN,G,PIVOT,N,1)
          IF(IERR.EQ.0) THEN
             CALL TSUTMD(AJA,VN,MAXM,M,N,VNP)
          ELSE
             CALL TSPRMV(VNS,VN,PBAR,N,1)
             CALL TSUTMD(AJA,VNS,MAXM,M+N,N,VNP)
          ENDIF
          DLT = ALPH*SQRT(ALPH)/DNRM2(N,VNP,1)**2
          IF(DLT.GT.STEPMX) THEN
             DLT = STEPMX
          ENDIF
        ENDIF

c form an orthonormal basis for the two-dimensional subspace

        CALL DCOPY(N,G,1,GBAR,1)
        CALL DSCAL(N,-ONE,GBAR,1)
        RES = DNRM2(N,DT,1)
        SUM = -DDOT(N,GBAR,1,DT,1)/RES**2
        CALL DAXPY(N,SUM,DT,1,GBAR,1)
        RESG = DNRM2(N,GBAR,1)
        IF(RESG.GT.ZERO) THEN
           RRESG = ONE/RESG
           CALL DSCAL(N,RRESG,GBAR,1)
        ENDIF
        RRES = ONE/RES
        CALL DSCAL(N,RRES,DT,1)

c compute Jacobian times DT

        CALL TSJMUV(ITN,METHOD,DT,CURPOS,PIVOT,PBAR,AJA,SHAT,
     +              FLAG,IERR,MAXM,MAXN,M,N,P,D,TEMP,VN,ADT)

c compute Jacobian times GBAR

        CALL TSJMUV(ITN,METHOD,GBAR,CURPOS,PIVOT,PBAR,AJA,SHAT,
     +              FLAG,IERR,MAXM,MAXN,M,N,P,D,TEMP,VNP,AG)

        IF(.NOT. NWTAKE) THEN

c compute SHAT times VN

            CALL TSUDQV(SHAT,VN,MAXN,N,P,CONST1)

c compute SHAT times VNP

            CALL TSUDQV(SHAT,VNP,MAXN,N,P,CONST2)
        ENDIF


 70     CONTINUE

c normalize DT

        IF(RES.LE.DLT) THEN
           DTAKEN = .TRUE.
           DO 80 I = 1,N
              D(I) = DT(I)*RES
 80        CONTINUE
           DLT = RES

        ELSE

c find the global minimizer of one-variable function in the
c interval (-dlt, dlt)

           CALL TSSMIN(ANLS,FQ,ADT,AG,CONST1,CONST2,DLT,MAXM,M,N,
     +                P,NWTAKE,IERR,EPSM,VN,VNP,VNS,OPTIM)

c compute the global step

           DO 90 I = 1,N
              D(I) = OPTIM*DT(I)+SQRT(DLT**2-OPTIM**2)*GBAR(I)
 90        CONTINUE

        ENDIF

c compute the tensor model residual

        CALL TSRSID(ITN,METHOD,FQ,D,CURPOS,PIVOT,PBAR,AJA,ANLS,
     +              SHAT,FLAG,NWTAKE,IERR,MAXM,MAXN,M,N,P,WRK1,VN,
     +              VNP,VNS,SCRES)

c check whether the global step is acceptable

        CALL TSTRUD(M,N,XC,FNORM,G,D,DTAKEN,STEPMX,STEPTL,DLT,
     +              MXTAKE,DXN,DFN,FVEC,SCRES,RETCD,XPLSP,FPLSP,
     +              TEMP,XPLS,FP,FPLS)

        IF(RETCD.GE.2) GO TO 70

       END

       SUBROUTINE TSBSLV(R,NR,M,N,B,Y)

       INTEGER NR,M,N
       DOUBLE PRECISION R(NR,N),B(N),Y(N)

C*********************************************************************
C THIS ROUTINE DOES A BACKWARD SOLVE.
C*********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       R  : UPPER TRIANGULAR MATRIX OBTAINED FROM A QR FACTORIZATION
C            OF AN M BY N MATRIX A. DIAG(R) IS STORED IN ROW M+2. THIS
C            IS THE STORAGE SCHEME USED IN STEWART, G. W., III(1973)
C            "INTRODUCTION TO MATRIX COMPUTATION", ACADEMIC PRESS,
C             NEW YORK
C       NR : LEADING DIMENSION OF MATRIX A
C       M  : ROW DIMENSION OF MATRIX A
C       N  : COLUMN DIMENSION OF MATRIX A
C       B  : RIGHT HAND SIDE
C
C       OUTPUT PARAMETERS :
C       -----------------
C
C       Y :  VECTOR SOLUTION ON EXIT
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DAXPY
C       TENSOLVE      ...  TSDLOD
C
C*********************************************************************

       INTEGER J
       DOUBLE PRECISION ZERO
       DATA ZERO/0.0D0/

c solve R Y = B

       Y(N) = B(N) / R(M+2,N)

       IF(N .GT. 2) THEN
          CALL TSDLOD(N-1,ZERO,Y,1)
          DO 20 J = N-1,2,-1
             CALL DAXPY(J,Y(J+1),R(1,J+1),1,Y,1)
             Y(J) = (B(J)-Y(J))/R(M+2,J)
 20       CONTINUE
          Y(1) = Y(1) + R(1,2) * Y(2)
          Y(1) = (B(1) - Y(1)) / R(M+2,1)
       ELSEIF(N .EQ. 2) THEN
          Y(1) = (B(1) - (R(1,2) * Y(2))) / R(M+2,1)
       ENDIF

       RETURN
       END

       SUBROUTINE TSCHKI(MAXM,MAXN,MAXP,M,N,LUNC,LNEM,LNEN,LIN,GRADTL,
     +                   STEPTL,FTOL,ITNLIM,JACFLG,METHOD,GLOBAL,
     +                   STEPMX,DLT,EPSM,MSG,TYPX,TYPF,DXN,DFN,
     +                   SQRN,TERMCD,IPR)

       INTEGER MAXM,MAXN,MAXP,M,N,LUNC,LNEM,LNEN,LIN,IPR,MSG,JACFLG
       INTEGER METHOD,GLOBAL,ITNLIM,SQRN,TERMCD
       DOUBLE PRECISION GRADTL,STEPTL,FTOL,STEPMX,DLT,EPSM
       DOUBLE PRECISION TYPX(N),TYPF(M),DXN(N),DFN(M)

C*********************************************************************
C THIS ROUTINE CHECKS INPUT FOR REASONABLENESS.
C*********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       MAXM : LEADING DIMENSION OF WORKSPACE WRKNEM
C              (SEE TOP OF THIS FILE )
C       MAXN : LEADING DIMENSION OF WORKSPACE WRKNEN
C              (SEE TOP OF THIS FILE )
C       MAXP : LEADING DIMENSION OF WORKSPACE WRKUNC
C              (SEE TOP OF THIS FILE )
C       M,N  : DIMENSIONS OF PROBLEM
C       IPR  : DEVICE TO WHICH TO SEND OUTPUT
C
C       INPUT-OUTPUT PARAMETERS :
C       ------------------------
C
C       GRADTL : TOLERANCE AT WHICH GRADIENT CONSIDERED CLOSE
C                ENOUGH TO ZERO TO TERMINATE ALGORITHM
C       STEPTL : TOLERANCE AT WHICH SUCCESSIVE ITERATES
C                CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM
C       FTOL   : TOLERANCE AT WHICH FUNCTION VALUE CONSIDERED
C                CLOSE ENOUGH TO ZERO
C       ITNLIM : MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
C       STEPMX : MAXIMUM STEP ALLOWED IN TRUST REGION
C       DLT    : TRUST RADIUS
C       JACFLG : JACOBIAN FLAG WITH THE FOLLOWING MEANINGS :
C                JACFLG = 1 : ANALYTIC JACOBIAN SUPPLIED
C                JACFLG = 0 : ANALYTIC JACOBIAN NOT SUPPLIED
C       METHOD : METHOD TO USE
C                METHOD = 0 : STANDARD METHOD IS USED
C                METHOD = 1 : TENSOR METHOD IS USED
C       GLOBAL : GLOBAL STRATEGY USED
C                GLOBAL = 0 : LINE SEARCH USED
C                GLOBAL = 1 : 2-DIMENSIONAL TRUST REGION USED
C       TYPX   : TYPICAL SIZE FOR EACH COMPONENT OF X
C       TYPF   : TYPICAL SIZE FOR EACH COMPONENT OF F
C       MSG    : MESSAGE TO INHIBIT CERTAIN AUTOMATIC CHECKS + OUTPUT
C
C       OUTPUT PARAMETERS :
C       -------------------
C
C       TERMCD: TERMINATION CODE
C       DXN   : DIAGONAL SCALING MATRIX FOR X
C       DFN   : DIAGONAL SCALING MATRIX FOR F
C       SQRN  : MAXIMUM COLUMN DIMENSION OF S AND FV
C
C       SUBPROGRAMS CALLED:
C
C       UNCMIN        ...  DPMEPS
C
C*********************************************************************

       INTEGER I,LEN
       DOUBLE PRECISION DPMEPS,ZERO,ONE,TWO,THREE,THOUS,TEMP
       INTRINSIC MOD,NINT,SQRT
       DATA ZERO,ONE,TWO,THREE,THOUS/0.0D0,1.0D0,2.0D0,3.0D0,1000.0D0/

c check that parameters only take on acceptable values
c if not, set them to default values

c set TERMCD to zero in case we abort prematuraly

       TERMCD = 0

c compute machine precision

       EPSM = DPMEPS()

c check dimensions of the problem

       IF(M.LE.0) THEN
          WRITE(IPR,601) M
          MSG = -1
          RETURN
       ENDIF

       IF(N.LE.0) THEN
          WRITE(IPR,602) N
          MSG = -2
          RETURN
       ENDIF

c check leading dimensions of the problem

        LEN = M+N+2
        IF(MAXM .LT. LEN) THEN
           WRITE(IPR,603) MAXM,LEN
           MSG = -3
           RETURN
        ENDIF

        LEN = N+2
        IF(MAXN .LT. LEN) THEN
           WRITE(IPR,604) MAXN,LEN
           MSG = -4
           RETURN
        ENDIF

        TEMP = SQRT(DBLE(N))
        SQRN = NINT(TEMP)

        IF(MAXP .LT. SQRN) THEN
           WRITE(IPR,605) MAXP,SQRN
           MSG = -5
           RETURN
        ENDIF

c check column dimensions of workspace arrays

       LEN = 2*SQRN+4
       IF(LUNC.LT.LEN) THEN
          WRITE(IPR,606) LUNC,LEN
          MSG = -6
          RETURN
       ENDIF

       LEN = N+2*SQRN+11
       IF(LNEM.LT.LEN) THEN
          WRITE(IPR,607) LNEM,LEN
          MSG = -7
          RETURN
       ENDIF

       LEN = 2*SQRN+9
       IF(LNEN.LT.LEN) THEN
          WRITE(IPR,608) LNEN,LEN
          MSG = -8
          RETURN
       ENDIF

       IF(LIN.LT.3) THEN
          WRITE(IPR,609) LIN
          MSG = -9
          RETURN
       ENDIF

c check JACFLG, METHOD, and GLOBAL

       IF(JACFLG.NE.1) JACFLG = 0

       IF(METHOD.NE.0 .AND. METHOD.NE.1) METHOD = 1

       IF(GLOBAL.NE.0 .AND. GLOBAL.NE.1) GLOBAL = 0

       IF(MOD(MSG/2,2).EQ.1 .AND. JACFLG.EQ.0) THEN
          WRITE(IPR,610) MSG,JACFLG
          MSG = -10
          RETURN
       ENDIF

c check scale matrices

       DO 10 I = 1,N
          IF(TYPX(I).EQ.ZERO) TYPX(I) = ONE
          IF(TYPX(I).LT.ZERO) TYPX(I) = -TYPX(I)
          DXN(I) = ONE/TYPX(I)
 10    CONTINUE

       DO 20 I = 1,M
          IF(TYPF(I).EQ.ZERO) TYPF(I) = ONE
          IF(TYPF(I).LT.ZERO) TYPF(I) = -TYPF(I)
          DFN(I) = ONE/TYPF(I)
 20    CONTINUE

c check gradient, step, and function tolerances

       TEMP = ONE/THREE
       IF(GRADTL.LT.ZERO) THEN
          GRADTL = EPSM**TEMP
       ENDIF

       IF(STEPTL.LT.ZERO) THEN
          STEPTL = EPSM**(TWO*TEMP)
       ENDIF

       IF(FTOL.LT.ZERO) THEN
          FTOL = EPSM**(TWO*TEMP)
       ENDIF

c check iteration limit

       IF(ITNLIM.LE.0) THEN
          ITNLIM = 150
       ENDIF

c check STEPMX and DLT

       IF(STEPMX.LT.ZERO) STEPMX = THOUS

       IF(DLT.LE.ZERO) THEN
          DLT = -ONE
          IF(DLT.GT.STEPMX) DLT = STEPMX
       ENDIF

 601  FORMAT('  TSCHKI     ILLEGAL DIMENSION M =',I5)
 602  FORMAT('  TSCHKI     ILLEGAL DIMENSION N =',I5)
 603  FORMAT('  TSCHKI     ILLEGAL DIMENSION MAXM =',I5,' < M+N+2 =',I5)
 604  FORMAT('  TSCHKI     ILLEGAL DIMENSION MAXN =',I5,' < N+2 =',I5)
 605  FORMAT('  TSCHKI     ILLEGAL DIMENSION MAXP =',I5,' <',
     +       '  NINT(SQRT (N)) =',I5)
 606  FORMAT('  TSCHKI     ILLEGAL DIMENSION LUNC =',I5,' <',
     +       '  2*NINT(SQRT (N))+4 =',I5)
 607  FORMAT('  TSCHKI     ILLEGAL DIMENSION LNEM =',I5,' <',
     +       '  N+2*NINT(SQRT (N))+11 =',I5)
 608  FORMAT('  TSCHKI     ILLEGAL DIMENSION LNEN =',I5,' <',
     +       '  2*NINT(SQRT (N))+9 =',I5)
 609  FORMAT('  TSCHKI     ILLEGAL DIMENSION LIN =',I5,' < 3')
 610  FORMAT('  TSCHKI     USER REQUESTS THAT ANALYTIC JACOBIAN BE',
     +       ' ACCEPTED AS PROPERLY CODED (MSG =',I5,')'/
     +       '  TSCHKI     BUT ANALYTIC JACOBIAN NOT SUPPLIED',
     +       ' (JACFLG =',I5,')')
      END

       SUBROUTINE TSCHKJ(AJANAL,XC,FC,NR,M,N,EPSM,DFN,DXN,
     +                   TYPX,IPR,FHAT,WRK1,FVEC,MSG)

       INTEGER NR,M,N,IPR,MSG
       DOUBLE PRECISION AJANAL(NR,N),XC(N),FC(M)
       DOUBLE PRECISION EPSM,DFN(M),DXN(N),TYPX(N)
       DOUBLE PRECISION FHAT(M),WRK1(M)
       EXTERNAL FVEC

C*********************************************************************
C THIS ROUTINE CHECKS THE ANALYTIC JACOBIAN AGAINST ITS FINITE
C DIFFERENCE APPROXIMATION.
C*********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       AJANAL : ANALYTIC JACOBIAN AT XC
C       XC   : CURRENT ITERATE
C       FC   : FUNCTION VALUE AT XC
C       NR   : LEADING DIMENSION OF AJANAL
C       M,N  : DIMENSIONS OF PROBLEM
C       EPSM : MACHINE PRECISION
C       DFN  : DIAGONAL SCALING MATRIX FOR F
C       DXN  : DIAGONAL SCALING MATRIX FOR X
C       TYPX : TYPICAL SIZE FOR EACH COMPONENT OF X
C       IPR  : DEVICE TO WHICH TO SEND OUTPUT
C       FHAT,WRK1 : WORKSPACE
C       FVEC  : SUBROUTINE TO EVALUATE THE USER'S FUNCTION
C
C       INPUT-OUTPUT PARAMETERS :
C       ------------------------
C
C       MSG : MESSAGE TO INHIBIT CERTAIN AUTOMATIC CHECKS + OUTPUT
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  IDAMAX
C       TENSOLVE      ...  TSUNSX,TSUNSF,TSSCLX,TSSCLF
C       USER          ...  FVEC
C
C*********************************************************************

       INTEGER I,J
       DOUBLE PRECISION NDIGIT,RNOISE,SQRNS,STEPSZ,XTMPJ,DINF,RSTPSZ
       DOUBLE PRECISION TOL,QUART,ONE,TEN
       INTEGER IDAMAX
       INTRINSIC ABS,MAX,SQRT
       DATA QUART,ONE,TEN/0.250D0,1.0D0,10.0D0/

c unscale XC and FC

       CALL TSUNSX(XC,DXN,N)
       CALL TSUNSF(FC,DFN,M)

c compute the finite difference Jacobian and check it against
c the analytic one

       NDIGIT = -LOG10(EPSM)
       RNOISE = MAX(TEN**(-NDIGIT),EPSM)
       SQRNS  = SQRT(RNOISE)
       TOL = EPSM**QUART

       DO 40 J = 1,N
          STEPSZ = SQRNS*MAX(ABS(XC(J)),ONE)
          XTMPJ = XC(J)
          XC(J) = XTMPJ+STEPSZ
          CALL FVEC(XC,FHAT,M,N)
          XC(J) = XTMPJ
          RSTPSZ = ONE/STEPSZ
          DO 10 I = 1,M
             WRK1(I) = (FHAT(I)-FC(I))*RSTPSZ
 10       CONTINUE
          DO 20 I = 1,M
             WRK1(I) = WRK1(I)*DFN(I)*TYPX(J)
 20       CONTINUE
          DINF = ABS(WRK1(IDAMAX(M,WRK1,1)))
          DO 30 I = 1,M
             IF(ABS(AJANAL(I,J)-WRK1(I)).GT.TOL*DINF) THEN
                WRITE(IPR,50)
                MSG = -11
                RETURN
             ENDIF
 30       CONTINUE
 40    CONTINUE

c scale back XC and FC

       CALL TSSCLX(XC,DXN,N)
       CALL TSSCLF(FC,DFN,M)

 50    FORMAT(/,'  TSCHKJ      PROBABLE ERROR IN CODING OF ANALYTIC',
     +        ' JACOBIAN')

       RETURN
       END

       SUBROUTINE TSCPMU(R,NR,N,EPSM,MU)

       INTEGER NR,N
       DOUBLE PRECISION R(NR,N),EPSM,MU

C*********************************************************************
C THIS ROUTINE COMPUTES A SMALL PERTURBATION MU. MU IS USED IN THE
C COMPUTATION OF THE LEVENBERG-MARQUARDT STEP.
C*********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       R  : UPPER TRIANGULAR MATRIX
C       NR : LEADING DIMENSION OF R
C       N  : COLUMN DIMENSION OF R
C       EPSM :  MACHINE PRECISION
C
C       OUTPUT PARAMETERS :
C       ------------------
C
C       MU = SQRT(L1 NORM OF R * INFINITY NORM OF R * N * EPSM * 100)
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DASUM
C
C*********************************************************************

       INTEGER I,J
       DOUBLE PRECISION AIFNRM,SUM,AL1NRM,ZERO,HUND
       DOUBLE PRECISION DASUM
       INTRINSIC ABS,MAX,SQRT
       DATA ZERO,HUND/0.0D0,100.0D0/

c compute the infinity norm of R

       AIFNRM = ZERO
       DO 20 I = 1,N
          SUM = ZERO
          DO 10 J = I,N
             SUM = SUM+ABS(R(I,J))
  10      CONTINUE
          AIFNRM = MAX(AIFNRM,SUM)
  20   CONTINUE

c compute the l1 norm of R

       AL1NRM = ZERO
       DO 40 J = 1,N
          SUM = DASUM(J,R(1,J),1)
          AL1NRM = MAX(AL1NRM,SUM)
 40    CONTINUE

c compute MU

       MU = SQRT(AIFNRM*AL1NRM*N*EPSM*HUND)

       RETURN
       END

       SUBROUTINE TSCPSS(S,MAXM,MAXN,M,N,P,METHOD,GLOBAL,EPSM,FCQ,
     +                   Y,W,FQT,AL2NRM,QHAT,ANLS,DN,FQQ,PTILDA,
     +                   CURPOS,PBAR,ZERO1,IERR,RESNEW,FLAG)

       INTEGER MAXM,MAXN,M,N,P,FLAG,ZERO1,GLOBAL,IERR
       DOUBLE PRECISION EPSM,RESNEW
       INTEGER METHOD,PTILDA(N),CURPOS(N),PBAR(N)
       DOUBLE PRECISION S(MAXN,P),FCQ(M)
       DOUBLE PRECISION Y(N),W(M),FQT(M),AL2NRM(N)
       DOUBLE PRECISION QHAT(MAXM,N),ANLS(MAXM,P)
       DOUBLE PRECISION DN(N),FQQ(M)

C**********************************************************************
C THIS ROUTINE COMPUTES THE STANDARD STEP.  NOTE THAT AT THIS STAGE
C THE JACOBIAN MATRIX (QHAT) HAS ALREADY BEEN FACTORED FROM COLUMNS 1
C TO N-P DURING THE TENSOR STEP COMPUTATION.  THIS ROUTINE FACTORS
C THE MATRIX QHAT FROM COLUMN N-P+1 TO N, THEREBY OBTAINING A QR
C FACTORIZATION OF THE FULL MATRIX QHAT, THEN COMPUTES THE STANDARD
C STEP BY PREMULTIPLYING THE RIGH-HAND SIDE FCQ BY AN ORTHOGONAL
C MATRIX AND BY PERFORMING A BACKWARD SOLVE.
C**********************************************************************
C
C
C       INPUT PARAMETERS :
C       -----------------
C
C       S    : FACTORED MATRIX OF PAST LINEARLY INDEPENDENT DIRECTIONS
C             (OBTAINED FROM TSQLFC SUBROUTINE)
C       MAXM : LEADING DIMENSION OF QHAT AND ANLS
C       MAXN : LEADING DIMENSION OF S
C       M,N  : DIMENSIONS OF PROBLEM
C       P    : COLUMN DIMENSION OF MATRIX S
C       METHOD : METHOD USED :
C                METHOD = 0 : STANDARD METHOD IS USED
C                METHOD = 1 : TENSOR METHOD IS USED
C       GLOBAL : GLOBAL STRATEGY USED
C                GLOBAL = 0 : LINE SEARCH IS USED
C                GLOBAL = 1 : 2-DIMENSIONAL TRUST REGION IS USED
C       EPSM   : MACHINE PRECISION
C       FCQ    : FUNCTION VALUE AT CURRENT ITERATE MULTIPLIED BY AN
C                ORTHOGONAL MATRIX
C       CURPOS : PIVOT VECTOR FOR THE FACTORIZATION OF QHAT FROM COLUMN
C                1 TO N-P
C       Y,W,FQT,AL2NRM : WORKING VECTORS
C
C
C       INPUT-OUTPUT PARAMETERS :
C       ------------------------
C
C       QHAT  : FACTORED MATRIX FROM COLUMN 1 TO N-P
C               ON ENTRY AND FACTORED MATRIX FROM 1 TO N ON EXIT
C       ANLS  : TENSOR TERM MATRIX ON ENTRY AND ANLS MULTIPLIED BY
C               ORTHOGONAL MATRICES ON EXIT (THIS IS PERFORMED IN THE
C               CASE WHERE THE GLOBAL STRATEGY USED IS THE 2-DIMENSIONAL
C               TRUST REGION)
C
C       OUTPUT PARAMETERS :
C       -------------------
C
C       DN    : STANDARD STEP
C       FQQ   : FUNCTION VALUE AT CURRENT ITERATE MULTIPLIED BY
C               ORTHOGONAL MATRICES (THIS IS USED IN THE CASE WHERE
C               THE GLOBAL STRATEGY USED IS THE 2-DIMENSIONAL
C               TRUST REGION)
C       PTILDA: PIVOT VECTOR FOR THE FACTORIZATION OF THE
C               MATRIX QHAT FROM N-P+1 TO N
C       PBAR  : PIVOT VECTOR FOR THE FACTORIZATION OF THE
C               TRANSFORMED MATRIX QHAT FROM 1 TO N
C               IN CASE OF SINGULARITY
C       ZERO1 : FIRST ZERO COLUMN OF MATRIX QHAT IN CASE OF SINGULARITY
C       IERR  : RETURNED CODE WITH THE FOLLOWING MEANING :
C               IERR = 1 : SINGULARITY OF JACOBIAN DETECTED
C               IERR = 0 : OTHERWISE
C       RESNEW: RESIDUAL OF THE STANDARD MODEL
C       FLAG  : RETURNED CODE WITH THE FOLLOWING MEANINGS :
C               FLAG = 0 : NO SINGULARITY DETECTED
C               FLAG = 1 : SINGULARITY DETECTED DURING QR FACTORIZATION
C                          OF QHAT FROM COLUMN 1 TO N-P
C               FLAG = 2 : SINGULARITY DETECTED DURING QR FACTORIZATION
C                          OF QHAT FROM COLUMN N-P+1 TO N
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DCOPY,DSCAL
C       TENSOLVE      ...  TSQRFC,TSQMUV,TSQMTS,TSSMRD,TSBSLV,TSPRMV
C       TENSOLVE      ...  TSDLOD,TSQMLV,TSCPMU
C
C **********************************************************************

       INTEGER ZEROTM,I,J
       DOUBLE PRECISION MU,ZERO,ONE
       DATA ZERO,ONE/0.0D0,1.0D0/

       FLAG = 0

c initialization

       CALL TSDLOD (M+N,ZERO,FQQ,1)

       CALL DCOPY(M,FCQ,1,W,1)
       CALL DSCAL(M,-ONE,W,1)

c if the Jacobian is singular then compute the Levenberg-Marquardt
c step (label 20)

       IF(IERR.EQ.1) THEN
          FLAG = 1
          GO TO 20
       ENDIF

c factor the matrix QHAT from column n-p+1 to n

       CALL TSQRFC(QHAT,MAXM,N,M,N-P+1,N,IERR,EPSM,AL2NRM,PTILDA,ZERO1)

       IF((M.EQ.N).AND.(IERR.EQ.0)) THEN
           ZEROTM = ZERO1-1
       ELSE
           ZEROTM = ZERO1
       ENDIF

c premultiply W by the orthogonal matrix resulting from the QR
c factorization of QHAT

       CALL TSQMUV(QHAT,W,FQQ,MAXM,M,N-P+1,ZEROTM,.FALSE.)

       IF(METHOD.EQ.1 .AND. GLOBAL.EQ.1) THEN

c premultiply ANLS by the orthogonal matrix resulting from the QR
c factorization of QHAT

          CALL TSQMTS(ANLS,QHAT,MAXM,M,N,M,P,N-P+1,FCQ,ZEROTM)
       ENDIF

       IF(IERR.EQ.1) THEN
          FLAG = 2
          GO TO 20
       ENDIF

c computate the residual of the standard model

       CALL TSSMRD(FQQ,RESNEW,DN,MU,IERR,M,N)

c if QHAT is nonsingular perform a backward solve to obtain Y

       CALL TSBSLV(QHAT,MAXM,M,N,FQQ,Y)

c pivot Y

       CALL TSPRMV(DN,Y,PTILDA,N,0)

       IF(N .NE. 1) THEN

          CALL TSPRMV(Y,DN,CURPOS,N,0)

c premultiply Y by the orthogonal matrix resulting from the QL
c factorization of S

          CALL TSQMLV(MAXN,N,P,S,Y,DN,.TRUE.)

       ENDIF

       IF(GLOBAL.EQ.1) THEN
          IERR = 0
          CALL DSCAL(M,-ONE,FQQ,1)
       ENDIF

       RETURN

 20    CONTINUE

c                    @   SINGULAR CASE   @

c solve ( QHAT-trans QHAT + MU I ) DN = -QHAT-trans W

c put the diagonal elements stored in row m+2 of QHAT into their
c propre positions and zero out the unwanted portions of QHAT

       DO 30 J = 1, ZERO1-1
          QHAT(J,J) = QHAT(M+2,J)
          CALL TSDLOD (M+N-J,ZERO,QHAT(J+1,J),1)
 30    CONTINUE

       DO 40 J = ZERO1, N
          CALL TSDLOD (M+N-ZERO1+1,ZERO,QHAT(ZERO1,J),1)
 40    CONTINUE

c compute a small perturbation MU

       CALL TSCPMU(QHAT,MAXM,N,EPSM,MU)

c form the augmented matrix QHAT by adding an (n*n) diag(MU) in
c the bottom

       DO 50 I = M+1,M+N
          QHAT(I,I-M) = MU
 50    CONTINUE

c factor the transformed matrix QHAT from 1 to n

       CALL TSQRFC(QHAT,MAXM,N,M+N,1,N,IERR,EPSM,AL2NRM,PBAR,ZERO1)

       IF(METHOD.EQ.1 .AND. GLOBAL.EQ.1) THEN

c premultiply ANLS by the orthogonal matrix resulting from the QR
c factorization of QHAT

         CALL TSQMTS(ANLS,QHAT,MAXM,M+N,N,M,P,1,FCQ,ZERO1)
       ENDIF

c compute the Levenberg-Marquardt step and the residual of the
c standard model

       IF(FLAG.EQ.1) THEN
          CALL TSQMUV(QHAT,W,FQQ,MAXM,M+N,1,N+1,.FALSE.)
          CALL TSBSLV(QHAT,MAXM,M+N,N,FQQ,Y)
          CALL TSPRMV(DN,Y,PBAR,N,0)
          CALL TSPRMV(Y,DN,CURPOS,N,0)
          CALL TSQMLV(MAXN,N,P,S,Y,DN,.TRUE.)
          CALL TSSMRD(FQQ,RESNEW,DN,MU,IERR,M,N)
          IF(GLOBAL.EQ.1) THEN
             IERR = 1
             CALL DSCAL(M+N,-ONE,FQQ,1)
          ENDIF
          RETURN
       ELSE
          CALL TSQMUV(QHAT,FQQ,FQT,MAXM,M+N,1,N+1,.FALSE.)
          CALL TSBSLV(QHAT,MAXM,M+N,N,FQT,DN)
          CALL TSPRMV(Y,DN,PBAR,N,0)
          CALL TSPRMV(DN,Y,PTILDA,N,0)
          CALL TSPRMV(Y,DN,CURPOS,N,0)
          CALL TSQMLV(MAXN,N,P,S,Y,DN,.TRUE.)
          CALL TSSMRD(FQT,RESNEW,DN,MU,IERR,M,N)
          IF(GLOBAL.EQ.1) THEN
             IERR = 1
             CALL DCOPY(M+N,FQT,1,FQQ,1)
             CALL DSCAL(M+N,-ONE,FQQ,1)
          ENDIF
       ENDIF

       END


      SUBROUTINE TSDLOD ( N, CONST, X, INCX )

      DOUBLE PRECISION   CONST
      INTEGER            INCX, N
      DOUBLE PRECISION   X(*)

C**********************************************************************
C THIS ROUTINE LOADS ELEMENTS OF X WITH CONST.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       ----------------
C
C      N     :  DIMENSION OF THE VECTOR X
C      CONST :  CONSTANT VALUE
C      INCX  :  INCREMENT
C
C      OUTPUT PARAMETERS :
C      ------------------
C
C      X     :  VECTOR WITH ELEMENTS EQUAL TO CONST
C
C**********************************************************************

      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
      INTEGER            IX

      IF (N .GT. 0) THEN
         IF (INCX .EQ. 1  .AND.  CONST .EQ. ZERO) THEN
            DO 10 IX = 1, N
               X(IX) = ZERO
   10       CONTINUE
         ELSE
            DO 20 IX = 1, 1 + (N - 1)*INCX, INCX
               X(IX) = CONST
   20       CONTINUE
         ENDIF
      ENDIF
      END


       SUBROUTINE TSD1SV(AJA,S,ANLS,FN,X,MAXM,MAXN,M,N,P,Q,EPSM,
     +                   WRK1,WRK2,WRK3,PIVOT,D1)

       INTEGER MAXM,MAXN,M,N,P,Q
       INTEGER PIVOT(N)
       DOUBLE PRECISION EPSM
       DOUBLE PRECISION AJA(MAXM,N),S(MAXN,P),ANLS(MAXM,P),FN(M),X(P)
       DOUBLE PRECISION WRK1(N),WRK2(N),WRK3(N),D1(N)

C*********************************************************************
C THIS ROUTINE SOLVES THE FIRST N-Q LINEAR EQUATIONS IN N-P UNKNOWNS
C OF THE TENSOR MODEL.
C*********************************************************************
C
C       INPUT PARAMETERS :
C       ----------------
C
C       AJA : JACOBIAN MATRIX AT CURRENT ITERATE
C       S   : MATRIX OF PAST LINEARLY INDEPENDENT DIRECTIONS
C       ANLS: TENSOR TERM MATRIX AT CURRENT ITERATE
C       FN  : FUNCTION VALUE AT CURRENT ITERATE
C        X  : SOLUTION OF THE LOWER M-N+Q QUADRATIC EQUATIONS IN P
C             UNKNOWNS OF THE TENSOR MODEL
C       MAXM: LEADING DIMENSION OF AJA AND ANLS
C       MAXN: LEADING DIMENSION OF S
C       M,N : DIMENSIONS OF PROBLEM
C       P   : COLUMN DIMENSION OF S AND ANLS
C       Q   : NUMERICAL RANK OF JACOBIAN :
C             Q > P : JACOBIAN IS SINGULAR
C             Q = P : OTHERWISE
C       EPSM: MACHINE PRECISION
C       WRK1,WRK2,WRK3 : WORKSPACE
C
C
C       OUTPUT PARAMETERS :
C       ------------------
C
C        PIVOT : PIVOT VECTOR
C        D1 : SOLUTION OF THE N-Q LINEAR EQUATIONS IN N-P UNKNOWNS OF
C             THE TENSOR MODEL
C 
C       SUBPROGRAMS CALLED:
C 
C       LEVEL 1 BLAS  ...  DCOPY
C       LEVEL 2 BLAS  ...  DGEMV
C       TENSOLVE      ...  TSDLOD,TSSTMX,TSBSLV,TSQRFC,TSPRMV
C       TENSOLVE      ...  TSFSLV,TSQMUV
C
C*********************************************************************

       INTEGER ZERO1,I,J,IERR,ICOL
       DOUBLE PRECISION EPSM1,ZERO,HALF,ALPHA,ONE
       DATA ZERO,ALPHA,HALF,ONE/0.0D0,1.0D-4,0.50D0,1.0D0/

c compute the top right (n-q) x p submatrix of AJA times X

       CALL DGEMV('N',N-Q,P,ONE,AJA(1,N-P+1),MAXM,
     +            X,1,ZERO,D1,1)

c compute S-trans times X

       CALL TSSTMX(S,X,MAXN,N,P,WRK3,WRK2)

c compute 0.5 * (S-trans times X)**2

       DO 10 I = 1, P
          WRK1(I) = HALF * WRK2(I)**2
  10   CONTINUE

c compute 0.5 * (top (n-q) x p submatrix of ANLS) *
c (S-trans times X)**2

       CALL DGEMV('N',N-Q,P,ONE,ANLS(1,1),MAXM,WRK1,1,ZERO,WRK2,1)

       DO 20 I = 1,N-Q
          WRK1(I) = -FN(I)-D1(I)-WRK2(I)
 20    CONTINUE

c if the Jacobian is nonsingular then solve for the first
c n-p components of the tensor step and return

       IF(P.EQ.Q) THEN
          CALL TSBSLV(AJA,MAXM,M,N-P,WRK1,D1)
          RETURN
       ENDIF

       CALL TSDLOD(Q-P,ZERO,WRK2(N-Q+1),1)

c copy top left (n-q) x (n-p) submatrix of AJA into bottom of AJA

       DO 30 J = 1,N-P
          CALL DCOPY(N-Q,AJA(1,J),1,AJA(M+3,J),1)
 30    CONTINUE

c copy the transpose of the top left (n-q) x (n-p) submatrix of AJA
c into top of AJA

       DO 50 J = 1,N-Q
          AJA(J,J) = AJA(M+2,J)
          DO 40 I = J+1,N-P
            AJA(I,J) = AJA(J,I)
 40       CONTINUE
 50    CONTINUE

c zero out the upper triangular (n-q) x (n-q) triangular part of
c the transpose of the top left (n-q) x (n-p) submatrix of AJA

       DO 60 J = 1,N-Q
          CALL TSDLOD(J-1,ZERO,AJA(1,J),1)
 60    CONTINUE

c factor the transpose of the top left (n-q) x (n-p) submatrix of AJA

       EPSM1 = EPSM*ALPHA

       CALL TSQRFC(AJA,MAXM,N-Q,N-P,1,N-Q,IERR,EPSM1,WRK3,PIVOT,ZERO1)

       IF(IERR .EQ. 0) THEN
          ICOL = N-Q
       ELSE
          ICOL = ZERO1-1
       ENDIF

       CALL TSPRMV(D1,WRK1,PIVOT,N-Q,0)

c solve for the first n-p components of the tensor step

       CALL TSFSLV(AJA,D1,MAXM,N-P,ICOL,WRK2)

       CALL TSQMUV(AJA,WRK2,D1,MAXM,N-P,1,ZERO1,.TRUE.)

c copy the (n-q) x (n-p) submatrix of AJA from bottom back to
c top of AJA

       DO 70 J = 1,N-P
          CALL DCOPY(N-Q,AJA(M+3,J),1,AJA(1,J),1)
 70    CONTINUE

       RETURN
       END

       SUBROUTINE TSDFCN(P,X,G,AJA,ANLS,S,FN,WRK1,WRK2,
     +                   WRK3,WRK4,WRK5,MAXM,MAXN,M,N,Q)

       INTEGER P,MAXM,MAXN,M,N,Q
       DOUBLE PRECISION X(P),G(P),AJA(MAXM,N),ANLS(MAXM,P),S(MAXN,P)
       DOUBLE PRECISION FN(M),WRK1(M),WRK2(P),WRK3(P),WRK4(M),WRK5(M)

C*********************************************************************
C THIS ROUTINE COMPUTES THE ANALYTIC GRADIENT OF THE FUNCTION GIVEN
C BY SUBROUTINE TSQFCN.
C*********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       P  : COLUMN DIMENSION OF ANLS AND S
C       X  : POINT AT WHICH GRADIENT IS EVALUATED
C       AJA: JACOBIAN MATRIX AT CURRENT ITERATE
C       ANLS : TENSOR TERM MATRIX AT CURRENT ITERATE
C       S  : MATRIX OF PAST LINEARLY INDEPENDENT DIRECTIONS
C       FN : FUNCTION VALUE AT CURRENT ITERATE
C       WRK1,WRK2,WRK3,WRK4,WRK5 : WORKSPACE
C       MAXM : LEADING DIMENSION OF AJA AND ANLS
C       MAXN : LEADING DIMENSION OF S
C       M,N  : DIMENSIONS OF PROBLEM
C       Q : NUMERICAL RANK OF JACOBIAN :
C           Q > P : JACOBIAN IS SINGULAR
C           Q = P : OTHERWISE
C
C
C       OUTPUT PARAMETERS :
C       -----------------
C
C       G : GRADIENT AT X
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DAXPY,DDOT
C       LEVEL 2 BLAS  ...  DGEMV
C       TENSOLVE      ...  TSSTMX,TSDLOD
C
C*********************************************************************

       INTEGER I,J,K,L
       DOUBLE PRECISION ZERO,HALF,ONE
       DOUBLE PRECISION DDOT
       DATA ZERO,HALF,ONE/0.0D0,0.50D0,1.0D0/

c compute the lower right (m-n+q) x p submatrix of AJA times X

       CALL DGEMV('N',M-N+Q,P,ONE,AJA(N-Q+1,N-P+1),MAXM,
     +            X,1,ZERO,WRK1,1)

c compute S-trans times X

       CALL TSSTMX(S,X,MAXN,N,P,WRK2,WRK3)

c compute 0.5 * (S-trans times X)**2

       DO 10 I = 1, P
          WRK2(I) = HALF * WRK3(I)**2
 10    CONTINUE

c compute 0.5 * (lower (m-n+q) x p submatrix of ANLS) *
c (S-trans times X)**2

       CALL DGEMV('N',M-N+Q,P,ONE,ANLS(N-Q+1,1),MAXM,
     +            WRK2,1,ZERO,WRK4,1)

       DO 20 I = 1,M-N+Q
          WRK4(I) = WRK4(I)+FN(N-Q+I)+WRK1(I)
 20    CONTINUE

c compute AJA-trans * WRK4

       CALL DGEMV('T',M-N+Q,P,ONE,AJA(N-Q+1,N-P+1),MAXM,
     +            WRK4,1,ZERO,WRK1,1)

c compute ANLS-trans * WRK4

       CALL DGEMV('T',M-N+Q,P,ONE,ANLS(N-Q+1,1),MAXM,
     +            WRK4,1,ZERO,WRK5,1)

c compute S * diag(S-trans * WRK3) * WRK5

       CALL TSDLOD(P,ZERO,WRK2,1)

       L = P+1
       DO 50 J = 1,P
          L = L-1
          WRK2(L) = S(N+2,L)
          DO 30 I = L+1,P
             WRK2(I) = S(N-P+J,I)
 30       CONTINUE
          DO 40 K = 1,P
             WRK2(K) = WRK2(K)*WRK3(K)
 40       CONTINUE
          G(J) = DDOT(P,WRK2,1,WRK5,1)
 50    CONTINUE

       CALL DAXPY(P,ONE,WRK1,1,G,1)

       RETURN
       END

       SUBROUTINE TSDFLT(M,N,ITNLIM,JACFLG,GRADTL,STEPTL,FTOL,METHOD,
     +                   GLOBAL,STEPMX,DLT,TYPX,TYPF,IPR,MSG)

       INTEGER M,N,ITNLIM,JACFLG,METHOD,GLOBAL,MSG,IPR
       DOUBLE PRECISION GRADTL,STEPTL,FTOL,STEPMX,DLT
       DOUBLE PRECISION TYPX(N),TYPF(M)

C*********************************************************************
C THIS ROUTINE SETS DEFAULT VALUES FOR EACH INPUT VARIABLE TO THE
C NONLINEAR EQUATION ALGORITHM.
C*********************************************************************
C
C       SUBPROGRAMS CALLED:
C
C       TENSOLVE      ...  TSDLOD
C       UNCMIN        ...  DPMEPS
C
C**********************************************************************

       DOUBLE PRECISION EPS,DPMEPS,ONE,TWO,THREE,THOUS
       DATA ONE,TWO,THREE,THOUS/1.0D0,2.0D0,3.0D0,1000.0D0/

       JACFLG = 0
       EPS = DPMEPS()
       GRADTL = EPS**(ONE/THREE)
       STEPTL = EPS**(TWO/THREE)
       FTOL = EPS**(TWO/THREE)
       ITNLIM = 150
       METHOD = 1
       GLOBAL = 0
       STEPMX = THOUS
       DLT = -ONE
       MSG = 0
       IPR = 6
       CALL TSDLOD(N,ONE,TYPX,1)
       CALL TSDLOD(M,ONE,TYPF,1)

       RETURN
       END

       SUBROUTINE TSDUMJ(X,AJA,NR,M,N)

       INTEGER NR, M, N
       DOUBLE PRECISION AJA(NR,N),X(N)

C*********************************************************************
C THIS IS A DUMMY ROUTINE TO PREVENT UNSATISFIED EXTERNAL DIAGNOSTIC
C WHEN SPECIFIC ANALYTIC JACOBIAN IS NOT SUPPLIED.
C*********************************************************************
C
C      INPUT PARAMETERS:
C      -----------------
C
C      X   : POINT AT WHICH JACOBIAN IS EVALUATED
C      AJA : JACOBIAN MATRIX
C      NR  : LEADING DIMENSION OF AJA
C      M,N : DIMENSIONS OF PROBLEM
C
C***********************************************************************

       RETURN
       END

       FUNCTION TSFAFA(ANLS,FQ,ADT,AG,CONST1,CONST2,ALPHA,DLT,
     +                 NR,M,N,P,NWTAKE,IERR,VN)

       INTEGER NR,M,N,P,IERR
       DOUBLE PRECISION ALPHA,DLT,TSFAFA
       DOUBLE PRECISION ADT(N),AG(N),CONST1(P),CONST2(P)
       DOUBLE PRECISION FQ(M),VN(M),ANLS(NR,P)
       LOGICAL NWTAKE

C********************************************************************
C THIS FUNCTION COMPUTES || F + J*D + 0.5*A*D**2 ||**2 IN THE
C L2 NORM SENS, WHERE D = ALPHA*DT + SQRT(DLT**2-ALPHA**2).
C********************************************************************
C
C
C   INPUT PARAMETERS
C   ----------------
C
C       ANLS   : TENSOR TERM MATRIX
C       FQ     : FUNCTION VALUE AT CURRENT ITERATE MULTIPLIED BY
C                ORTHOGONAL MATRICES
C       ADT    : JACOBIAN MATRIX TIMES DT
C        AG    : JACOBIAN MATRIX TIMES GBAR (SEE SUBROUTINE TS2DTR)
C       CONST1 : SHAT-TRANS TIMES DT
C       CONST2 : SHAT-TRANS TIMES GBAR
C       ALPHA  : POINT AT WHICH TO EVALUATE THE FUNCTION TSFAFA
C       DLT    : CURRENT TRUST RADIUS
C       NR     : LEADING DIMENSION OF THE JACOBIAN
C       M,N    : DIMENSIONS OF THE PROBLEM
C       P      : COLUMN DIMENSION OF THE MATRICES SHAT AND ANLS
C       NWTAKE : LOGICAL VARIABLE WITH THE FOLLOWING MEANINGS:
C                NWTAKE = .TRUE.  : STANDARD STEP TAKEN
C                NWTAKE = .FALSE. : TENSOR STEP TAKEN
C       IERR   : RETURN CODE FROM QRP FACTORIZATION ROUTINE:
C                IERR = 0 : NO SINGULARITY OF JACOBIAN DETECTED
C                IERR = 1 : SINGULARITY OF JACOBIAN DETECTED
C
C
C       OUTPUT PARAMETERS
C       -----------------
C
C       VN     : F + J*D + 0.5*A*D**2
C       TSFAFA :  || F + J*D + 0.5*A*D**2 ||**2
C                WHERE D = ALPHA*DT + SQRT(DLT**2-ALPHA**2)
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DDOT
C       TENSOLVE      ...  TSMAFA
C
C********************************************************************

        INTEGER LEN
        DOUBLE PRECISION DDOT
        DOUBLE PRECISION HALF
        DATA HALF/0.50D0/

        CALL TSMAFA(ANLS,FQ,ADT,AG,CONST1,CONST2,ALPHA,DLT,
     +              NR,M,N,P,NWTAKE,IERR,VN)

        LEN = M
        IF(IERR.GT.0) LEN = M + N

        TSFAFA = HALF*DDOT(LEN,VN,1,VN,1)

        RETURN
        END

       SUBROUTINE TSFDFJ(XC,FC,NR,M,N,EPSM,FVEC,FHAT,AJA)

       INTEGER NR,M,N
       DOUBLE PRECISION EPSM
       DOUBLE PRECISION AJA(NR,N),FHAT(M),XC(N),FC(M)
       EXTERNAL FVEC

C***********************************************************************
C THIS ROUTINE COMPUTES THE FINITE DIFFERENCE JACOBIAN AT THE CURRENT
C ITERATE XC.
C***********************************************************************
C
C    INPUT PARAMETERS :
C    ----------------
C
C    XC   : CURRENT ITERATE
C    FC   : FUNCTION VALUE AT XC
C    NR   : LEADING DIMENSION OF AJA
C    M,N  : DIMENSIONS OF PROBLEM
C    EPSM : MACHINE PRECISION
C    FVEC : SUBROUTINE TO EVALUATE THE USER'S FUNCTION
C    FHAT : WORKSPACE
C
C    OUTPUT PARAMETERS :
C    --------------------
C
C    AJA : FINITE DIFFERENCE JACOBIAN AT XC
C
C    SUBPROGRAMS CALLED:
C
C    USER   ...  FVEC
C
C***********************************************************************

       INTEGER I,J
       DOUBLE PRECISION NDIGIT,RNOISE,STEPSZ,XTMPJ
       DOUBLE PRECISION SQRTR,RSTPSZ,ONE,TEN
       INTRINSIC ABS,MAX,SQRT
       DATA ONE,TEN/1.0D0,10.0D0/

       NDIGIT = -LOG10(EPSM)
       RNOISE = MAX(TEN**(-NDIGIT),EPSM)
       SQRTR = SQRT(RNOISE)

       DO 20 J = 1,N
          STEPSZ = SQRTR*MAX(ABS(XC(J)),ONE)
          XTMPJ = XC(J)
          XC(J) = XTMPJ+STEPSZ
          CALL FVEC(XC,FHAT,M,N)
          XC(J) = XTMPJ
          RSTPSZ = ONE/STEPSZ
          DO 10 I = 1,M
             AJA(I,J) = (FHAT(I)-FC(I))*RSTPSZ
 10       CONTINUE
 20    CONTINUE

       RETURN
       END

       SUBROUTINE TSFRMT(SHAT,S,AJA,FV,FN,MAXM,MAXN,MAXP,M,N,P,IDP,
     +                   AM,X,B,SCALE,A)

       INTEGER MAXM,MAXN,MAXP,M,N,P
       INTEGER IDP(P)
       DOUBLE PRECISION A(MAXM,P),SHAT(MAXN,P),S(MAXN,P),AJA(MAXM,N)
       DOUBLE PRECISION FV(MAXM,P),FN(M),AM(MAXP,P),X(P),B(P),SCALE(P)

C*********************************************************************
C THIS ROUTINE FORM THE TENSOR TERM MATRIX OF THE TENSOR MODEL.
C*********************************************************************
C
C       INPUT PARAMETERS :
C       ----------------
C
C       SHAT: MATRIX OF PAST LINEARLY INDEPENDENT DIRECTIONS
C       S   : MATRIX OF PREVIOUS DIRECTIONS
C       AJA : JACOBIAN MATRIX AT CURRENT ITERATE
C       FV  : MATRIX OF PAST FUNCTION VALUES
C       FN  : FUNCTION VALUE AT CURRENT ITERATE
C       MAXM: LEADING DIMENSION OF AJA, ANLS, AND FV
C       MAXN: LEADING DIMENSION OF S AND SHAT
C       MAXP: LEADING DIMENSION OF AM
C       M   : ROW DIMENSION OF MATRICES A,FV,AND AJA
C       N   : COLUMN DIMENSION OF JACOBIAN MATRIX
C       P   : COLUMN DIMENSION OF MATRIX SHAT
C       IDP : VECTOR WHICH KEEPS TRACK OF LINEARLY INDEPENDENT
C             DIRECTION POSITIONS WITHIN THE MATRIX S
C       AM,X,B,SCALE,: WORKSPACE
C
C       OUTPUT PARAMETERS :
C       ------------------
C
C       A   : TENSOR TERM MATRIX
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DDOT,DNRM2,DSCAL
C       UNCMIN        ...  CHOLDC,LLTSLV
C
C*********************************************************************

       INTEGER I,J,JJ
       DOUBLE PRECISION SUM,SC,TOL,DIAGMX,ADDMAX
       DOUBLE PRECISION ZERO,ONE,TWO
       DOUBLE PRECISION DDOT,DNRM2
       DATA ZERO,ONE,TWO/0.0D0,1.0D0,2.0D0/

c scale the matrix SHAT and save scaling in SCALE

       DO 10 J = 1,P
          SC = ONE/DNRM2(N,SHAT(1,J),1)
          CALL DSCAL(N,SC,SHAT(1,J),1)
          SCALE(J) = SC**2
 10    CONTINUE

c form the matrix AM = (Si Sj)**2

       DO 30 J = 1,P
          JJ = IDP(J)
          DO 20 I = 1,P
             AM(I,J) = DDOT(N,S(1,IDP(I)),1,S(1,JJ),1)**2
 20      CONTINUE
 30   CONTINUE

c scale the matrix AM

       DO 50 I = 1,P
          DO 40 J = 1,P
             AM(I,J) = SCALE(I)*SCALE(J)*AM(I,J)
 40      CONTINUE
 50   CONTINUE

c perform a Cholesky decomposition of AM

       TOL = ZERO
       DIAGMX = ZERO
       CALL CHOLDC(MAXP,P,AM,DIAGMX,TOL,ADDMAX)

c form the tensor term matrix A

       DO 70 I = 1,M
          DO 60 J = 1,P
             JJ = IDP(J)
             SUM = DDOT(N,AJA(I,1),MAXM,S(1,JJ),1)
             B(J) = TWO*(FV(I,JJ) - FN(I) - SUM)
             B(J) = SCALE(J)*B(J)
 60       CONTINUE

c solve AM*X = B

          CALL LLTSLV(MAXP,P,AM,X,B)

c copy X into row i of A

          CALL DCOPY(P,X,1,A(I,1),MAXM)

 70    CONTINUE

       RETURN
       END

       SUBROUTINE TSFSCL(X,DX,DF,M,N,FVEC,F)

       INTEGER M,N
       DOUBLE PRECISION X(N),DX(N),F(M),DF(M)
       EXTERNAL FVEC

C********************************************************************
C THIS ROUTINE EVALUATES THE FUNCTION AT THE CURRENT ITERATE X THEN
C SCALES ITS VALUE.
C********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       X  : CURRENT ITERATE
C       DX : DIAGONAL SCALING MATRIX FOR X
C       DF : DIAGONAL SCALING MATRIX FOR F
C       M,N :  DIMENSIONS OF PROBLEM
C       FVEC : SUBROUTINE TO EVALUATE FUNCTION
C
C
C       OUTPUT PARAMETERS :
C       -----------------
C
C       F  : SCALED FUNCTION VALUE AT CURRENT ITERATE X
C
C       SUBPROGRAMS CALLED:
C
C       TENSOLVE      ...  TSUNSX,TSSCLF,TSSCLX
C       USER          ...  FVEC
C
C********************************************************************

       CALL TSUNSX(X,DX,N)
       CALL FVEC(X,F,M,N)
       CALL TSSCLF(F,DF,M)
       CALL TSSCLX(X,DX,N)

       RETURN
       END

       SUBROUTINE TSFSLV(L,B,NR,M,N,Y)

       INTEGER NR,M,N
       DOUBLE PRECISION B(N),L(NR,N),Y(N)

C********************************************************************
C THIS ROUTINE DOES A FORWARD SOLVE.
C********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       L   : THE TRANSPOSE OF THE UPPER TRIANGULAR MATRIX OBTAINED
C             FROM A QR FACTORIZATION OF AN M BY N MATRIX A. DIAG(L)
C             IS STORED IN ROW M+2. THIS IS THE STORAGE SCHEME USED
C             IN STEWART, G. W., III(1973) "INTRODUCTION TO MATRIX
C             COMPUTATION", ACADEMIC PRESS,NEW YORK
C       B   : RIGHT HAND SIDE
C       NR  : LEADING DIMENSION OF MATRIX A
C        M  : ROW DIMENSION OF MATRIX A
C        N  : COLUMN DIMENSION OF MATRIX A
C
C       OUTPUT PARAMETERS :
C       ------------------
C
C        Y  : VECTOR SOLUTION ON EXIT
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DDOT
C
C*********************************************************************

       INTEGER J
       DOUBLE PRECISION S
       DOUBLE PRECISION DDOT

c solve L Y = B

       Y(1) = B(1) / L(M+2,1)
       IF(N .GT. 1) THEN
          S = L(1,2) * Y(1)
          Y(2) = (B(2) - S) / L(M+2,2)
          DO 10 J = 3,N
             S = DDOT(J-1,L(1,J),1,Y,1)
             Y(J) = (B(J) - S) / L(M+2,J)
 10       CONTINUE
       ENDIF

       RETURN
       END

      SUBROUTINE TSJMUV(ITN,METHOD,V,CURPOS,PIVOT,PBAR,AJA,SHAT,
     +                  FLAG,IERR,MAXM,MAXN,M,N,P,WRK1,WRK2,VN,AV)

      INTEGER MAXM,MAXN,M,N,P,IERR,ITN,METHOD,FLAG
      INTEGER CURPOS(N),PIVOT(N),PBAR(N)
      DOUBLE PRECISION V(N),WRK1(N),WRK2(N),VN(M),AJA(MAXM,N)
      DOUBLE PRECISION AV(N),SHAT(MAXN,P)

C****************************************************************
C THIS ROUTINE CALCULATES THE PRODUCT JACOBIAN TIMES A VECTOR.
C****************************************************************
C
C       INPUT PARAMETERS
C       ----------------
C
C        ITN    : CURRENT ITERATION NUMBER
C        METHOD : METHOD TO BE USED
C        V      : VECTOR TO BE MULTIPLIED BY AJA
C        CURPOS : PIVOT VECTOR (USED DURING THE FACTORIZATION OF AJA
C                 FROM COLUMN 1 TO N-P)
C        PIVOT  : PIVOT VECTOR (USED DURING THE FACTORIZATION OF AJA
C                 FROM COLUMN N-P+1 TO N)
C        PBAR   : PIVOT VECTOR (USED DURING THE FACTORIZATION OF AJA
C                 IF IT IS SINGULAR
C        AJA    : JACOBIAN MATRIX AT CURRENT ITERATE
C        SHAT   : MATRIX OF LINEARLY INDEPENDENT DIRECTIONS AFTER
C                 A QL FACTORIZATION
C        FLAG   : RETURN CODE WITH THE FOLLOWING MEANINGS:
C                FLAG = 0 : NO SINGULARITY DETECTED DURING FACTORIZATION
C                           OF THE JACOBIAN FROM COLUMN 1 TO N
C                FLAG = 1 : SINGULARITY DETECTED DURING FACTORIZATION
C                           OF THE JACOBIAN FROM COLUMN 1 TO N-P
C                FLAG = 2 : SINGULARITY DETECTED DURING FACTORIZATION
C                           OF THE JACOBIAN FROM COLUMN N-P+1 TO N
C        IERR   : RETURN CODE FROM QRP FACTORIZATION ROUTINE:
C                IERR = 0 : NO SINGULARITY OF JACOBIAN DETECTED
C                IERR = 1 : SINGULARITY OF JACOBIAN DETECTED
C        MAXM   : LEADING DIMENSION OF AJA
C        MAXN   : LEADING DIMENSION OF SHAT
C        M,N    : DIMENSIONS OF THE PROBLEM
C        P      : COLUMN DIMENSION OF THE MATRICES SHAT AND ANLS
C        WRK1,WRK2,VN : WORKSPACE VECTORS
C
C       OUTPUT PARAMETERS
C       -----------------
C
C        AV     : JACOBIAN TIMES V
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DCOPY
C       TENSOLVE      ...  TSPRMV,TSQMLV,TSUTMD
C
C **********************************************************************

      INTEGER LEN
      IF(ITN .EQ. 1 .OR. METHOD .EQ. 0) THEN
         CALL TSPRMV(WRK1,V,PIVOT,N,1)
         IF(IERR .EQ. 1) THEN
           CALL TSPRMV(WRK2,WRK1,PBAR,N,1)
           CALL DCOPY(N,WRK2,1,WRK1,1)
         ENDIF
      ELSEIF(N .EQ. 1) THEN
           VN(1) = V(1)
      ELSE
           CALL TSQMLV(MAXN,N,P,SHAT,V,VN,.FALSE.)
           CALL TSPRMV(WRK2,VN,CURPOS,N,1)
           IF(FLAG .EQ. 0) THEN
              CALL TSPRMV(WRK1,WRK2,PIVOT,N,1)
           ELSEIF(FLAG .EQ. 1) THEN
              CALL TSPRMV(WRK1,WRK2,PBAR,N,1)
           ELSEIF(FLAG .EQ. 2 ) THEN
              CALL TSPRMV(WRK1,WRK2,PIVOT,N,1)
              CALL TSPRMV(WRK2,WRK1,PBAR,N,1)
              CALL DCOPY(N,WRK2,1,WRK1,1)
           ENDIF
      ENDIF

      LEN = M
      IF(IERR .GT. 0) LEN = M + N

      CALL TSUTMD(AJA,WRK1,MAXM,LEN,N,AV)

      RETURN
      END

       SUBROUTINE TSJQTP(Q,MAXM,MAXN,N,M,P,WRK1,WRK2,AJA)

       INTEGER MAXM,MAXN,N,M,P
       DOUBLE PRECISION AJA(MAXM,N),Q(MAXN,P),WRK1(N),WRK2(N)

C**********************************************************************
C THIS ROUTINE GETS J*(Q-TRANS) BY COMPUTING EACH ROW OF THE
C RESULTING MATRIX AS FOLLOWS : (J*Q-TRANS)I-TH ROW<--Q*(J)I-TH ROW.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       Q    : RESULTING MATRIX FROM A QL FACTORIZATION
C       MAXM : LEADING DIMENSION OF AJA
C       MAXN : LEADING DIMENSION OF Q
C       M,N  : DIMENSIONS OF PROBLEM
C       P    : COLUMN DIMENSION OF MATRIX Q
C       WRK1,WRK2: WORKING VECTOR
C
C       INPUT-OUTPUT PARAMETERS :
C       ------------------------
C
C       AJA : JACOBIAN MATRIX ON ENTRY AND JACOBIAN MULTIPLIED BY THE
C             ORTHOGONAL MATRIX Q ON EXIT
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DCOPY
C       TENSOLVE      ...  TSQMLV
C
C**********************************************************************

       INTEGER I

       DO 30 I = 1,M

c copy the i-th row of AJA into WRK1

          CALL DCOPY(N,AJA(I,1),MAXM,WRK1,1)

          CALL TSQMLV(MAXN,N,P,Q,WRK1,WRK2,.FALSE.)

c form the i-th row of AJA*(Q-trans)

          CALL DCOPY(N,WRK2,1,AJA(I,1),MAXM)

 30    CONTINUE

       RETURN
       END

       SUBROUTINE TSLMIN(XC,XP,P1,Q,ANLS,FQ,ADT,AG,CONST1,CONST2,
     +                   DLT,NR,M,N,P,NWTAKE,IERR,TOL,VN,VNP,VNS,XPLUS)

        INTEGER NR,M,N,P,IERR
        DOUBLE PRECISION XC,XP,XPLUS,P1,Q,DLT,TOL
        DOUBLE PRECISION ADT(N),AG(N),VN(M),VNP(M),VNS(M)
        DOUBLE PRECISION ANLS(NR,P),FQ(M),CONST1(P),CONST2(P)
        LOGICAL NWTAKE

C***********************************************************************
C THIS ROUTINE FINDS A LOCAL MINIMIZER OF A ONE-VARIABLE FUNCTION IN AN
C INTERVAL [XC XP].
C***********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       XC,XP  : LOWER AND UPPER BOUND OF INTERVAL IN WHICH THE SEARCH
C                IS PERFORMED
C       P1,Q   : FIRST DERIVATIVES OF THE ONE-VARIABLE FUNCTION
C       ANLS   : TENSOR TERM MATRIX
C       FQ     : FUNCTION VALUE AT CURRENT ITERATE MULTIPLIED BY
C                ORTHOGONAL MATRICES
C       ADT    : JACOBIAN TIMES THE STEP DT (SEE SUBROUTINE TS2DTR)
C       AG     : JACOBIAN TIMES THE GRADIENT G (SEE SUBROUTINE TS2DTR)
C       CONST1 : SHAT-TRANS * DT  (SEE SUBROUTINE TS2DTR)
C       CONST2 : SHAT-TRANS * GBAR (SEE SUBROUTINE TS2DTR)
C       DLT    : TRUST RADIUS
C       NR     : LEADING DIMENSION OF ANLS MATRIX
C       M,N    : DIMENSIONS OF PROBLEM
C       P      : COLUMN DIMENSION OF MATRIX ANLS
C       NWTAKE : LOGICAL VARIABLE WITH THE FOLLOWING MEANINGS:
C                NWTAKE = .TRUE.  : STANDARD STEP TAKEN
C                NWTAKE = .FALSE. : TENSOR STEP TAKEN
C       IERR   : RETURN CODE FROM QRP FACTORIZATION ROUTINE:
C                IERR = 0 : NO SINGULARITY OF JACOBIAN DETECTED
C                IERR = 1 : OTHERWISE
C       TOL    : SMALL TOLERANCE
C       VN,VNP,VNS : WORKING VECTORS
C
C
C       OUTPUT PARAMETERS :
C       -----------------
C
C       XPLUS  :  LOCAL MINIMIZER OF THE ONE-VARIABLE FUNCTION
C
C       SUBPROGRAMS CALLED :
C
C       TENSOLVE      ...  TSMSDA,TSFAFA,TSLMSP,TSMFDA
C
C***********************************************************************

        INTEGER ITERCD,RETCD,ITNCNT
        DOUBLE PRECISION ALEFT,ARIGHT,T,E,TSMSDA,S,SINIT,TSFAFA,TSMFDA
        DOUBLE PRECISION ZERO,OTT,TWO,SMALL
        LOGICAL SKIP
        INTRINSIC ABS,MIN,MAX
        DATA ZERO,OTT,TWO,SMALL/0.0D0,1.0D-4,2.0D0,2.0D-20/

        RETCD = 0
        ALEFT = MIN(XC,XP)
        ARIGHT = MAX(XC,XP)
        ITNCNT = 0
        T = ABS(XC-XP)
        SKIP = .FALSE.

c compute the second derivative value at the current point

        E = TSMSDA(ANLS,FQ,ADT,AG,CONST1,CONST2,XC,DLT,
     +             NR,M,N,P,NWTAKE,IERR,SKIP,VN,VNP,VNS)

 10     IF(E.GT.ZERO) THEN
           S = -P1/E
           IF(ABS(S).GT.TWO*T) THEN
              IF (S.LT.ZERO) THEN
                  S = -TWO*T
              ELSE
                  S = TWO*T
              ENDIF
           ENDIF
        ELSE
           IF (P1.GT.ZERO) THEN
               S = -T
           ELSE
               S = T
           ENDIF
        ENDIF

        IF(XC+S.GT.ARIGHT) S = ARIGHT-XC
        IF(XC+S.LT.ALEFT)  S = ALEFT-XC
        SINIT = ABS(S)

 20     CONTINUE

c compute a next iterate XPLUS

        IF (TSFAFA(ANLS,FQ,ADT,AG,CONST1,CONST2,XC+S,DLT,
     +      NR,M,N,P,NWTAKE,IERR,VN).GT.Q + OTT*S*P1) THEN
            S = S/2
            IF(ABS(S).LT.SMALL*SINIT.OR.S.EQ.ZERO) THEN
               RETCD = 1
            ELSE
               GO TO 20
            ENDIF
        ENDIF

        XPLUS = XC+S
        ITNCNT = ITNCNT+1

c check stopping criteria

        CALL TSLMSP(XC,XPLUS,ITNCNT,RETCD,ITERCD,ANLS,ADT,AG,
     +              CONST1,CONST2,DLT,NR,M,N,P,NWTAKE,IERR,TOL,VN,VNP)

        IF(ITERCD.GT.0) RETURN

c update XC

        XC = XPLUS

c compute function and derivative values at the new point

        Q = TSFAFA(ANLS,FQ,ADT,AG,CONST1,CONST2,XC,DLT,
     +             NR,M,N,P,NWTAKE,IERR,VN)
        P1 = TSMFDA(ANLS,ADT,AG,CONST1,CONST2,XC,DLT,
     +              NR,M,N,P,NWTAKE,IERR,VN,VNP)
        SKIP = .TRUE.
        E = TSMSDA(ANLS,FQ,ADT,AG,CONST1,CONST2,XC,DLT,
     +             NR,M,N,P,NWTAKE,IERR,SKIP,VN,VNP,VNS)
        GO TO 10

        END

       SUBROUTINE TSLMSP(XC,XP,ITNCNT,RETCD,ITERCD,ANLS,ADT,AG,CONST1,
     +                  CONST2,DLT,NR,M,N,P,NWTAKE,IERR,TOL,VN,VNP)

       INTEGER NR,M,N,P,IERR,RETCD,ITERCD,ITNCNT
       DOUBLE PRECISION XC,XP,DLT,TOL
       DOUBLE PRECISION ADT(N),AG(N),CONST1(P)
       DOUBLE PRECISION CONST2(P),VN(M),VNP(M),ANLS(NR,P)
       LOGICAL NWTAKE

C***********************************************************************
C THIS ROUTINE CHECKS THE STOPPING CRITERIA FOR A LOCAL MINIMIZER.
C***********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       XC       : CURRENT ITERATE (FROM SEARCH SUBROUTINE)
C       XP       : NEXT ITERATE (FROM SEARCH SUBROUTINE)
C       ITNCNT   : ITERATION LIMIT
C       RETCD    : RETURN CODE FROM LINE SEARCH
C       DLT      : TRUST RADIUS
C       AJA      : JACOBIAN AT THE CURRENT ITERATE
C       NR       : LEADING DIMENSION OF THE JACOBIAN MATRIX
C       M,N      : DIMENSIONS OF THE PROBLEM
C       P        : COLUMN DIMENSION OF THE MATRICES SHAT AND ANLS
C       NWTAKE   : LOGICAL VARIABLE WITH THE FOLLOWING MEANINGS :
C                  NWTAKE = .TRUE.  : STANDARD STEP TAKEN
C                  NWTAKE = .FALSE. : TENSOR STEP TAKEN
C       IERR     : RETURN CODE FROM THE QRP FACTORIZATION ROUTINE :
C                  IERR = 0 : NO SINGULARITY OF JACOBIAN DETECTED
C                  IERR = 1 : OTHERWISE
C       TOL      : SMALL TOLERANCE
C       METHOD   : METHOD TO USE
C                = 0   : STANDARD METHOD USED
C                = 1   : TENSOR METHOD USED
C       VN,VNP  : WORKING VECTORS
C
C
C       OUTPUT PARAMETERS :
C       ------------------
C
C       ITERCD  : RETURN CODE WITH FOLLOWING MEANINGS :
C                 ITERCD = 1 : FIRST DERIVATIVE AT THE CURRENT POINT
C                              CLOSE TO 0
C                 ITERCD = 2 : SUCCESSIVE ITERATES WITHIN TOLERANCE
C                 ITERCD = 3 : LINE SEARCH FAILED TO LOCATE A POINT
C                              LOWER THAT THE CURRENT POINT
C                 ITERCD = 4 : ITERATION LIMIT EXCEEDED
C
C***********************************************************************

        DOUBLE PRECISION TSMFDA,GRDT,ZERO
        INTRINSIC ABS,SQRT
        DATA ZERO/0.0D0/

        GRDT = SQRT(TOL)
        ITERCD = 0

        IF(RETCD.EQ.1) THEN
           ITERCD = 3
        ELSEIF(ABS(TSMFDA(ANLS,ADT,AG,CONST1,CONST2,XP,DLT,
     +         NR,M,N,P,NWTAKE,IERR,VN,VNP)) .LT. GRDT) THEN
               ITERCD = 1
        ELSEIF(XP.NE.ZERO .AND. ABS(XP-XC)/ABS(XP).LE.TOL) THEN
               ITERCD = 2
        ELSEIF(ITNCNT.GE.150) THEN
               ITERCD = 4
        ENDIF

        RETURN
        END

       SUBROUTINE TSLSCH(M,N,XC,D,G,STEPTL,DX,DF,FVEC,
     +                   MXTAKE,STEPMX,XP,FP,FCNORM,FPNORM,RETCD)

       INTEGER M,N,RETCD
       DOUBLE PRECISION STEPTL,FCNORM,FPNORM
       DOUBLE PRECISION XC(N)
       DOUBLE PRECISION D(N),G(N),XP(N),FP(M)
       DOUBLE PRECISION DX(N),DF(M),STEPMX
       LOGICAL MXTAKE
       EXTERNAL FVEC

C**********************************************************************
C THIS ROUTINE FINDS A NEXT ITERATE USING A STANDARD LINE SEARCH METHOD.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       M,N : DIMENSIONS OF PROBLEM
C       XC  : CURRENT ITERATE
C       D   : SEARCH DIRECTION
C       G   : GRADIENT AT CURRENT ITERATE
C       STEPTL : RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES
C                ARE CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM
C       DX  : DIAGONAL SCALING MATRIX FOR X
C       DF  : DIAGONAL SCALING MATRIX FOR F
C       FVEC: SUBROUTINE TO EVALUATE THE FUNCTION
C       STEPMX: MAXIMUM ALLOWABLE STEP SIZE
C
C
C       OUTPUT PARAMETERS :
C       -----------------
C
C       MXTAKE: BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED
C       XP : NEXT ITARATE
C       FP : FUNCTION VALUE AT NEXT ITERATE
C       FCNORM : 0.5 * || F(XC) ||**2
C       FPNORM : 0.5 * || F(XP) ||**2
C       RETCD : RETURN CODE WITH THE FOLLOWING MEANING :
C                RETCD = 0 : SATISFACTORY LOCATION OF A NEW ITERATE
C                RETCD = 1 : NO SATISFACTORY POINT FOUND SUFFICIENTLY
C                            DISTINCT FROM X
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DDOT,DNRM2
C       TENSOLVE      ...  TSFSCL
C       USER          ...  FVEC
C
C**********************************************************************

       INTEGER I
       DOUBLE PRECISION ALPHA,SLOPE,RELENG
       DOUBLE PRECISION TEMP1,TEMP2,ALMDA,TEMP,ALMDAT,ALMDAM
       DOUBLE PRECISION SLN,SCL
       DOUBLE PRECISION DDOT,DNRM2
       DOUBLE PRECISION ZERO,TENTH,HALF,Z99,ONE,TWO,TEN
       INTRINSIC ABS
       PARAMETER (ALPHA = 1.0D-4)
       DATA ZERO,TENTH,HALF,Z99,ONE,TWO,TEN/0.0D0,0.10D0,0.50D0,0.99D0,
     + 1.0D0,2.0D0,10.0D0/

       MXTAKE = .FALSE.
       SLN = DNRM2(N,D,1)
       IF(SLN .GT. STEPMX) THEN

c step longer than maximum allowed

          SCL = STEPMX/SLN
          CALL DSCAL(N,SCL,D,1)
          SLN = STEPMX
       ENDIF

c compute SLOPE  =  G-trans * D

       SLOPE = DDOT(N,G,1,D,1)

c initialization of RETCD

       RETCD = 0

c compute the smallest value allowable for the damping
c parameter ALMDA, i.e, ALMDAM

       RELENG = ZERO
       DO 20 I = 1,N
         TEMP1 = MAX(ABS(XC(I)), ONE)
         TEMP2 = ABS(D(I))/TEMP1
         RELENG = MAX(RELENG,TEMP2)
 20    CONTINUE
       ALMDAM = STEPTL/RELENG
       ALMDA = ONE

 40    CONTINUE

c compute the next iterate XP

       DO 50 I = 1,N
          XP(I) = XC(I)+ALMDA*D(I)
 50    CONTINUE

c evaluate the objective function at XP and its residual

       CALL TSFSCL(XP,DX,DF,M,N,FVEC,FP)

       FPNORM = HALF*DNRM2(M,FP,1)**2

c test whether the full step produces enough decrease in
c the l2 norm of the objective function. If not update ALMDA
c and compute a new step

       IF (FPNORM.GT.(FCNORM + (ALPHA* ALMDA * SLOPE))) THEN
          ALMDAT = ((-ALMDA**2)*SLOPE)/(TWO*(FPNORM-FCNORM-ALMDA*SLOPE))
          TEMP = ALMDA/TEN
          ALMDA = MAX(TEMP,ALMDAT)
          IF(ALMDA.LT.ALMDAM) THEN
             RETCD = 1
             RETURN
          ENDIF
          GO TO 40
       ELSE
          IF(ALMDA.EQ.TENTH .AND. SLN.GT.Z99*STEPMX) MXTAKE=.TRUE.
       ENDIF

       RETURN
       END

       SUBROUTINE TSMAFA(ANLS,F,ADT,AG,CONST1,CONST2,ALPHA,DLT,
     +                   NR,M,N,P,NWTAKE,IERR,VN)

       INTEGER NR,M,N,P,IERR
       DOUBLE PRECISION ALPHA,DLT
       DOUBLE PRECISION ADT(N),AG(N),CONST1(P)
       DOUBLE PRECISION CONST2(P),F(M),VN(M),ANLS(NR,P)
       LOGICAL NWTAKE

C***********************************************************************
C THIS ROUTINE COMPUTES THE VECTOR VN = F(XC) + J(XC)*D + 0.5*A*D**2,
C WHERE D = ALPHA*DT + SQRT(DLT**2-ALPHA**2).
C***********************************************************************
C
C
C      INPUT PARAMETERS :
C      -----------------
C
C      ANLS  : TENSOR TERM MATRIX
C       ADT  : JACOBIAN MATRIX TIMES DT (SEE SUBROUTINE TS2DTR)
C        AG  : JACOBIAN MATRIX TIMES GBAR (SEE SUBROUTINE TS2DTR)
C      CONST1: SHAT-TRANS * DT (SEE SUBROUTINE TS2DTR)
C      CONST2: SHAT-TRABS * GBAR (SEE SUBROUTINE TS2DTR)
C      ALPHA : POINT AT WHICH DERIVATIVE IS EVALUATED
C        DLT : CURRENT TRUST RADIUS
C        NR  : LEADING DIMENSION OF ANLS
C        M,N : DIMENSIONS OF THE PROBLEM
C        P   : COLUMN DIMENSION OF THE MATRIX ANLS
C        NWTAKE : LOGICAL VARIABLE WITH THE FOLLOWING MEANINGS
C                 NWTAKE = .TRUE.  : STANDARD STEP TAKEN
C                 NWTAKE = .FALSE. : TENSOR STEP TAKEN
C       IERR : RETURN CODE FROM THE QRP FACTORIZATION ROUTINE :
C              IERR = 0 : NO SINGULARITY OF JACOBIAN DETECTED
C              IERR = 1 : SINGULARITY OF JACOBIAN DETECTED
C
C       OUTPUT PARAMETERS :
C       -------------------
C
C       VN  : F + J*D + 0.5*A*D**2, WHERE
C             D = ALPHA*DT + SQRT(DLT**2-ALPHA**2)
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DAXPY
C       TENSOLVE      ...  TSDLOD
C
C*******************************************************************

       INTEGER I,J,LEN
       DOUBLE PRECISION EXPR,CONST,ZERO,HALF
       INTRINSIC SQRT
       DATA ZERO,HALF/0.0D0,0.50D0/

       EXPR = SQRT(DLT**2 - ALPHA**2)
       DO 10 I = 1,N
          VN(I) = ALPHA*ADT(I) + EXPR*AG(I)
 10    CONTINUE

       CALL TSDLOD (M,ZERO,VN(N+1),1)

       LEN = M
       IF(IERR .GT. 0) LEN = M + N

       DO 30 I = 1, LEN
          VN(I) = VN(I) + F(I)
 30    CONTINUE

       IF(NWTAKE) RETURN
       DO 70 J = 1,P
          CONST = HALF*(ALPHA*CONST1(J) + EXPR*CONST2(J))**2
          CALL DAXPY(LEN,CONST,ANLS(1,J),1,VN,1)
 70    CONTINUE

       RETURN
       END

       SUBROUTINE TSMDLS(AJA,SHAT,ANLS,XC,M,N,MAXM,MAXN,P,DT,G,
     +                   DX,DF,FVEC,METHOD,STEPTL,GLOBAL,STEPMX,
     +                   EPSM,FQ,WRK1,WRK2,WRK3,WRK4,DN,FQQ,PIVOT,
     +                   CURPOS,PBAR,MXTAKE,XP,FP,FCNORM,FPNORM,
     +                   ZERO1,RETCD,IERR)

        INTEGER M,N,MAXM,MAXN,P,METHOD,GLOBAL,ZERO1,RETCD,IERR
        INTEGER PIVOT(N),PBAR(N),CURPOS(N)
        DOUBLE PRECISION STEPTL,STEPMX,EPSM,FCNORM,FPNORM
        DOUBLE PRECISION AJA(MAXM,N),SHAT(MAXN,P),ANLS(MAXM,P)
        DOUBLE PRECISION XC(N),DT(N),G(N),DX(N),DF(M),FQ(M)
        DOUBLE PRECISION WRK1(M),WRK2(M),WRK3(M),WRK4(N)
        DOUBLE PRECISION DN(N),FQQ(M),XP(N),FP(M)
        LOGICAL MXTAKE
        EXTERNAL FVEC

C**********************************************************************
C THIS ROUTINE FINDS A NEXT ITERATE USING A LINE SEARCH METHOD.  IT
C TRIES THE FULL TENSOR STEP FIRST. IF THIS IS NOT SUCCESSFUL THEN
C IT COMPUTES THE STANDARD DIRECTION AND COMPUTES A STEP IN THAT
C DIRECTION. NEXT, IF THE TENSOR DIRECTION IS DESCENT, IT COMPUTES
C A STEP IN THE TENSOR DIRECTION.  THE ITERATE THAT PRODUCES
C THE LOWER RESIDUAL IS THE NEXT ITERATE FOR THE NONLINEAR ALGORITHM.
C**********************************************************************
C
C   INPUT PARAMETERS
C   ----------------
C
C       AJA    : JACOBIAN AT CURRENT ITERATE
C       SHAT   : MATRIX OF PAST LINEARLY INDEPENDENT DIRECTIONS
C                AFTER A QL FACORIZATION
C       ANLS   : TENSOR TERM MATRIX
C       XC     : CURRENT ITERATE
C       M,N    : DIMENSIONS OF THE PROBLEM
C       MAXM   : LEADING DIMENSION OF AJA AND ANLS
C       MAXN   : LEADING DIMENSION OF SHAT
C       P      : COLUMN DIMENSION OF THE MATRICES SHAT AND ANLS
C       DT     : TENSOR STEP
C       G      : GRADIENT AT CURRENT ITERATE
C       DX     : DIAGONAL SCALING MATRIX FOR X
C       DF     : DIAGONAL SCALING MATRIX FOR F
C       GBAR   : STEEPEST DESCENT DIRECTION (= -G)
C       METHOD : METHOD TO USE
C                = 0  : STANDARD METHOD USED
C                = 1  : TENSOR METHOD USED
C       STEPTL : STEP TOLERANCE
C       GLOBAL : GLOBAL STRATEGY USED
C                   =  0 : LINE SEARCH IS USED
C                   =  1 : 2-DIMENSIONAL TRUST REGION IS USED
C       STEPMX : MAXIMUM ALLOWABLE STEP SIZE
C       EPSM   : MACHINE PRECISION
C       FQ     : FUNCTION VALUE AT CURRENT ITERATE MULTIPLIED BY AN
C                ORTHOGOL MATRIX
C       WRK1,WRK2,WRK3,WRK4 : WORKING VECTORS
C
C
C       OUTPUT PARAMETERS
C       -----------------
C
C       DN     : NEWTON STEP
C       FQQ    : FQ MULTIPLIED BY AN ORTHOGONAL MATRIX
C       CURPOS : PIVOT VECTOR (USED DURING THE FACTORIZATION OF THE
C                JACOBIAN FROM COLUMN 1 TO N-P)
C       PIVOT  : PIVOT VECTOR (USED DURING THE FACTORIZATION OF THE
C                JACOBIAN FROM COLUMN N-P+1 TO N)
C       PBAR   : PIVOT VECTOR (USED DURING THE FACTORIZATION OF THE
C                JACOBIAN IF IT IS SINGULAR
C       MXTAKE : BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED
C       XP     : NEXT ITERATE
C       FP     : FUNCTION VALUE AT NEXT ITERATE
C       FCNORM :  0.5 * || F(XC) ||**2
C       FPNORM :  0.5 * || F(XP) ||**2
C        ZERO1 : FIRST ZERO COLUMN OF THE JACOBIAN IN CASE OF
C                SINGULARITY
C       RETCD  : RETURN CODE WITH THE FOLLOWING MEANING :
C                RETCD  =  0 : SATISFACTORY LOCATION OF A NEW ITERATE
C                RETCD  =  1 : NO SATISFACTORY POINT FOUND SUFFICIENTLY
C                              DISTINCT FROM X
C       IERR   : RETURN CODE FROM THE QRP FACTORIZATION ROUTINE
C                IERR = 0 : NO SINGULARITY OF JACOBIAN DETECTED
C                IERR = 1 : SINGULARITY OF JACOBIAN DETECTED
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DCOPY,DDOT,DNRM2
C       TENSOLVE      ...  TSFSCL,TSCPSS,TSLSCH
C
C***********************************************************************

        INTEGER I,FLAG,RETCD1
        DOUBLE PRECISION ALPHA,SLOPE,RELENG
        DOUBLE PRECISION TEMP1,TEMP2,ALMDA,RESNEW,F1N,DTNORM,GNORM
        DOUBLE PRECISION SLN,SCL
        DOUBLE PRECISION BETA,TEMP,ALMDAT,ALMDAM
        DOUBLE PRECISION DDOT,DNRM2
        DOUBLE PRECISION ZERO,TENTH,HALF,Z99,ONE,TWO,TEN
        INTRINSIC ABS
        PARAMETER (ALPHA = 1.0D-4)
        DATA ZERO,TENTH,HALF,Z99,ONE,TWO,TEN/0.0D0,0.10D0,0.50D0,0.99D0,
     +  1.0D0,2.0D0,10.0D0/

        MXTAKE = .FALSE.
        SLN = DNRM2(N,DT,1)
        IF(SLN .GT. STEPMX) THEN

c step longer than maximum allowed

           SCL = STEPMX/SLN
           CALL DSCAL(N,SCL,DT,1)
           SLN = STEPMX
        ENDIF

c compute SLOPE = G-Trans * DT

        SLOPE = DDOT(N,G,1,DT,1)

c initialization of RETCD

        RETCD = 0

c compute the smallest value allowable for the damping
c parameter ALMDA, i.e, ALMDAM

        RELENG = ZERO
        DO 20 I = 1,N
           TEMP1 = MAX(ABS(XC(I)), ONE)
           TEMP2 = ABS(DT(I))/TEMP1
           RELENG = MAX(RELENG, TEMP2)
 20     CONTINUE
        ALMDAM = STEPTL/RELENG
        ALMDA = ONE

c compute the next iterate XP

        DO 30 I = 1,N
           XP(I) = XC(I)+ALMDA*DT(I)
 30     CONTINUE

c evaluate the objective function at XP and its residual

        CALL TSFSCL(XP,DX,DF,M,N,FVEC,FP)

        FPNORM = HALF*DNRM2(M,FP,1)**2

c test whether the full tensor step produces enough decrease in the
c l2 norm of of the objective function

        IF (FPNORM.LT.(FCNORM + (ALPHA* ALMDA * SLOPE))) RETURN

c compute the standard direction

        CALL TSCPSS(SHAT,MAXM,MAXN,M,N,P,METHOD,GLOBAL,EPSM,FQ,
     +              WRK1,WRK2,WRK3,WRK4,AJA,ANLS,DN,FQQ,PIVOT,
     +              CURPOS,PBAR,ZERO1,IERR,RESNEW,FLAG)

c compute a step in the standard direction

        CALL TSLSCH(M,N,XC,DN,G,STEPTL,DX,DF,FVEC,
     +              MXTAKE,STEPMX,WRK1,WRK2,FCNORM,F1N,RETCD1)

c test whether the tensor direction is descent

        DTNORM = DNRM2(N,DT,1)
        GNORM = DNRM2(N,G,1)
        IF(M.GT.N) THEN
           BETA = TENTH
        ELSE
           BETA = ALPHA
        ENDIF
        TEMP1 = -BETA*DTNORM*GNORM

c compute a step in the tensor direction

        IF(SLOPE .LE. TEMP1) THEN
 50      CONTINUE
         ALMDAT = ((-ALMDA**2)*SLOPE)/(TWO*(FPNORM-FCNORM-ALMDA*SLOPE))
           TEMP = ALMDA/TEN
           ALMDA = MAX(TEMP, ALMDAT)
           IF(ALMDA.LT.ALMDAM) THEN
              IF(RETCD1. EQ. 1) THEN
                 RETCD = 1
                 GO TO 70
              ENDIF
           ENDIF
           DO 60 I = 1,N
              XP(I) = XC(I)+ALMDA*DT(I)
 60        CONTINUE
           CALL TSFSCL(XP,DX,DF,M,N,FVEC,FP)
           FPNORM = HALF*DNRM2(M,FP,1)**2
           IF (FPNORM .GT.(FCNORM + (ALPHA* ALMDA * SLOPE))) GO TO 50
           IF(ALMDA.EQ.TENTH .AND. SLN.GT.Z99*STEPMX) MXTAKE=.TRUE.
 70        CONTINUE

c select the next iterate that produces the lower function value

           IF(F1N .LT. FPNORM) THEN
              CALL DCOPY(N,WRK1,1,XP,1)
              CALL DCOPY(M,WRK2,1,FP,1)
              FPNORM  =  F1N
           ENDIF
        ELSE
           CALL DCOPY(N,WRK1,1,XP,1)
           CALL DCOPY(M,WRK2,1,FP,1)
           FPNORM = F1N
        ENDIF

        RETURN
        END

       FUNCTION TSMFDA(ANLS,ADT,AG,CONST1,CONST2,ALPHA,DLT,
     +                 NR,M,N,P,NWTAKE,IERR,VN,VNP)

       INTEGER NR,M,N,P,IERR
       DOUBLE PRECISION ALPHA,DLT,TSMFDA
       DOUBLE PRECISION ADT(N),AG(N),CONST1(P),CONST2(P),VN(M),VNP(M)
       DOUBLE PRECISION ANLS(NR,P)
       LOGICAL NWTAKE

C***********************************************************************
C THIS FUNCTION COMPUTES THE DERIVATIVE OF || F + J*D + 0.5*A*D**2 ||**2
C IN THE L2 NORM SENS, WHERE D = ALPHA*DT + SQRT(DLT**2-ALPHA**2).
C***********************************************************************
C
C
C   INPUT PARAMETERS
C   ----------------
C
C       ANLS   : TENSOR MATRIX
C       FQ     : FUNCTION VALUE AT CURRENT ITERATE MULTIPLIED BY
C                ORTHOGONAL MATRICES
C       ADT    : JACOBIAN MATRIX TIMES DT (SEE SUBROUTINE TS2DTR)
C       AG     : JACOBIAN MATRIX TIMES GBAR (SEE SUBROUTINE TS2DTR)
C       CONST1 : SHAT-TRANS TIMES DT (SEE SUBROUTINE TS2DTR)
C       CONST2 : SHAT-TRANS TIMES GBAR (SEE SUBROUTINE TS2DTR)
C       ALPHA  : POINT AT WHICH TO EVALUATE THE DERIVATIVE OF FUNCTION
C       DLT    : CURRENT TRUST RADIUS
C       NR     : LEADING DIMENSION OF THE JACOBIAN
C       M,N    : DIMENSIONS OF THE PROBLEM
C       P      : COLUMN DIMENSION OF THE MATRICES SHAT AND ANLS
C       NWTAKE : LOGICAL VARIABLE WITH THE FOLLOWING MEANINGS:
C                NWTAKE = .TRUE.  : STANDARD STEP TAKEN
C                NWTAKE = .FALSE. : TENSOR STEP TAKEN
C       IERR   : RETURN CODE FROM QRP FACTORIZATION ROUTINE:
C                IERR=0 : NO SINGULARITY OF JACOBIAN DETECTED
C                IERR=1 : SINGULARITY OF JACOBIAN DETECTED
C
C
C       OUTPUT PARAMETERS
C       -----------------
C
C
C       VN     : F + J*D + 0.5*A*D**2
C       VNP    : DERIVATIVE IN ALPHA OF F + J*D + 0.5*A*D**2
C       TSMFDA : DERIVATIVE IN ALPHA OF || F + J*D + 0.5*A*D**2 ||**2
C                WHERE D = ALPHA*DT + SQRT(DLT**2-ALPHA**2)
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DDOT
C       TENSOLVE      ...  TSMFDV
C
C***********************************************************************

       INTEGER LEN
       DOUBLE PRECISION DDOT

       CALL TSMFDV(ANLS,ADT,AG,CONST1,CONST2,ALPHA,DLT,
     +                   NR,M,N,P,NWTAKE,IERR,VNP)

       LEN = M
       IF(IERR.GT.0) LEN = M + N

       TSMFDA = DDOT(LEN,VNP,1,VN,1)

       RETURN
       END

       SUBROUTINE TSMFDV(ANLS,ADT,AG,CONST1,CONST2,ALPHA,DLT,
     +                   NR,M,N,P,NWTAKE,IERR,VNP)

       INTEGER NR,M,N,P,IERR
       DOUBLE PRECISION ALPHA,DLT
       DOUBLE PRECISION ADT(N),AG(N),CONST1(P)
       DOUBLE PRECISION CONST2(P),VNP(M),ANLS(NR,P)
       LOGICAL NWTAKE

C***********************************************************************
C THIS ROUTINE COMPUTES THE DERIVATIVE IN ALPHA OF THE VECTOR
C VN = F(XC) + J(XC)*D + 0.5*A*D**2, WHERE D = ALPHA*DT +
C SQRT(DLT**2-ALPHA**2).
C***********************************************************************
C
C
C      INPUT PARAMETERS :
C      -----------------
C
C      ANLS  : TENSOR TERM MATRIX
C       ADT  : JACOBIAN MATRIX TIMES DT (SEE SUBROUTINE TS2DTR)
C        AG  : JACOBIAN MATRIX TIMES GBAR (SEE SUBROUTINE TS2DTR)
C      CONST1: SHAT-TRANS TIMES DT (SEE SUBROUTINE TS2DTR)
C      CONST2: SHAT-TRANS TIMES GBAR (SEE SUBROUTINE TS2DTR)
C      ALPHA : POINT AT WHICH DERIVATIVE IS EVALUATED
C        DLT : CURRENT TRUST RADIUS
C        NR  : LEADING DIMENSION OF ANLS
C        M,N : DIMENSIONS OF THE PROBLEM
C        P   : COLUMN DIMENSION OF THE MATRIX ANLS
C        NWTAKE : LOGICAL VARIABLE WITH THE FOLLOWING MEANINGS :
C                NWTAKE = .TRUE.  : STANDARD STEP TAKEN
C                NWTAKE = .FALSE. : TENSOR STEP TAKEN
C        IERR : RETURN CODE FROM THE QRP FACTORIZATION ROUTINE
C               IERR = 0 : NO SINGULARITY OF JACOBIAN DETECTED
C               IERR = 1 : SINGULARITY OF JACOBIAN DETECTED
C
C
C       OUTPUT PARAMETERS :
C       -------------------
C
C       VNP  : THE DERIVATIVE IN ALPHA OF VN = F(XC) + J(XC)*D +
C              0.5*A*D**2, WHERE D = ALPHA*DT +  SQRT(DLT**2-ALPHA**2)
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DAXPY
C       TENSOLVE      ...  TSDLOD
C
C*******************************************************************

       INTEGER I,J,LEN
       DOUBLE PRECISION QUANT1,QUANT2,EXPR,CONST
       DOUBLE PRECISION ZERO,HALF,TWO
       INTRINSIC SQRT
       DATA ZERO,HALF,TWO/0.0D0,0.50D0,2.0D0/

       QUANT1 = SQRT(DLT**2 - ALPHA**2)
       EXPR = - ALPHA/QUANT1

       DO 10 I = 1,N
          VNP(I) = ADT(I) + EXPR*AG(I)
 10    CONTINUE

       CALL TSDLOD (M,ZERO,VNP(N+1),1)

       IF(NWTAKE) RETURN

       QUANT2 = QUANT1 - ALPHA**2/QUANT1

       LEN = M
       IF(IERR.GT.0) LEN = M + N

       DO 30 J = 1,P
           CONST = HALF*(TWO*ALPHA*(CONST1(J)**2 - CONST2(J)**2)
     +             +TWO*CONST1(J)*CONST2(J)*QUANT2)
           CALL DAXPY(LEN,CONST,ANLS(1,J),1,VNP,1)
 30    CONTINUE

       RETURN
       END

       SUBROUTINE TSMGSA(S,NR,N,SQRN,ITN,SHAT,P,IDP)

       INTEGER NR,N,SQRN,ITN,P
       INTEGER IDP(SQRN)
       DOUBLE PRECISION S(NR,SQRN),SHAT(NR,SQRN)

C*********************************************************************
C THIS ROUTINE FINDS A SET OF LINEARLY INDEPENDENT DIRECTIONS USING
C THE MODIFIED GRAM-SCHMIDT ALGORITHM.
C*********************************************************************
C
C       INPUT PARAMETERS :
C       ---------------
C
C       S   : MATRIX OF PAST DIRECTIONS
C       NR  : LEADING DIMENSION OF THE MATRICES S AND SHAT
C       N   : ROW DIMENSION OF MATRIX S AND SHAT
C       SQRN: MAXIMUM COLUMN DIMENSION OF SHAT
C       ITN : CURRENT ITERATION NUMBER
C
C       OUTPUT PARAMETERS :
C       -------------------
C
C       SHAT: MATRIX OF LINEARLY INDEPENDENT DIRECTIONS
C       P   : COLUMN DIMENSION OF THE MATRIX SHAT
C       IDP : VECTOR THAT KEEPS TRACK OF THE INDICES CORRESPONDING TO
C             THE LINEARLY INDEPENDENT DIRECTIONS IN THE MATRIX S
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DAXPY,DCOPY,DDOT,DNRM2
C
C*********************************************************************

       INTEGER J,K,L
       DOUBLE PRECISION TOL,TJ,SJ,SUM,RTJS,ONE,TWO
       DOUBLE PRECISION DNRM2,DDOT
       INTRINSIC SQRT
       DATA ONE,TWO/1.0D0,2.0D0/

       IF(SQRN.LT.ITN) THEN
          K = SQRN
       ELSE
          K = ITN-1
       ENDIF

       TOL = SQRT(TWO)/TWO

       DO 10 J = 1,K
         CALL DCOPY(N,S(1,J),1,SHAT(1,J),1)
 10    CONTINUE

       P = 0
       DO 30 J = 1,K
          TJ = DNRM2(N,SHAT(1,J),1)
          SJ = DNRM2(N,S(1,J),1)
          IF(TJ/SJ.GT.TOL) THEN
              P = P+1
              IDP(P) = J
              RTJS = ONE/TJ**2
              DO 20 L = J+1,K
                 SUM = -RTJS*DDOT(N,SHAT(1,L),1,SHAT(1,J),1)
                 CALL DAXPY(N,SUM,SHAT(1,J),1,SHAT(1,L),1)
 20           CONTINUE
           ENDIF
 30    CONTINUE

       DO 40 J = 1,P
          CALL DCOPY(N,S(1,IDP(J)),1,SHAT(1,J),1)
 40    CONTINUE

       RETURN
       END

       FUNCTION TSMSDA(ANLS,FQ,ADT,AG,CONST1,CONST2,ALPHA,
     +                 DLT,NR,M,N,P,NWTAKE,IERR,SKIP,VN,VNP,VNS)

       INTEGER NR,M,N,P,IERR
       DOUBLE PRECISION ALPHA,DLT,TSMSDA
       DOUBLE PRECISION ADT(N),AG(N),VN(M),VNP(M)
       DOUBLE PRECISION VNS(M),ANLS(NR,P),FQ(M)
       DOUBLE PRECISION CONST1(P),CONST2(P)
       LOGICAL NWTAKE

C***********************************************************************
C THIS FUNCTION COMPUTES THE SECOND DERIVATIVE OF || F + J*D +
C 0.5*A*D**2 ||**2 IN THE L2 NORM SENS, WHERE D = ALPHA*DT +
C SQRT(DLT**2-ALPHA**2).
C***********************************************************************
C
C
C   INPUT PARAMETERS
C   ----------------
C
C       ANLS   : TENSOR TERM MATRIX AT CURRENT ITERATE
C       FQ     : FUNCTION VALUE AT CURRENT ITERATE MULTIPLIED BY
C                ORTHOGONAL MATRICES
C       ADT    : JACOBIAN MATRIX TIMES DT (SEE SUBROUTINE TS2DTR)
C        AG    : JACOBIAN MATRIX TIMES GBAR (SEE SUBROUTINE TS2DTR)
C       CONST1 : SHAT-TRANS TIMES DT (SEE SUBROUTINE TS2DTR)
C       CONST2 : SHAT-TRANS TIMES GBAR (SEE SUBROUTINE TS2DTR)
C       ALPHA  : POINT AT WHICH TO EVALUATE THE SECOND DERIVATIVE OF
C                FUNCTION
C       DLT    : CURRENT TRUST RADIUS
C       NR     : LEADING DIMENSION OF THE JACOBIAN
C       M,N    : DIMENSIONS OF THE PROBLEM
C       P      : COLUMN DIMENSION OF THE MATRICES SHAT AND ANLS
C       NWTAKE : LOGICAL VARIABLE WITH THE FOLLOWING MEANINGS:
C                NWTAKE = .TRUE.  : STANDARD STEP TAKEN
C                NWTAKE = .FALSE. : TENSOR STEP TAKEN
C       IERR   : RETURN CODE FROM QRP FACTORIZATION ROUTINE
C                IERR = 0 : NO SINGULARITY OF JACOBIAN DETECTED
C                IERR = 1 : SINGULARITY OF JACOBIAN DETECTED
C
C
C       OUTPUT PARAMETERS
C       -----------------
C
C       VN     : F + J*D + 0.5*A*D**2
C       VNP    : DERIVATIVE IN ALPHA OF F + J*D + 0.5*A*D**2
C       VNS    : SECOND DERIVATIVE IN ALPHA OF F + J*D + 0.5*A*D**2
C       TSMSDA : SECOND DERIVATIVE IN ALPHA OF || F + J*D +
C                0.5*A*D**2 ||**2
C                WHERE D=ALPHA*DT + SQRT(DLT**2-ALPHA**2)
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DDOT
C       TENSOLVE      ...  TSMAFA,TSMFDV,TSMSDV
C
C***********************************************************************

       INTEGER LEN
       DOUBLE PRECISION DDOT
       LOGICAL SKIP

       IF(.NOT. SKIP) THEN
          CALL TSMAFA(ANLS,FQ,ADT,AG,CONST1,CONST2,ALPHA,DLT,
     +                NR,M,N,P,NWTAKE,IERR,VN)
          CALL TSMFDV(ANLS,ADT,AG,CONST1,CONST2,ALPHA,DLT,
     +                NR,M,N,P,NWTAKE,IERR,VNP)
       ENDIF

       CALL TSMSDV(ANLS,AG,CONST1,CONST2,ALPHA,DLT,
     +             NR,M,N,P,NWTAKE,IERR,VNS)

       LEN = M
       IF(IERR.GT.0) LEN = M + N

       TSMSDA = DDOT(LEN,VNP,1,VNP,1)+DDOT(M,VNS,1,VN,1)

       RETURN
       END

       SUBROUTINE TSMSDV(ANLS,AG,CONST1,CONST2,ALPHA,DLT,
     +                   NR,M,N,P,NWTAKE,IERR,VNS)

       INTEGER NR,M,N,P,IERR
       DOUBLE PRECISION ALPHA,DLT
       DOUBLE PRECISION AG(N),CONST1(P)
       DOUBLE PRECISION CONST2(P),VNS(M),ANLS(NR,P)
       LOGICAL NWTAKE

C***********************************************************************
C THIS ROUTINE COMPUTES THE SECOND DERIVATIVE IN ALPHA OF THE VECTOR
C VN = F(XC) + J(XC)*D + 0.5*A*D**2, WHERE D = ALPHA*DT +
C SQRT(DLT**2-ALPHA**2).
C***********************************************************************
C
C
C      INPUT PARAMETERS :
C      -----------------
C
C      ANLS  : TENSOR TERM MATRIX
C       ADT  : JACOBIAN MATRIX TIMES DT (SEE SUBROUTINE TS2DTR)
C        AG  : JACOBIAN MATRIX TIMES GBAR (SEE SUBROUTINE TS2DTR)
C      CONST1: SHAT-TRANS * DT (SEE SUBROUTINE TS2DTR)
C      CONST2: SHAT-TRABS * GBAR (SEE SUBROUTINE TS2DTR)
C      ALPHA : POINT AT WHICH DERIVATIVE IS EVALUATED
C        DLT : CURRENT TRUST RADIUS
C        NR  : LEADING DIMENSION OF ANLS
C        M,N : DIMENSIONS OF THE PROBLEM
C        P   : COLUMN DIMENSION OF THE MATRIX ANLS
C        NWTAKE : LOGICAL VARIABLE WITH THE FOLLOWING MEANINGS :
C                 NWTAKE = .TRUE.  : STANDARD STEP TAKEN
C                 NWTAKE = .FALSE. : TENSOR STEP TAKEN
C        IERR : RETURN CODE FROM THE QRP FACTORIZATION ROUTINE :
C               IERR = 0 : NO SINGULARITY OF JACOBIAN DETECTED
C               IERR = 1 : SINGULARITY OF JACOBIAN DETECTED
C
C       OUTPUT PARAMETERS :
C       -------------------
C
C       VNP  : THE SECOND DERIVATIVE IN ALPHA OF VN = F(XC) + J(XC)*D
C              + 0.5*A*D**2, WHERE D = ALPHA*DT +  SQRT(DLT**2-ALPHA**2)
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DAXPY
C       TENSOLVE      ...  TSDLOD
C
C*******************************************************************

       INTEGER I,J,LEN
       DOUBLE PRECISION QUANT1,EXPR,CONST,QUANT2
       DOUBLE PRECISION ZERO,HALF,ONEPF,TWO,THREE
       INTRINSIC SQRT
       DATA ZERO,HALF,ONEPF,TWO,THREE/0.0D0,0.50D0,1.50D0,2.0D0,3.0D0/

       QUANT1 = DLT**2 - ALPHA**2
       EXPR = -DLT**2 * SQRT(QUANT1) / QUANT1**2
       DO 10 I = 1,N
          VNS(I) =  EXPR*AG(I)
 10    CONTINUE

       CALL TSDLOD (M,ZERO,VNS(N+1),1)

       IF(NWTAKE) RETURN

       QUANT2 = -THREE*ALPHA/SQRT(QUANT1)-ALPHA**3/QUANT1**ONEPF

       LEN = M
       IF(IERR .GT. 0) LEN = M + N

       DO 30 J = 1,P
          CONST = HALF*(TWO*(CONST1(J)**2 - CONST2(J)**2)
     +            +TWO*CONST1(J)*CONST2(J)*QUANT2)
          CALL DAXPY(LEN,CONST,ANLS(1,J),1,VNS,1)
 30    CONTINUE

       RETURN
       END

       SUBROUTINE TSMSLV(AJA,S,ANLS,FC,P,MAXM,MAXN,SQRN,M,N,EPSM,
     +                   METHOD,GLOBAL,WRK1,WRK2,WRK3,WRK4,X,TYPXU,
     +                   XPLS,GPLS,A,WRK,CURPOS,PBAR,PIVOT,FQ,FQQ,
     +                   DN,DT,RESTNS,RESNEW,ITRMCD,FLAG,ZERO1,IERR)

       INTEGER MAXM,MAXN,M,N,P,GLOBAL,ZERO1,FLAG
       INTEGER ITRMCD,IERR,MSG,ITNLIM,IPR,METHOD,SQRN
       INTEGER PIVOT(N),PBAR(N),CURPOS(N)
       DOUBLE PRECISION EPSM,RESTNS,RESNEW
       DOUBLE PRECISION AJA(MAXM,N),S(MAXN,P),ANLS(MAXM,P),FQ(M),FQQ(M)
       DOUBLE PRECISION WRK1(M),WRK2(M),WRK3(M),WRK4(M),DN(N),DT(N)
       DOUBLE PRECISION FC(M),X(P),TYPXU(P),XPLS(P),GPLS(P),A(SQRN,P)
       DOUBLE PRECISION WRK(SQRN,P)

C**********************************************************************
C THIS ROUTINE FINDS THE TENSOR AND STANDARD STEPS.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       ---------------
C
C       S    : MATRIX OF PAST LINEARLY INDEPENDENT DIRECTIONS
C       P    : COLUMN DIMENSION OF MATRICES ANLS AND S
C       MAXM : LEADING DIMENSION OF AJA AND ANLS
C       MAXN : LEADING DIMENSION OF S
C       SQRN : LEADING DIMENSION OF MATRICES A AND WRK
C       M,N  : DIMENSIONS OF PROBLEM
C       EPSM : MACHINE PRECISION
C       AJA  : JACOBIAN AT CURRENT POINT XC
C       ANLS : TENSOR TERM MATRIX AT XC
C       FC   : FUCTION VALUE XC
C       X    : ESTIMATE TO A ROOT OF FCN (USED BY UNCMIN)
C       TYPXU: TYPICAL SIZE FOR EACH COMPONENT OF X (USED BY UNCMIN)
C       A    : WORKSPACE FOR HESSIAN (OR ESTIMATE) (USED BY UNCMIN)
C       WRK  : WORKSPACE (USED BY UNCMIN)
C       METHOD : METHOD TO USE
C                METHOD = 0 : STANDARD METHOD IS USED
C                METHOD = 1 : TENSOR METHOD IS USED
C       GLOBAL : GLOBAL STRATEGY USED
C       WRK1,WRK2,WRK3,WRK4,FQ,FQQ,WRK3 : WORKSPACE
C
C       OUTPUT PARAMETERS :
C       ------------------
C
C       DN : STANDARD STEP
C       DT : TENSOR STEP
C       FLAG : RETURNED CODE WITH THE FOLLOWING MEANING :
C              FLAG = 0 : NO SINGULARITY DETECTED WHEN FACTORIZING AJA
C              FLAG = 1 : SINGULARITY DETECTED WHEN FACTORIZING AJA
C                         FROM 1 TO N-P
C              FLAG = 2 : SINGULARITY DETECTED WHEN FACTORIZING AJA
C                         FROM N-P TO N
C       IERR   : RETURNED CODE WITH THE FOLLOWING MEANING :
C                IERR = 0 : NO SINGULARITY DETECTED WHEN FACTORIZING AJA
C                IERR = 1 : SINGULARITY DETECTED WHEN FACTORIZING AJA
C       XPLS : LOCAL MINIMUM OF OPTIMIZATION FUNCTION FCN (USED
C              BY UNCMIN)
C       FPLS : FUNCTION VALUE AT SOLUTION OF OPTIMIZATION FUNCTION FCN
C              (USED IN UNCMIN)
C       GPLS : GRADIENT AT SOLUTION XPLS (USED BY UNCMIN)
C       CURPOS,PIVOT,PBAR : PIVOT VECTORS
C       RESTNS : TENSOR RESIDUAL
C       RESNEW : STANDARD RESIDUAL
C       ITRMCD : TERMINATION CODE (FOR UNCMIN)
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DCOPY
C       TENSOLVE      ...  TSQLFC,QTRNS,TSQRFC,TSQMTS,TSQMUV,TSSQP1
C       TENSOLVE      ...  TSDLOD,TSQ1P1,TSD1SV,TSPRMV,TSQMLV,TSCPSS
C       UNCMIN        ...  DFAUT,OPTIF9
C
C**********************************************************************

      INTEGER Q,METH,IEXP,NDIGIT,IAGFLG,IAHFLG
      DOUBLE PRECISION ROOT,TYPFU,DLT,GRADLT,STEPMX,STEPTL,FPLS
      DOUBLE PRECISION ZERO,ONE,TWO
      INTRINSIC SQRT
      EXTERNAL TSQFCN,TSDFCN,D2FCN
      DATA ZERO,ONE,TWO/0.0D0,1.0D0,2.0D0/

      ITRMCD = 0
      IF(N .EQ. 1) THEN
         S(2,1) = ONE
         S(3,1) = ONE
         CURPOS(1) = 1
         CALL DCOPY(M,FC,1,FQ,1)
      ELSE

c perform a QL decomposition of S

         CALL TSQLFC(S,MAXN,N,P,IERR)

c compute AJA times Q-trans

         CALL TSJQTP(S,MAXM,MAXN,N,M,P,WRK1,FQ,AJA)

c perform a QR factorization of AJA

         CALL TSQRFC(AJA,MAXM,N,M,1,N-P,IERR,EPSM,WRK1,CURPOS,ZERO1)

         IF(IERR.EQ.1) THEN
            Q = N-ZERO1+1
         ELSE
            Q = P
         ENDIF
         CALL TSQMTS(ANLS,AJA,MAXM,M,N,M,P,1,WRK1,ZERO1)

         CALL TSQMUV(AJA,FC,FQ,MAXM,M,1,ZERO1,.FALSE.)
      ENDIF

c minimize the lower m-n+q quadratic equations in p unknowns
c of the tensor model. The minimization is performed analytically
c if p=1,q>1, or p=1,q=1,m>n, or n=1,m>n. Otherwise an unconstrained
c minimization software, UNCMIN, is used.

      IF((P.EQ.1.AND.Q.GT.1).OR.(P.EQ.1 .AND. Q.EQ.1 .AND. M.GT.N)
     +   .OR. (N .EQ. 1 .AND. M .GT. N)) THEN
         CALL TSSQP1(AJA,ANLS,S,FQ,MAXM,MAXN,M,N,Q,ROOT,RESTNS)
         XPLS(1) = ROOT
      ELSEIF((M.EQ.N).AND.(P.EQ.1).AND.(Q.EQ.1) .OR.
     +      (M.EQ.1.AND.N.EQ.1)) THEN
         CALL TSQ1P1(AJA,ANLS,S,FQ,MAXM,MAXN,N,ROOT,RESTNS)
         XPLS(1) = ROOT
      ELSE
         CALL DFAUT(P,TYPXU,TYPFU,METH,IEXP,MSG,NDIGIT,ITNLIM,
     +              IAGFLG,IAHFLG,IPR,DLT,GRADLT,STEPMX,STEPTL)

         IAGFLG = 1
         IAHFLG = 0
         IEXP = 0
         METH = 2
         MSG = 9

         CALL TSDLOD (P,ZERO,X,1)

         CALL OPTIF9(SQRN,P,X,TSQFCN,TSDFCN,D2FCN,TYPXU,TYPFU,METH,IEXP,
     +               MSG,NDIGIT,ITNLIM,IAGFLG,IAHFLG,IPR,DLT,GRADLT,
     +               STEPMX,STEPTL,XPLS,FPLS,GPLS,ITRMCD,A,WRK,AJA,
     +               ANLS,S,FQ,WRK1,FQQ,WRK2,WRK3,WRK4,MAXM,MAXN,M,N,Q)

c compute the tensor residual

         RESTNS = SQRT(TWO*FPLS)
       ENDIF

       CALL DCOPY(P,XPLS,1,WRK4(N-P+1),1)

       IF(N .EQ. 1) THEN
          DT(1) = WRK4(1)
       ELSE

c compute the first n-p components of the tensor step

          CALL TSD1SV(AJA,S,ANLS,FQ,XPLS,MAXM,MAXN,M,N,P,Q,EPSM,
     +                WRK1,FQQ,WRK2,PIVOT,WRK3)
          CALL TSPRMV(WRK4,WRK3,CURPOS,N-P,0)

c premultiply the tensor step by the orthogonal matrix resulting
c from the QL factorization of S

          CALL TSQMLV(MAXN,N,P,S,WRK4,DT,.TRUE.)
       ENDIF

c compute the standard step if needed

       IF(GLOBAL .EQ. 1 .OR. (M .GT. N .AND. GLOBAL .EQ. 0)) THEN
          CALL TSCPSS(S,MAXM,MAXN,M,N,P,METHOD,GLOBAL,EPSM,FQ,
     +                WRK1,WRK2,WRK3,WRK4,AJA,ANLS,DN,FQQ,PIVOT,
     +                CURPOS,PBAR,ZERO1,IERR,RESNEW,FLAG)
       ENDIF

       RETURN
       END

       SUBROUTINE TSNECI(MAXM,MAXN,MAXP,X0,M,N,TYPX,TYPF,ITNLIM,
     +                   JACFLG,GRADTL,STEPTL,FTOL,METHOD,GLOBAL,
     +                   STEPMX,DLT,IPR,WRKUNC,LUNC,WRKNEM,LNEM,
     +                   WRKNEN,LNEN,IWRKN,LIN,FVEC,JAC,MSG,
     +                   XP,FP,GP,TERMCD)

       INTEGER MAXM,MAXN,M,N,MAXP,JACFLG,ITNLIM,TERMCD,METHOD
       INTEGER MSG,GLOBAL,IPR,LUNC,LNEM,LNEN,LIN
       INTEGER IWRKN(MAXN,LIN)
       DOUBLE PRECISION STEPTL,GRADTL,FTOL,STEPMX,DLT
       DOUBLE PRECISION XP(N),FP(M),GP(N),X0(N)
       DOUBLE PRECISION WRKUNC(MAXP,LUNC)
       DOUBLE PRECISION WRKNEM(MAXM,LNEM)
       DOUBLE PRECISION WRKNEN(MAXN,LNEN)
       DOUBLE PRECISION TYPX(N),TYPF(M)
       EXTERNAL FVEC,JAC

C
C**********************************************************************
C THIS ROUTINE PROVIDES A COMPLETE INTERFACE TO THE NONLINEAR EQUATION/
C NONLINEAR LEAST SQUARES PACKAGE. THE USER HAS FULL CONTROL OVER
C THE OPTIONS.
C**********************************************************************
C
C       SUBPROGRAMS CALLED:
C
C       TENSOLVE      ...  TSCHKI,TSNESV
C
C**********************************************************************

       INTEGER SQRN
       DOUBLE PRECISION EPSM

c check input parameters

       CALL TSCHKI(MAXM,MAXN,MAXP,M,N,LUNC,LNEM,LNEN,LIN,GRADTL,STEPTL,
     +            FTOL,ITNLIM,JACFLG,METHOD,GLOBAL,STEPMX,DLT,EPSM,
     +            MSG,TYPX,TYPF,WRKNEN(1,2),WRKNEM(1,2),SQRN,
     +            TERMCD,IPR)
       IF(MSG.LT.0) RETURN

c call nonlinear equations/nonlinear least squares solver

       CALL TSNESV(MAXM,MAXN,MAXP,X0,M,N,TYPX,TYPF,ITNLIM,JACFLG,
     +             GRADTL,STEPTL,FTOL,METHOD,GLOBAL,STEPMX,DLT,IPR,
     +             WRKUNC(1,1),WRKUNC(1,2),WRKUNC(1,3),WRKUNC(1,4),
     +             WRKUNC(1,5),WRKUNC(1,SQRN+5),WRKNEM(1,2),WRKNEM(1,3),
     +             WRKNEM(1,4),WRKNEM(1,5),WRKNEM(1,6),WRKNEM(1,7),
     +             WRKNEM(1,8),WRKNEM(1,9),WRKNEM(1,10),WRKNEM(1,11),
     +             WRKNEM(1,12),WRKNEM(1,SQRN+12),WRKNEM(1,2*SQRN+12),
     +             WRKNEN(1,2),WRKNEN(1,3),WRKNEN(1,4),WRKNEN(1,5),
     +             WRKNEN(1,6),WRKNEN(1,7),WRKNEN(1,8),WRKNEN(1,9),
     +             WRKNEN(1,10),WRKNEN(1,SQRN+10),IWRKN(1,1),IWRKN(1,2),
     +             IWRKN(1,3),EPSM,SQRN,FVEC,JAC,MSG,XP,FP,GP,
     +             TERMCD)

       RETURN
       END

       SUBROUTINE TSNESI(MAXM,MAXN,MAXP,X0,M,N,WRKUNC,LUNC,WRKNEM,
     +                   LNEM,WRKNEN,LNEN,IWRKN,LIN,FVEC,MSG,XP,
     +                   FP,GP,TERMCD)

       INTEGER MAXM,MAXN,M,N,MAXP,JACFLG,ITNLIM,TERMCD,METHOD
       INTEGER GLOBAL,MSG,IPR,LUNC,LNEM,LNEN,LIN
       INTEGER IWRKN(MAXN,LIN)
       DOUBLE PRECISION STEPTL,GRADTL,FTOL,STEPMX,DLT
       DOUBLE PRECISION XP(N),FP(M),GP(N),X0(N)
       DOUBLE PRECISION WRKUNC(MAXP,LUNC)
       DOUBLE PRECISION WRKNEM(MAXM,LNEM)
       DOUBLE PRECISION WRKNEN(MAXN,LNEN)
       EXTERNAL TSDUMJ,FVEC

C**********************************************************************
C THIS ROUTINE PROVIDES A SIMPLE INTERFACE TO THE NONLINEAR EQUATION/
C NONLINEAR LEAST SQUARES PROBLEMS PACKAGE.  THE USER HAS NO CONTROL
C OVER THE PACKAGE OPTIONS.
C**********************************************************************
C
C       SUBPROGRAMS CALLED:
C
C       TENSOLVE      ...  TSDFLT,TSCHKI,TSNESV
C
C**********************************************************************

       INTEGER SQRN
       DOUBLE PRECISION EPSM

c set default values for each variable to the nonlinear equations/
c nonlinear least squares solver

       CALL TSDFLT(M,N,ITNLIM,JACFLG,GRADTL,STEPTL,FTOL,METHOD,GLOBAL,
     +             STEPMX,DLT,WRKNEN(1,1),WRKNEM(1,1),IPR,MSG)

c check input parameters

       CALL TSCHKI(MAXM,MAXN,MAXP,M,N,LUNC,LNEM,LNEN,LIN,GRADTL,STEPTL,
     +             FTOL,ITNLIM,JACFLG,METHOD,GLOBAL,STEPMX,DLT,EPSM,
     +             MSG,WRKNEN(1,1),WRKNEM(1,1),WRKNEN(1,2),WRKNEM(1,2),
     +             SQRN,TERMCD,IPR)
       IF(MSG.LT.0) RETURN

c call nonlinear equations/nonlinear least squares solver

       CALL TSNESV(MAXM,MAXN,MAXP,X0,M,N,WRKNEN(1,1),WRKNEM(1,1),ITNLIM,
     +            JACFLG,GRADTL,STEPTL,FTOL,METHOD,GLOBAL,STEPMX,DLT,
     +            IPR,WRKUNC(1,1),WRKUNC(1,2),WRKUNC(1,3),WRKUNC(1,4),
     +            WRKUNC(1,5),WRKUNC(1,SQRN+5), WRKNEM(1,2),WRKNEM(1,3),
     +            WRKNEM(1,4),WRKNEM(1,5),WRKNEM(1,6),WRKNEM(1,7),
     +            WRKNEM(1,8),WRKNEM(1,9),WRKNEM(1,10),WRKNEM(1,11),
     +            WRKNEM(1,12),WRKNEM(1,SQRN+12),WRKNEM(1,2*SQRN+12),
     +            WRKNEN(1,2),WRKNEN(1,3),WRKNEN(1,4),WRKNEN(1,5),
     +            WRKNEN(1,6),WRKNEN(1,7),WRKNEN(1,8),WRKNEN(1,9),
     +            WRKNEN(1,10),WRKNEN(1,SQRN+10),IWRKN(1,1),IWRKN(1,2),
     +            IWRKN(1,3),EPSM,SQRN,FVEC,TSDUMJ,MSG,XP,FP,GP,
     +            TERMCD)

       RETURN
       END

       SUBROUTINE TSNESV(MAXM,MAXN,MAXP,XC,M,N,TYPX,TYPF,ITNLIM,
     +                  JACFLG,GRADTL,STEPTL,FTOL,METHOD,GLOBAL,
     +                  STEPMX,DLT,IPR,X,TYPXU,XPLS,GPLS,A,WRK,DFN,
     +                  WRK1,WRK2,WRK3,WRK4,WRK5,FQ,FQQ,FC,FHAT,
     +                  ANLS,FV,AJA,DXN,DN,DT,DF,D,GBAR,DBAR,DBARP,
     +                  S,SHAT,CURPOS,PIVOT,PBAR,EPSM,SQRN,FVEC,
     +                  JAC,MSG,XP,FP,GP,TERMCD)

        INTEGER MAXM,MAXN,MAXP,M,N,SQRN,TERMCD
        INTEGER ITNLIM,JACFLG,METHOD,GLOBAL,MSG,IPR
        INTEGER PBAR(N),CURPOS(N),PIVOT(N)
        DOUBLE PRECISION GRADTL,STEPTL,FTOL,STEPMX,DLT,FPLS,EPSM
        DOUBLE PRECISION TYPXU(SQRN),XPLS(SQRN),GPLS(SQRN),A(MAXP,SQRN)
        DOUBLE PRECISION WRK(MAXP,SQRN),X(SQRN),AJA(MAXM,N),S(MAXN,SQRN)
        DOUBLE PRECISION ANLS(MAXM,SQRN),SHAT(MAXN,SQRN),FV(MAXM,SQRN)
        DOUBLE PRECISION XC(N),FC(M),XP(N),FP(M),DN(N),DT(N),DF(N)
        DOUBLE PRECISION D(N),WRK1(M),WRK2(M),WRK3(M),WRK4(M)
        DOUBLE PRECISION WRK5(M),FQQ(M),FQ(M),GP(N),FHAT(M)
        DOUBLE PRECISION GBAR(N),DBAR(N),DBARP(N)
        DOUBLE PRECISION TYPX(N),TYPF(M),DXN(N),DFN(M)
        EXTERNAL FVEC,JAC

C**********************************************************************
C THIS IS THE DRIVER FOR NONLINEAR EQUATIONS/NONLINEAR LEAST SQUARES
C PROBLEMS.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       MAXM   : LEADING DIMENSION OF AJA, ANLS, AND FV
C       MAXN   : LEADING DIMENSION OF S AND SHAT
C       XC     : INITIAL ESTIMATE OF SOLUTION
C       M,N    : DIMENSIONS OF PROBLEM
C       TYPX   : TYPICAL SIZE FOR EACH COMPONENT OF X
C       TYPF   : TYPICAL SIZE FOR EACH COMPONENT OF F
C       ITNLIM : MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
C       JACFLG : JACOBIAN FLAG WITH THE FOLLOWING MEANINGS:
C                JACFLG = 1 IF ANALYTIC JACOBIAN SUPPLIED
C                JACFLG = 0 IF ANALYTIC JACOBIAN NOT SUPPLIED
C       GRADTL : TOLERANCE AT WHICH GRADIENT IS CONSIDERED CLOSE ENOUGH
C                TO ZERO TO TERMINATE ALGORITHM
C       STEPTL : TOLERANCE AT WHICH SUCCESSIVE ITERATES ARE CONSIDERED
C                CLOSE ENOUGH TO TERMINATE ALGORITHM
C       FTOL : TOLERANCE AT WHICH FUNCTION VALUE IS CONSIDERED CLOSE
C                ENOUGH TO ZERO
C       METHOD : METHOD TO USE
C                METHOD = 0 : STANDARD METHOD IS USED
C                METHOD = 1 : TENSOR METHOD IS USED
C       GLOBAL : GLOBAL STRATEGY TO USE
C                GLOBAL = 0 : LINE SEARCH
C                GLOBAL = 1 : 2-DIMENSIONAL TRUST REGION
C       STEPMX : MAXIMUM ALLOWABLE STEP SIZE
C       DLT    : TRUST REGION RADIUS
C       IPR    : DEVICE TO WHICH TO SEND OUTPUT
C       X      : ESTIMATE TO A ROOT OF FCN ( USED BY UNCMIN)
C       TYPXU  : TYPICAL SIZE FOR EACH COMPONENT OF X (USED BY UNCMIN)
C       XPLS   : LOCAL MINIMUM OF OPTIMIZATION FUNCTION FCN USED BY
C                UNCMIN
C       GPLS   : GRADIENT AT SOLUTION XPLS (USED BY UNCMIN)
C       A      : WORKSPACE FOR HESSIAN (OR ESTIMATE) (USED BY UNCMIN)
C       WRK    : WORKSPACE (USED BY UNCMIN)
C       WRK1,WRK2,WRK3,WRK4,WRK5,FQ,FQQ:  WORKSPACE
C       FC     : FUNCTION VALUE AT CURRENT ITERATE
C       FHAT   : WORKSPACE
C       DFN    : DIAGONAL SCALING MATRIX FOR F
C       ANLS   : TENSOR TERM MATRIX
C       FV     : WORKSPACE USED TO STORE PAST FUNCTION VALUES
C       AJA    : JACOBIAN MATRIX
C       DN     : STANDARD STEP
C       DT     : TENSOR STEP
C       DF,D,GBAR,DBAR,DBARP : WORKSPACE
C       DXN    : DIAGONAL SCALING MATRIX FOR X
C       S      : MATRIX OF PREVIOUS DIRECTIONS
C       SHAT   : MATRIX OF PAST LINEARLY INDEPENDENT DIRECTIONS
C       CURPOS,PIVOT,PBAR : PIVOT VECTORS
C       SQRN   : MAXIMUM COLUMN DIMENSION OF ANLS, S, AND SHAT
C       EPSM   : MACHINE PRECISION
C       FVEC   : NAME OF SUBROUTINE TO EVALUATE FUNCTION
C       JAC    : (OPTIONAL) NAME OF SUBROUTINE TO EVALUATE JACOBIAN.
C                MUST BE DECLARED EXTERNAL IN CALLING ROUTINE
C
C
C       INPUT-OUTPUT PARAMETERS :
C       ------------------------
C
C       MSG : MESSAGE TO INHIBIT CERTAIN AUTOMATIC CHECKS + OUTPUT
C
C
C       OUTPUT PARAMETERS :
C       -----------------
C
C       XP : SOLUTION TO THE SYSTEM OF NONLINEAR EQUATIONS
C       FP : FUNCTION VALUE AT THE SOLUTION
C       GP : GRADIENT AT THE SOLUTION
C       TERMCD : TERMINATION CODE
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DCOPY,DNRM2
C       LEVEL 2 BLAS  ...  DGEMV
C       TENSOLVE      ...  TSSCLX,TSFSCL,TSSCLJ,TSCHKJ,TSNSTP,TSSSTP,
C       TENSOLVE      ...  TSLSCH,TS2DTR,TSRSLT,TSMGSA,TSFRMT,TSMSLV,
C       TENSOLVE      ...  TSSLCT,TSMDLS,TSUPSF
C
C*********************************************************************

        INTEGER P,ITN,I,J,FLAG,RETCD,ZERO1,IERR,ITRMCD,ICSCMX
        DOUBLE PRECISION FNORM,RESTNS,RESNEW
        DOUBLE PRECISION ZERO,HALF,ONE
        DOUBLE PRECISION DNRM2
        LOGICAL NWTAKE,MXTAKE
        DATA ZERO,HALF,ONE/0.0D0,0.50D0,1.0D0/

c-----------------
c initialization
c-----------------

        ITN = 0
        NWTAKE = .TRUE.

        CALL TSSCLX(XC,DXN,N)

        IF(MOD(MSG/8,2).NE.1) THEN
           WRITE(IPR,896)
           WRITE(IPR,900) (TYPX(I),I = 1,N)
           WRITE(IPR,897)
           WRITE(IPR,900) (DXN(I),I = 1,N)
           WRITE(IPR,898)
           WRITE(IPR,900) (TYPF(I),I = 1,M)
           WRITE(IPR,899)
           WRITE(IPR,900) (DFN(I),I = 1,M)
           WRITE(IPR,901) JACFLG
           WRITE(IPR,902) METHOD
           WRITE(IPR,903) GLOBAL
           WRITE(IPR,904) ITNLIM
           WRITE(IPR,905) EPSM
           WRITE(IPR,906) STEPTL
           WRITE(IPR,907) GRADTL
           WRITE(IPR,908) FTOL
           WRITE(IPR,909) STEPMX
           WRITE(IPR,910) DLT
        ENDIF

c evaluate analytic or finite difference Jacobian and check analytic
c Jacobian, if requested

        CALL TSFSCL(XC,DXN,DFN,M,N,FVEC,FC)
        CALL TSSCLJ(XC,DXN,TYPX,FC,DFN,FHAT,MAXM,M,N,
     +              EPSM,JACFLG,FVEC,JAC,AJA)
        IF(JACFLG.EQ.1) THEN
          IF(MOD(MSG/2,2).EQ.0) THEN
             CALL TSCHKJ(AJA,XC,FC,MAXM,M,N,EPSM,DFN,DXN,TYPX,
     +                   IPR,FHAT,WRK1,FVEC,MSG)
             IF(MSG.LT.0) RETURN
          ENDIF
        ENDIF

c compute the gradient at the current iterate XC

        CALL DGEMV('T',M,N,ONE,AJA,MAXM,FC,1,ZERO,GP,1)

c compute the residual of FC

        FNORM = HALF*DNRM2(M,FC,1)**2

c check stopping criteria for input XC

        CALL TSNSTP(GP,XC,FC,XC,STEPTL,GRADTL,RETCD,FTOL,ITN,
     +              ITNLIM,ICSCMX,MXTAKE,M,N,MSG,IPR,FNORM,TERMCD)

        IF(TERMCD.GT.0) THEN
           FPLS = FNORM
           GO TO 120
        ENDIF

c---------------
c iteration 1
c---------------

        ITN = 1

c compute the standard step

        CALL DCOPY(M,FC,1,FHAT,1)

        CALL TSSSTP(AJA,FHAT,M,N,MAXM,EPSM,GLOBAL,WRK1,WRK2,WRK3,
     +              DN,FQQ,PIVOT,PBAR,IERR)

c choose next iterate XP by a global strategy

        IF(GLOBAL.EQ.0) THEN
          CALL TSLSCH(M,N,XC,DN,GP,STEPTL,DXN,DFN,FVEC,
     +                MXTAKE,STEPMX,XP,FP,FNORM,FPLS,RETCD)
        ELSE
          DO 20 I = 1,N
            DO 10 J = 1,SQRN
              SHAT(I,J) = ZERO
   10       CONTINUE
   20     CONTINUE
          CALL TS2DTR(AJA,SHAT,ANLS,DN,GP,GBAR,XC,METHOD,
     +                NWTAKE,STEPMX,STEPTL,EPSM,MXTAKE,DLT,
     +                FQQ,MAXM,MAXN,M,N,SQRN,CURPOS,PIVOT,
     +                PBAR,ITN,IERR,FLAG,DXN,DFN,FVEC,DBAR,
     +                DBARP,D,FHAT,WRK1,WRK2,WRK3,WRK4,WRK5,
     +                XPLS,GPLS,FNORM,XP,FP,FPLS,RETCD)
        ENDIF

        IF(MOD(MSG/8,2).EQ.0) CALL TSRSLT(N,XC,FNORM,GP,0,IPR)

c evaluate the Jacobian at the new iterate XP

        CALL TSSCLJ(XP,DXN,TYPX,FP,DFN,FHAT,MAXM,M,N,EPSM,JACFLG,
     +              FVEC,JAC,AJA)

c compute the gradient at the new iterate XP

        CALL DGEMV('T',M,N,ONE,AJA,MAXM,FP,1,ZERO,GP,1)

c check stopping criteria for the new iterate XP

        CALL TSNSTP(GP,XP,FP,XC,STEPTL,GRADTL,RETCD,FTOL,ITN,
     +              ITNLIM,ICSCMX,MXTAKE,M,N,MSG,IPR,FPLS,TERMCD)

        IF(TERMCD.GT.0) GO TO 120
        IF(MOD(MSG/16,2).EQ.1) CALL TSRSLT(N,XP,FPLS,GP,ITN,IPR)

c update S and FV

        DO 40 I = 1,N
           S(I,1) = XC(I)-XP(I)
 40     CONTINUE
        CALL DCOPY(M,FC,1,FV(1,1),1)

c update XC and FC

        CALL DCOPY(N,XP,1,XC,1)
        CALL DCOPY(M,FP,1,FC,1)
        FNORM = FPLS

c---------------
c iteration > 1
c---------------

 80    ITN = ITN+1

c if the standard method is selected then compute the standard step

              IF(METHOD.EQ.0) THEN
                 CALL DCOPY(M,FC,1,FHAT,1)
                 CALL TSSSTP(AJA,FHAT,M,N,MAXM,EPSM,GLOBAL,WRK1,WRK2,
     +                       WRK3,DF,FQQ,PIVOT,PBAR,IERR)
              ENDIF

c if the tensor method is selected then form the tensor model

              IF(METHOD.EQ.1) THEN

c select the past linearly independent directions

                 CALL TSMGSA(S,MAXN,N,SQRN,ITN,SHAT,P,CURPOS)

c form the tensor term

                 CALL TSFRMT(SHAT,S,AJA,FV,FC,MAXM,MAXN,MAXP,M,N,P,
     +                       CURPOS,A,X,XPLS,GPLS,ANLS)

c solve the tensor model for the tensor step DT and compute DN
c as a by-product if the global strategy selected is the
c two-dimensional trust region or M > N

                 CALL TSMSLV(AJA,SHAT,ANLS,FC,P,MAXM,MAXN,SQRN,M,N,
     +                       EPSM,METHOD,GLOBAL,WRK1,WRK2,WRK3,WRK4,
     +                       X,TYPXU,XPLS,GPLS,A,WRK,CURPOS,PBAR,PIVOT,
     +                       FQ,FQQ,DN,DT,RESTNS,RESNEW,ITRMCD,FLAG,
     +                       ZERO1,IERR)

c decide which step to use (DN or DT)

                 IF(GLOBAL.EQ.1 .OR. (M.GT.N .AND. GLOBAL.EQ.0)) THEN
                    CALL TSSLCT(RESTNS,RESNEW,ITRMCD,FC,M,N,DN,DT,GP,
     +                          DF,NWTAKE)
                 ENDIF

              ENDIF

c choose the next iterate XP by a global strategy

             IF(GLOBAL.EQ.0) THEN
                IF(METHOD.EQ.0) THEN
                  CALL TSLSCH(M,N,XC,DF,GP,STEPTL,DXN,DFN,FVEC,
     +                        MXTAKE,STEPMX,XP,FP,FNORM,FPLS,RETCD)
                 ELSEIF(M.EQ.N) THEN
                CALL TSMDLS(AJA,SHAT,ANLS,XC,M,N,MAXM,MAXN,P,DT,GP,
     +                      DXN,DFN,FVEC,METHOD,STEPTL,GLOBAL,STEPMX,
     +                      EPSM,FQ,WRK1,WRK2,WRK3,WRK4,DN,FQQ,PIVOT,
     +                      CURPOS,PBAR,MXTAKE,XP,FP,FNORM,FPLS,
     +                      ZERO1,RETCD,IERR)
                 ELSE
                  CALL TSLSCH(M,N,XC,DF,GP,STEPTL,DXN,DFN,FVEC,
     +                        MXTAKE,STEPMX,XP,FP,FNORM,FPLS,RETCD)
                 ENDIF
              ELSE
                  CALL TS2DTR(AJA,SHAT,ANLS,DF,GP,GBAR,XC,
     +                 METHOD,NWTAKE,STEPMX,STEPTL,EPSM,MXTAKE,
     +                 DLT,FQQ,MAXM,MAXN,M,N,P,CURPOS,PIVOT,
     +                 PBAR,ITN,IERR,FLAG,DXN,DFN,FVEC,DBAR,
     +                 DBARP,D,FHAT,WRK1,WRK2,WRK3,WRK4,WRK5,
     +                 XPLS,GPLS,FNORM,XP,FP,FPLS,RETCD)
              ENDIF

c evaluate the Jacobian at the new iterate XP

              CALL TSSCLJ(XP,DXN,TYPX,FP,DFN,FHAT,MAXM,M,N,EPSM,
     +                    JACFLG,FVEC,JAC,AJA)

c evaluate the gradient at the new iterate XP

              CALL DGEMV('T',M,N,ONE,AJA,MAXM,FP,1,ZERO,GP,1)

c check stopping criteria for the new iterate XP

              CALL TSNSTP(GP,XP,FP,XC,STEPTL,GRADTL,RETCD,FTOL,ITN,
     +                    ITNLIM,ICSCMX,MXTAKE,M,N,MSG,IPR,FPLS,TERMCD)

              IF(TERMCD.GT.0) GO TO 120
              IF(MOD(MSG/16,2).EQ.1) CALL TSRSLT(N,XP,FPLS,GP,ITN,IPR)

c if tensor method is selected then update the matrices S and FV

              IF(METHOD.EQ.1) THEN
                 CALL TSUPSF(FC,XC,XP,SQRN,ITN,MAXM,MAXN,M,N,WRK1,S,FV)
              ENDIF

c update XC, FC, and FNORM

               CALL DCOPY(N,XP,1,XC,1)
               CALL DCOPY(M,FP,1,FC,1)
               FNORM = FPLS
        GO TO 80

c termination

 120    IF(MOD(MSG/8,2).EQ.0) THEN
           IF(ITN.NE.0) THEN
              CALL TSRSLT(N,XP,FPLS,GP,ITN,IPR)
           ELSE
              FPLS = HALF*DNRM2(M,FC,1)**2
              CALL TSRSLT(N,XC,FPLS,GP,ITN,IPR)
           ENDIF
        ENDIF

 896    FORMAT('  TSNESV      TYPICAL X')
 897    FORMAT('  TSNESV      DIAGONAL SCALING MATRIX FOR X')
 898    FORMAT('  TSNESV      TYPICAL F')
 899    FORMAT('  TSNESV      DIAGONAL SCALING MATRIX FOR F')
 900    FORMAT(100('  TSNESV     ',3(D20.13,3X),/))
 901    FORMAT('  TSNESV      JACOBIAN FLAG      = ',I1)
 902    FORMAT('  TSNESV      METHOD USED        = ',I1)
 903    FORMAT('  TSNESV      GLOBAL STRATEGY    = ',I1)
 904    FORMAT('  TSNESV      ITERATION LIMIT    =',I5)
 905    FORMAT('  TSNESV      MACHINE EPSILON    =',D20.13)
 906    FORMAT('  TSNESV      STEP TOLERANCE     =',D20.13)
 907    FORMAT('  TSNESV      GRADIENT TOLERANCE =',D20.13)
 908    FORMAT('  TSNESV      FUNCTION TOLERANCE =',D20.13)
 909    FORMAT('  TSNESV      MAXIMUM STEP SIZE  =',D20.13)
 910    FORMAT('  TSNESV      TRUST REG RADIUS   =',D20.13)
        END

       SUBROUTINE TSNSTP(G,XPLUS,FPLUS,XC,STEPTL,GRADTL,RETCD,FTOL,ITN,
     +                  ITNLIM,ICSCMX,MXTAKE,M,N,MSG,IPR,FNORM,TERMCD)

       INTEGER M,N,ITN,ITNLIM,MSG,IPR,TERMCD,RETCD,ICSCMX
       DOUBLE PRECISION STEPTL,GRADTL,FTOL,FNORM
       DOUBLE PRECISION XPLUS(N),FPLUS(M),G(N),XC(N)
       LOGICAL MXTAKE

C**********************************************************************
C THIS ROUTINE DECIDES WHETHER TO TERMINATE THE NONLINEAR ALGORITHM.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       ------------------
C
C       G     : GRADIENT AT XC
C       XPLUS : NEW ITERATE
C       FPLUS : FUNCTION VALUE AT XPLUS
C       XC    : CURRENT ITERATE
C       STEPTL: STEP TOLERANCE
C       GRADTL: GRADIENT TOLERANCE
C       RETCD : RETURN CODE WITH THE FOLLOWING MEANINGS :
C               RETCD = 0 : SUCCESSFUL GLOBAL STRATEGY
C               RETCD = 1 : UNSUCCESSFUL GLOBAL STRATEGY
C       FTOL  : FUNCTION TOLERANCE
C       ITN   : ITERATION NUMBER
C       ITNLIM: ITERATION NUMBER LIMIT
C       ICSCMX: NUMBER CONSECUTIVE STEPS .GE. STEPMX
C       MXTAKE: BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH
C       M     : DIMENSION OF FPLUS
C       N     : DIMENSION OF G, XC, AND XPLUS
C       MSG   : MESSAGE TO INHIBIT CERTAIN AUTOMATIC CHECKS + OUTPUT
C       IPR   : DEVICE TO WHICH TO SEND OUTPUT
C
C
C       OUTPUT PARAMETERS :
C       ------------------
C
C       TERMCD: RETURN CODE WITH THE FOLLOWING MEANINGS :
C              TERMCD = 0 NO TERMINATION CRITERION SATISFIED
C
C              TERMCD > 0 : SOME TERMINATION CRITERION SATISFIED
C              TERMCD = 1 : NORM OF SCALED FUNCTION VALUE IS LESS THAN
C              FTOL
C              TERMCD = 2 :  GRADIENT TOLERANCE REACHED
C              TERMCD = 3 : SCALED DISTANCE BETWEEN LAST TWO STEPS
C              LESS THAN STEPTL
C              TERMCD = 4 : UNSUCCESSFUL GLOBAL STRATEGY
C              TERMCD = 5 : ITERATION LIMIT EXCEEDED
C              TERMCD = 6 : FIVE CONSECUTIVE STEPS OF LENGTH STEPMX
C                           HAVE BEEN TAKEN
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  IDAMAX
C
C**********************************************************************

       INTEGER I
       DOUBLE PRECISION RES,D,RGX,RELGRD,RSX,RELSTP,ZERO,ONE
       INTEGER IDAMAX
       INTRINSIC ABS,MAX
       DATA ZERO,ONE/0.0D0,1.0D0/

c check whether scaled function is within tolerance

       RES = ABS(FPLUS(IDAMAX(M,FPLUS,1)))
       IF(RES.LE.FTOL) THEN
          TERMCD = 1
          IF(MOD(MSG/8,2).EQ.0) THEN
             WRITE(IPR,701)
          ENDIF
          RETURN
       ENDIF

c check whether scaled gradient is within tolerance

       D = ONE/MAX(FNORM, DBLE(N/2))
       RGX = ZERO
       DO 200 I = 1,N
        RELGRD = ABS(G(I)) * MAX(ABS(XPLUS(I)), ONE)*D
        RGX = MAX(RGX,RELGRD)
 200   CONTINUE
       IF(RGX.LE.GRADTL) THEN
          TERMCD = 2
          IF(MOD(MSG/8,2).EQ.0) THEN
             WRITE(IPR,702)
          ENDIF
          RETURN
       ENDIF

       IF(ITN.EQ.0) RETURN

       IF(RETCD.EQ.1) THEN
          TERMCD = 4
          IF(MOD(MSG/8,2).EQ.0)  THEN
          WRITE(IPR,703)
          ENDIF
          RETURN
       ENDIF

c check whether relative step length is within tolerance

       RSX = ZERO
       DO 300 I = 1,N
          RELSTP = ABS(XPLUS(I) - XC(I))/MAX(XPLUS(I), ONE)
          RSX = MAX(RSX, RELSTP)
 300   CONTINUE
       IF(RSX.LE.STEPTL) THEN
         TERMCD = 3
         IF(MOD(MSG/8,2).EQ.0) THEN
         WRITE(IPR,704)
         ENDIF
         RETURN
       ENDIF

c check iteration limit

       IF(ITN.GE.ITNLIM) THEN
          TERMCD = 5
          IF(MOD(MSG/8,2).EQ.0) THEN
             WRITE(IPR,705)
          ENDIF
       ENDIF

c check number of consecutive steps .ge. stepmx

       IF(MXTAKE) THEN
          ICSCMX = ICSCMX+1
          IF(ICSCMX.GE.5) THEN
             TERMCD = 6
             IF(MOD(MSG/8,2).EQ.0) THEN
                WRITE(IPR,706)
             ENDIF
          ENDIF
       ELSE
          ICSCMX=0
       ENDIF

 701   FORMAT(/,'  TSNSTP      FUNCTION VALUE CLOSE TO ZERO')
 702   FORMAT(/,'  TSNSTP      RELATIVE GRADIENT CLOSE TO ZERO')
 703   FORMAT(/,'  TSNSTP      LAST GLOBAL STEP FAILED TO LOCATE A',/
     +        '  TSNSTP      POINT LOWER THAN THE CURRENT ITERATE')
 704   FORMAT(/,'  TSNSTP      SUCCESSIVE ITERATES WITHIN TOLERANCE',/
     +        '  TSNSTP      CURRENT ITERATE IS PROBABLY SOLUTION')
 705   FORMAT(/,'  TSNSTP      ITERATION LIMIT EXCEEDED',/
     +        '  TSNSTP      ALGORITHM FAILED')
 706   FORMAT(/,'  TSNSTP      MAXIMUM STEP SIZE EXCEEDED 5',
     +        ' CONSECUTIVE TIMES',/
     +        '  TSNSTP      EITHER THE FUNCTION IS UNBOUNDED BELOW',/
     +        '  TSNSTP      BECOMES ASYMPTOTIC TO A FINITE VALUE',/
     +        '  TSNSTP      FROM ABOVE IN SOME DIRECTION',/
     +        '  TSNSTP      OR STEPMX IS TOO SMALL')

       RETURN
       END

      SUBROUTINE TSPRMV(X,Y,PIVOT,N,JOB)

      INTEGER N,JOB
      INTEGER PIVOT(N)
      DOUBLE PRECISION X(N),Y(N)

C**********************************************************************
C THIS SUBROUTINE PERFORMS A VECTOR PERMUTATION.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       Y :  VECTOR TO TSPRMV
C       PIVOT :  PIVOT VECTOR
C       N :  DIMENSION OF THE VECTORS Y AND PIVOT
C
C       OUTPUT PARAMETERS :
C       -------------------
C
C        X : PIVOTED VECTOR
C
C**********************************************************************

       INTEGER I

       IF(JOB .EQ. 0) THEN

c permute Y

          DO 20 I = 1,N
             X(PIVOT(I)) = Y(I)
 20       CONTINUE
       ELSE

c reverse permute of Y

          DO 30 I = 1,N
             X(I) = Y(PIVOT(I))
 30       CONTINUE

       ENDIF

       RETURN
       END

       SUBROUTINE TSRSLT(N,XP,FVAL,GP,ITN,IPR)

       INTEGER N,ITN,IPR
       DOUBLE PRECISION FVAL
       DOUBLE PRECISION XP(N),GP(N)

C**********************************************************************
C THIS ROUTINE PRINTS INFORMATION.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       M,N  : DIMENSIONS OF PROBLEM
C       XP   : NEXT ITERATE
C       FVAL : SUM OF SQUARES OF F(XP)
C       GP   : GRADIENT AT XP
C       ITN  : ITERATION NUMBER
C       IPR  : DEVICE TO WHICH TO SEND OUTPUT
C
C**********************************************************************

       INTEGER I

       WRITE(IPR,801) ITN
       WRITE(IPR,802)
       WRITE(IPR,803) (XP(I),I = 1,N)
       WRITE(IPR,804)
       WRITE(IPR,805) FVAL
       WRITE(IPR,806)
       WRITE(IPR,807) (GP(I),I = 1,N)

 801   FORMAT(/,'  TSRSLT    ITERATION K   =',I5)
 802   FORMAT('  TSRSLT    X(K)')
 803   FORMAT(100('  TSRSLT    ',3(D20.13,3X),/))
 804   FORMAT('  TSRSLT    FUNCTION AT X(K)')
 805   FORMAT('  TSRSLT       ',D20.13)
 806   FORMAT('  TSRSLT    GRADIENT AT X(K)')
 807   FORMAT(100('  TSRSLT    ',3(D20.13,3X),/))

       RETURN
       END

       SUBROUTINE TSQ1P1(AJA,ANLS,S,F,MAXM,MAXN,N,ROOT,RESTNS)

       INTEGER MAXM,MAXN,N
       DOUBLE PRECISION ROOT,RESTNS
       DOUBLE PRECISION AJA(MAXM,N),S(MAXN,*),F(N),ANLS(MAXM,*)

C**********************************************************************
C THIS ROUTINE SOLVES THE LOWER M-N+Q QUADRATIC EQUATIONS IN P UNKNOWNS
C OF THE TENSOR MODEL WHEN Q = 1 AND P = 1.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       AJA  : JACOBIAN MATRIX AT CURRENT ITERATE
C       ANLS : TENSOR TERM MATRIX AT CURRENT ITERATE
C       S    : MATRIX OF PAST LINEARLY INDEPENDENT DIRECTIONS
C       F    : FUNCTION VALUE AT CURRENT ITERATE MULTIPIED BY AN
C              ORTHOGONAL MATRIX
C       MAXM : LEADING DIMENSION OF AJA AND ANLS
C       MAXN : LEADING DIMENSION OF S
C       N    : COLUMN DIMENSION OF AJA
C
C       OUTPUT PARAMETERS :
C       -----------------
C
C       ROOT   : SOLUTION TO THE SYSTEM
C       RESTNS : TENSOR RESIDUAL
C
C**********************************************************************

       DOUBLE PRECISION DELTA,T1,T2,ZERO,HALF,ONE,TWO
       INTRINSIC ABS,SQRT
       DATA ZERO,HALF,ONE,TWO/0.0D0,0.50D0,1.0D0,2.0D0/

c find the roots of the equation:
c F(N) + AJA(N,N)*D + 0.5*ANLS(N,1)*(S(N+2,1)*D)**2

       T1 = AJA(N,N)
       T2 = ANLS(N,1) * S(N+2,1)**2
       IF(ANLS(N,1).EQ.ZERO) THEN
          ROOT = -F(N)/T1
       ELSE
          DELTA = T1**2 - TWO*F(N)*T2
          IF(DELTA.GE.ZERO) THEN
             ROOT = (-T1+SIGN(ONE,T1) * SQRT(DELTA))/T2
          ELSE
             ROOT = -T1/T2
          ENDIF
       ENDIF

c compute tensor residual

       RESTNS = ABS(F(N)+AJA(N,N)*ROOT+HALF*ANLS(N,1)*(S(N+2,1)**2)*
     +              (ROOT**2))
       RETURN
       END

        SUBROUTINE TSQFCN(P,X,SUM,AJA,ANLS,S,F,WRK1,WRK2,WRK3,
     +                    WRK4,WRK5,MAXM,MAXN,M,N,Q)

        INTEGER MAXM,MAXN,M,N,P,Q
        DOUBLE PRECISION X(P),AJA(MAXM,N),ANLS(MAXM,P),S(MAXN,P)
        DOUBLE PRECISION F(M),WRK1(M),WRK2(P),WRK3(P),WRK4(M),WRK5(M)

C*********************************************************************
C THIS ROUTINE IS USED TO EVALUATE THE RESIDUAL OF THE LAST M-N+P
C QUADRATIC EQUATIONS IN P UNKNOWNS OF THE TENSOR MODEL. NOTE THAT
C THIS ROUTINE IS CALLED BY UNCMIN TO SOLVE THE NONLINEAR LEAST SQUARES
C PART OF THE TENSOR MODEL.
C*********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       P : DIMENSION OF THE PROBLEM SOLVED BY UNCMIN
C       AJA  : JACOBIAN MATRIX AT CURRENT ITERATE
C       ANLS : TENSOR TERM MATRIX AT CURRENT ITERATE
C       S    : MATRIX OF PAST LINEARLY INDEPENDENT DIRECTIONS
C       F    : FUNCTION VALUE AT CURRENT ITERATE MULTIPLIED BY AN
C              ORTHOGONAL MATRIX
C       WRK1,WRK2,WRK3,WRK4,WRK5 : WORKING VECTORS
C       MAXM : LEADING DIMENSION OF AJA AND ANLS
C       MAXN : LEADING DIMENSION OF S
C       M,N  : DIMENSION OF PROBLEM
C       Q    : NUMERICAL RANK OF JACOBIAN :
C           Q > P : JACOBIAN IS SINGULAR
C           Q = P : OTHERWISE
C
C       INPUT-OUTPUT PARAMETERS :
C       -----------------------
C
C       X : NULL VECTOR ON ENTRY AND APPROXIMATION OF THE SOLUTION
C           TO THE SYSTEM OF M-N+Q QUADRATIC EQUATIONS IN P UNKNOWNS
C           OF THE TENSOR MODEL ON EXIT
C
C       OUTPUT PARAMETERS :
C       -----------------
C
C       SUM : RESIDUAL OF THE LAST M-N+P QUADRATIC EQUATIONS IN P
C             UNKNOWNS OF THE TENSOR MODEL
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DNRM2
C       LEVEL 2 BLAS  ...  DGEMV
C       TENSOLVE      ...  TSSTMX
C
C*********************************************************************

       INTEGER I
       DOUBLE PRECISION SUM,ZERO,HALF,ONE
       DOUBLE PRECISION DNRM2
       DATA ZERO,HALF,ONE/0.0D0,0.50D0,1.0D0/

c compute the lower right (m-n+q) x p submatrix of AJA times X

       CALL DGEMV('N',M-N+Q,P,ONE,AJA(N-Q+1,N-P+1),MAXM,
     +            X,1,ZERO,WRK1,1)

c compute S-trans times X

       CALL TSSTMX(S,X,MAXN,N,P,WRK2,WRK3)

c compute 0.5 * (S-trans times X)**2

       DO 10 I = 1, P
          WRK2(I) = HALF * WRK3(I)**2
 10    CONTINUE

c compute 0.5 * (down (m-n+q) x p submatrix of ANLS) *
c (S-trans times X)**2

       CALL DGEMV('N',M-N+Q,P,ONE,ANLS(N-Q+1,1),MAXM,
     +            WRK2,1,ZERO,WRK4,1)

       DO 20 I = 1,M-N+Q
          WRK5(I) = WRK4(I)+F(N-Q+I)+WRK1(I)
 20    CONTINUE

       SUM = HALF*DNRM2(M-N+Q,WRK5,1)**2

       RETURN
       END

      SUBROUTINE TSQLFC(QL,NR,M,N,IERR)

      INTEGER NR,M,N,IERR
      DOUBLE PRECISION QL(NR,N)

C**********************************************************************
C THIS ROUTINE PERFORMS A QL DECOMPOSITION.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       ----------------
C
C        NR  : LEADING DIMENSION OF MATRIX QL
C        M   : ROW DIMENSION OF QL
C        N   : COLUMN DIMENSION OF QL
C
C       INPUT-OUTPUT PARAMETERS :
C       -----------------------
C
C        QL : INPUT MATRIX ON ENTRY AND FACTORED MATRIX ON EXIT
C
C       OUTPUT PARAMETERS :
C       ------------------
C
C        IERR : RETURN CODE WITH THE FOLLOWING MEANINGS :
C               IERR = 1 : SINGULARITY DETECTED
C               IERR = 0 : OTHERWISE
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DAXPY,DDOT,DSCAL
C
C**********************************************************************

       INTEGER I,J,K
       DOUBLE PRECISION NU,SIGMA,TAU,RNU,ZERO,ONE
       INTRINSIC ABS,MAX
       DOUBLE PRECISION DDOT,DNRM2
       DATA ZERO,ONE/0.0D0,1.0D0/

c zero out rows m+1 and m+2 of matrix QL

       DO 10 J = 1,N
          QL(M+1,J) = ZERO
          QL(M+2,J) = ZERO
 10    CONTINUE

       IERR = 0

       K = 1

 20    CONTINUE
       IF((K.LT.M).AND.(K.LE.N)) THEN

c find NU = max element of col K on or above l-diagonal

         NU = ZERO
         DO 40 I = 1,M+1-K
            NU = MAX(NU,ABS(QL(I,K)))
 40      CONTINUE

         IF(NU.NE.ZERO) THEN

c normalize col K on or above l-diagonal

            RNU = ONE/NU
            CALL DSCAL(M-K+1,RNU,QL(1,K),1)

c code to find SIGMA = SGN(QL(M+1-K,K))*l2-norm of col K on or
c above l-diagonal

            SIGMA = DNRM2(M-K+1,QL(1,K),1)
            SIGMA = SIGN(SIGMA,QL(M+1-K,K))

c store last element(1st in normal QR) of U-vector in QL(M+1-K,K)

            QL(M+1-K,K) = QL(M+1-K,K)+SIGMA

c store value of <U,U>/2 in QL(M+1,K)

            QL(M+1,K) = SIGMA*QL(M+1-K,K)
            IF(QL(M+1,K).EQ.ZERO) THEN
               IERR = 1
               RETURN
            ENDIF

c store L(M+1-K,K) in QL(M+2,K)

            QL(M+2,K) = -NU*SIGMA

c code to get (I-2U*UT/<U,U>)*QL for kth iteration

            IF(K.LT.N) THEN
               DO 50 J = K+1,N

c loop to get TAU = <U,J-TH COL OF QL>

               TAU = DDOT(M-K+1,QL(1,K),1,QL(1,J),1)
               TAU = -TAU/QL(M+1,K)

c loop to get (I-2U*UT/<U,U>)*j-th col of QL

               CALL DAXPY(M-K+1,TAU,QL(1,K),1,QL(1,J),1)

 50            CONTINUE
            ENDIF
            K = K+1
        ELSE
            IERR = 1
            RETURN
        ENDIF

        GOTO 20

      ENDIF

      IF(M.EQ.N) QL(M+2,N) = QL(1,N)

      IF(QL(M+2,N).EQ.ZERO) THEN
         IERR = 1
      ENDIF

      RETURN
      END

       SUBROUTINE TSQMLV(NR,N,P,Q,V,W,TRANS)

       INTEGER NR,N,P
       DOUBLE PRECISION Q(NR,P),V(N),W(N)

C**********************************************************************
C THIS ROUTINE MULTIPLYS AN ORTHOGONAL MATRTIX Q OR ITS TRANSPOSE BY
C A VECTOR.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       ----------------
C
C       NR : LEADING DIMENSION OF MATRIX Q
C       N  : DIMENSION OF VECTORS V AND W
C       P  : COLUMN DIMENSION OF MATRIX Q
C       Q  : ORTHOGONAL MATRIX (OBTAINED FROM TSQLFC SUBROUTINE)
C       V  : VECTOR TO BE MULTIPLIED BY Q
C       TRANS : BOOLEAN PARAMETER:
C               = TRUE  : PERFORM Q-TRANS*V
C               = FALSE : PERFORM Q*V
C
C       OUTPUT PARAMETERS :
C       -----------------
C
C       W  : VECTOR Q*V OR Q-TRANS*V ON EXIT
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DAXPY,DCOPY,DDOT
C
C**********************************************************************

       INTEGER J,K
       DOUBLE PRECISION TAU,CONST
       LOGICAL TRANS
       DOUBLE PRECISION DDOT

       CALL DCOPY(N,V,1,W,1)

       DO 40 J = 1,P
          IF(TRANS) THEN
             K = P+1-J
          ELSE
             K = J
          ENDIF
          TAU = DDOT(N-K+1,Q(1,K),1,W,1)
          CONST = -TAU/Q(N+1,K)
          CALL DAXPY(N-K+1,CONST,Q(1,K),1,W,1)
 40    CONTINUE

       RETURN
       END

      SUBROUTINE TSQMTS(ANLS,QHAT,NR,MJ,N,M,P,LB,WRK1,ZERO1)

      INTEGER NR,M,N,P,MJ,LB,ZERO1
      DOUBLE PRECISION ANLS(NR,P),QHAT(NR,N),WRK1(M)

C**********************************************************************
C THIS ROUTINE MULTIPLIES AN ORTHOGONAL MATRIX QHAT BY THE TENSOR
C MATRIX ANLS.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       ----------------
C
C       QHAT :  ORTHOGONAL MATRIX (OBTAINED FROM TSQRFC SUBROUTINE)
C       NR : LEADIND DIMENSION OF MATRIX QHAT
C       MJ : ROW DIMENSION OF QHAT
C       N  : COLUMN DIMENSION OF MATRIX QHAT
C       M  : ROW DIMENSION OF MATRIX ANLS
C       P  : COLUMN DIMENSION OF MATRIX ANLS
C       LB : STARTING COLUMN FROM WHICH QR DECOMPOSITION WAS PERFORMED
C       WRK1 : WORKING VECTOR
C       ZERO1: FIRST ZERO COLUMN OF MATRIX QHAT IN CASE OF SINGULARITY
C
C        INPUT-OUTPUT PARAMETERS :
C       ------------------------
C
C       ANLS : MATRIX TO BE MULTIPLIED BY AN ORTHOGONAL MATRIX
C       ON ENTRY AND THE MATRIX QHAT*ANLS ON EXIT
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DCOPY
C       TENSOLVE      ...  TSQMUV
C
C**********************************************************************

       INTEGER J

       DO 40 J = 1,P
          CALL TSQMUV(QHAT,ANLS(1,J),WRK1,NR,MJ,LB,ZERO1,.FALSE.)
          CALL DCOPY(M,WRK1,1,ANLS(1,J),1)
 40    CONTINUE

       RETURN
       END

       SUBROUTINE TSQMUV(Q,V,W,NR,M,LB,ZERO1,TRANSP)

       INTEGER NR,M,LB,ZERO1
       DOUBLE PRECISION Q(NR,*),V(M),W(M)
       LOGICAL TRANSP

C**********************************************************************
C THIS SUBROUTINE MULTIPLIES AN ORTHOGONAL MATRIX BY A VECTOR.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       Q  : FACTORED MATRIX (OBTAINED FROM SUBROUTINE TSQRFC)
C       V  : VECTOR TO BE MULTIPLIED BY THE ORTHOGONAL MATRIX Q
C       NR : LEADING DIMENSION OF MATRIX Q
C       M  : ROW DIMENSION OF MATRIX Q
C       LB : STARTING COLUMN FROM WHICH QR DECOMPOSITION WAS PERFORMED
C       ZERO1: FIRST ZERO COLUMN OF THE MATRIX Q
C       TRANSP : BOOLEAN PARAMETER :
C                = TRUE  : PERFORM Q-TRANS*V
C                = FALSE : PERFORM Q*V
C
C       OUTPUT PARAMETERS :
C       -----------------
C
C       W : Q*V OR Q-TRANS*V ON EXIT
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DAXPY,DCOPY,DDOT
C
C********************************************************************
C
       INTEGER UB,A,B,C,K
       DOUBLE PRECISION TAU,CONST
       DOUBLE PRECISION DDOT

c copy the vector V to W

       CALL DCOPY(M,V,1,W,1)

       UB = ZERO1-1
       IF(TRANSP) THEN
         A = UB
         B = LB
         C = -1
       ELSE
         A = LB
         B = UB
         C = 1
       ENDIF

       DO 50 K = A,B,C
          TAU = DDOT(M-K+1,Q(K,K),1,W(K),1)
          CONST = -TAU/Q(M+1,K)
          CALL DAXPY(M-K+1,CONST,Q(K,K),1,W(K),1)
 50    CONTINUE

       RETURN
       END

       SUBROUTINE TSQRFC(QR,NR,N,M,LB,UB,IERR,EPSM,AL2NRM,CURPOS,ZERO1)

       INTEGER NR,N,M,LB,UB,IERR,ZERO1
       INTEGER CURPOS(N)
       DOUBLE PRECISION QR(NR,N),AL2NRM(M),EPSM

C**********************************************************************
C THIS ROUTINE PERFORMS COLUMN-PIVOTED QR DECOMPOSITION ON AN M*N
C MATRIX. THE DECOMPOSITION IS DONE BETWEEN COLS LB AND UB.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       NR : LEADING DIMENSION OF MATRIX QR
C       N  : COLUMN DIMENSION OF MATRIX QR
C       M  : ROW DIMENSION OF MATRIX QR
C       LB,UB : SUBSPACE OF QR DECOMPOSITION
C       EPSM  : MACHINE PRECISION
C       AL2NRM: WORKING VECTOR
C
C       INPUT-OUTPUT PARAMETERS :
C       ------------------------
C       QR  : INPUT MATRIX ON ENTRY AND FACTORED MATRIX ON EXIT
C
C       OUTPUT PARAMETERS :
C       ------------------
C
C       IERR : RETURN CODE WITH TH FOLLOWING MEANINGS:
C              IERR  =  1 : SINGULARITY DETECTED
C              IERR  =  0 : OTHERWISE
C       CURPOS :  PIVOT VECTOR
C       ZERO1  :  FIRST ZERO COLUMN OF MATRIX QR IN CASE OF SINGULARITY
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DAXPY,DDOT,DNRM2,DSCAL,DSWAP,IDAMAX
C
C **********************************************************************

       INTEGER COLPIV,I,J,K,L
       DOUBLE PRECISION COLMAX,SIGMA,TAU,AMAX
       DOUBLE PRECISION NU,RNU,ZERO,ONE
       DOUBLE PRECISION DDOT,DNRM2
       INTEGER IDAMAX
       INTRINSIC ABS,SIGN
       DATA ZERO,ONE/0.0D0,1.0D0/

c zero rows m+1 and m+2 of QR matrix

       DO 10 J = 1,N
          CURPOS(J) = J
 10    CONTINUE

       DO 20 J = LB,UB
          QR(M+1,J) = ZERO
          QR(M+2,J) = ZERO
 20    CONTINUE

       IERR = 0
       ZERO1 = UB+1
       K = LB

c  get L2NORM**2 of columns (LB to UB)

       DO 30 J = K,UB
          AL2NRM(J) = DNRM2(M-K+1,QR(K,J),1)**2
 30    CONTINUE

 40    CONTINUE
       IF((K.LT.M).AND.(K.LE.UB)) THEN

         AMAX = ZERO
         DO 60 J = K,UB
            IF(AL2NRM(J).GE.AMAX) THEN
               AMAX = AL2NRM(J)
               COLPIV = J
            ENDIF
 60      CONTINUE

         IF(AMAX.EQ.ZERO) THEN
            IERR = 1
            ZERO1 = K
            RETURN
         ENDIF

         IF(K.EQ.LB) THEN
            COLMAX = AMAX
         ENDIF

         IF(AL2NRM(COLPIV).LE.EPSM*COLMAX) THEN
            IERR = 1
            ZERO1 = K
            RETURN
         ENDIF

         IF(COLPIV .NE. K) THEN
            CALL DSWAP(M+2,QR(1,COLPIV),1,QR(1,K),1)
            L = CURPOS(K)
            CURPOS(K) = CURPOS(COLPIV)
            CURPOS(COLPIV) = L
            CALL DSWAP(1,AL2NRM(COLPIV),1,AL2NRM(K),1)
         ENDIF

c find NU = max element of col K on or below diagonal

         L = IDAMAX(M-K+1,QR(K,K),1) + K - 1
         NU = ABS(QR(L,K))

         IF(NU.EQ.ZERO) THEN
            IERR = 1
            ZERO1 = K
            RETURN
         ENDIF

c normalize col K on or below diagonal

         RNU = ONE/NU
         CALL DSCAL(M-K+1,RNU,QR(K,K),1)

c code to find SIGMA = SGN(QR(K,K))*l2-norm of col K on or
c below diagonal

         SIGMA = DNRM2(M-K+1,QR(K,K),1)
         SIGMA = SIGN(SIGMA,QR(K,K))

c store 1st element of U-vector in QR(K,K)

         QR(K,K) = QR(K,K)+SIGMA

c store value of <U,U>/2 in QR(M+1,K)

         QR(M+1,K) = SIGMA*QR(K,K)
         IF(QR(M+1,K).EQ.ZERO) THEN
            IERR = 1
            ZERO1 = K
            RETURN
         ENDIF

c store R(K,K) in QR(M+2,K)

         QR(M+2,K) = -NU*SIGMA
         IF(QR(M+2,K).EQ.ZERO) THEN
            IERR = 1
            ZERO1 = K
            RETURN
         ENDIF

c code to get (I-2U*UT/<U,U>)*QR for kth iteration

         IF(K.LT.N) THEN
            DO 130 J = K+1,N

c loop to get UT*J-TH col of QR

               TAU = DDOT(M-K+1,QR(K,K),1,QR(K,J),1)
               TAU = -TAU/QR(M+1,K)

c loop to get (I-2U*UT/<U,U>)*j-th col of QR

               CALL DAXPY(M-K+1,TAU,QR(K,K),1,QR(K,J),1)

 130        CONTINUE
         ENDIF

c update l2norm**2 (K+1 to UB)

         DO 140 I = K+1,UB
           AL2NRM(I) = AL2NRM(I)-QR(K,I)**2
 140     CONTINUE

         K = K+1
         GOTO 40

       ELSE

         IF(LB.EQ.UB) COLMAX = AL2NRM(LB)

       ENDIF

       IF(M.EQ.UB) QR(M+2,UB) = QR(M,M)
       IF(ABS(QR(M+2,UB)).LE.EPSM*COLMAX) THEN
          IERR = 1
          ZERO1 = UB
       ENDIF

       RETURN
       END

      SUBROUTINE TSRSID(ITN,METHOD,FQ,D,CURPOS,PIVOT,PBAR,AJA,ANLS,
     +                  SHAT,FLAG,NWTAKE,IERR,MAXM,MAXN,M,N,P,WRK1,
     +                  VN,VNP,VNS,SCRES)

      INTEGER MAXM,MAXN,M,N,P,IERR,FLAG,ITN,METHOD
      INTEGER CURPOS(N),PIVOT(N),PBAR(N)
      DOUBLE PRECISION SCRES,D(N),VN(M),VNP(M),VNS(M),AJA(MAXM,N)
      DOUBLE PRECISION WRK1(M),SHAT(MAXN,P),FQ(M)
      DOUBLE PRECISION ANLS(MAXM,P)
      LOGICAL NWTAKE

C**********************************************************************
C THIS ROUTINE COMPUTES || F + J*D + 0.5*A*D**2 ||**2 IN THE L2
C NORM SENS AT A GIVEN STEP D.
C**********************************************************************
C
C       INPUT PARAMETERS
C       ----------------
C
C        ITN   : CURRENT ITERATION NUMBER
C        METHOD: METHOD TO BE USED
C        FQ    : FUNCTION VALUE AT CURRENT ITERATE MULTIPLIED BY
C                ORTHOGONAL MATRICES
C        D     : STEP AT WHICH TO EVALUATE THE TENSOR MODEL
C        CURPOS: PIVOT VECTOR (USED DURING THE FACTORIZATION OF AJA
C                FROM COLUMN 1 TO N-P)
C        PIVOT : PIVOT VECTOR ( USED DURING THE FACTORIZATION OF AJA
C                FROM COLUMN N-P+1 TO N)
C        PBAR  : PIVOT VECTOR (USED DURING THE FACTORIZATION OF AJA
C                IF IT IS SINGULAR
C        AJA   : JACOBIAN MATRIX AT CURRENT ITERATE
C        ANLS  : TENSOR TERM MATRIX AT CURRENT ITERATE
C        SHAT  : MATRIX OF PAST LINEARLY INDEPENDENT DIRECTIONS AFTER
C                A QL FACTORIZATION
C        FLAG  : RETURN CODE WITH THE FOLLOWING MEANINGS:
C                FLAG = 0 : NO SINGULARITY DETECTED DURING FACTORIZATION
C                           OF THE JACOBIAN FROM COLUMN 1 TO N
C                FLAG = 1 : SINGULARITY DETECTED DURING FACTORIZATION
C                           OF THE JACOBIAN FROM COLUMN 1 TO N-P
C                FLAG = 2 : SINGULARITY DETECTED DURING FACTORIZATION
C                           OF THE JACOBIAN FROM COLUMN N-P+1 TO N
C        NWTAKE: LOGICAL VARIABLE WITH THE FOLLOWING MEANINGS:
C                NWTAKE = .TRUE.  : NEWTON STEP TAKEN
C                NWTAKE = .FALSE. : TENSOR STEP TAKEN
C        IERR  : RETURN CODE FROM QRP FACTORIZATION ROUTINE:
C                IERR = 0 : NO SINGULARITY DETECTED
C                IERR = 1 : SINGULARITY DETECTED
C        MAXM   : LEADING DIMENSION OF AJA AND ANLS
C        MAXN   : LEADING DIMENSION OF SHAT
C        M,N    : DIMENSIONS OF THE PROBLEM
C        P      : COLUMN DIMENSION OF THE MATRICES SHAT AND ANLS
C        WRK1,VN,VNP,VNS : WORKSPACE VECTORS
C
C       OUTPUT PARAMETERS
C       -----------------
C
C        SCRES :  || F + J*D + 0.5*A*D**2 ||**2
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DNRM2
C       LEVEL 2 BLAS  ...  DGEMV
C       TENSOLVE      ...  TSDLOD,TSJMUV,TSUDQV
C
C **********************************************************************

      INTEGER I,J,LEN
      DOUBLE PRECISION ZERO,HALF,ONE
      DOUBLE PRECISION DNRM2
      DATA ZERO,HALF,ONE/0.0D0,0.50D0,1.0D0/

      CALL TSJMUV(ITN,METHOD,D,CURPOS,PIVOT,PBAR,AJA,SHAT,FLAG,
     +            IERR,MAXM,MAXN,M,N,P,VN,VNP,VNS,WRK1)

      CALL TSDLOD (M,ZERO,WRK1(N+1),1)

      LEN = M
      IF(IERR .GT. 0) LEN = M + N

      DO 20 I = 1, LEN
         VN(I) = WRK1(I) + FQ(I)
 20   CONTINUE

      IF( .NOT. NWTAKE) THEN
        CALL TSUDQV(SHAT,VNS,MAXN,N,P,VNP)
        DO 30 J = 1, P
           VNP(J) = VNP(J) ** 2
 30     CONTINUE
        CALL DGEMV('N',LEN,P,HALF,ANLS,MAXM,VNP,1,ONE,VN,1)
      ENDIF

      SCRES = DNRM2(LEN,VN,1)

      RETURN
      END

       SUBROUTINE TSSCLF(F,DF,M)

       INTEGER M
       DOUBLE PRECISION F(M),DF(M)

C*******************************************************************
C THIS ROUTINE SCALES A FUNCTION VALUE F.
C*******************************************************************
C
C       INPUT PARAMETERS :
C       ------------------
C
C       DF : DIAGONAL SCALING MATRIX FOR F
C       M  : DIMENSION OF F
C
C       INPUT-OUTPUT PARAMETERS :
C       ------------------
C
C       F  : UNSCALED FUNCTION VALUE ON ENTRY AND SCALED FUNCTION
C            VALUE ON EXIT
C
C*********************************************************************

       INTEGER I

       DO 10 I = 1,M
          F(I) = DF(I)*F(I)
 10    CONTINUE

       RETURN
       END

       SUBROUTINE TSSCLJ(X,DX,TYPX,F,DF,FHAT,NR,M,N,EPSM,JACFLG,
     +                   FVEC,JAC,AJA)

       INTEGER NR,M,N,JACFLG
       DOUBLE PRECISION EPSM
       DOUBLE PRECISION X(N),DX(N),TYPX(N),F(M),DF(M)
       DOUBLE PRECISION AJA(NR,N),FHAT(M)
       EXTERNAL FVEC,JAC

C**********************************************************************
C THIS ROUTINE COMPUTES THE JACOBIAN MATRIX AT THE CURRENT ITERATE
C X AND SCALES ITS VALUE.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       X    : SCALED CURRENT ITERATE
C       DX   : DIAGONAL SCALING MATRIX FOR X
C       F    : SCALED FUNCTION VALUE AT X
C       DF   : DIAGONAL SCALING MATRIX FOR F
C       FHAT : WORKSPACE ARRAY
C       NR   : LEADING DIMENSION OF AJA
C       M,N  : DIMENSIONS OF PROBLEM
C       EPSM : MACHINE PRECISION
C       JACFLG : JACOBIAN FLAG
C       FVEC : SUBROUTINE TO EVALUATE FUNCTION
C       JAC  : SUBROUTINE TO EVALUATE ANALYTIC JACOBIAN
C
C
C       INPUT-OUTPUT PARAMETERS :
C       ------------------------
C
C       AJA  : SCALED JACOBIAN AT CURRENT ITERATE
C
C       SUBPROGRAMS CALLED:
C
C       TENSOLVE      ...  TSUNSX,TSUNSF,TSFDFJ,TSSCLF,TSSCLX
C       USER          ...  FVEC,JAC
C
C********************************************************************

       INTEGER I,J

c unscale X AND F

       CALL TSUNSX(X,DX,N)
       CALL TSUNSF(F,DF,M)

c compute the finite difference or analytic Jacobian at X

       IF(JACFLG.EQ.0) THEN
          CALL TSFDFJ(X,F,NR,M,N,EPSM,FVEC,FHAT,AJA)
       ELSE
          CALL JAC(X,AJA,NR,M,N)
       ENDIF

c scale the Jacobian matrix

       DO 20 J = 1,N
          DO 10 I = 1,M
             AJA(I,J) = AJA(I,J)*DF(I)*TYPX(J)
 10       CONTINUE
 20    CONTINUE

c scale back X AND F

       CALL TSSCLF(F,DF,M)
       CALL TSSCLX(X,DX,N)

       RETURN
       END

       SUBROUTINE TSSCLX(X,DX,N)

       INTEGER N
       DOUBLE PRECISION X(N),DX(N)

C**********************************************************************
C THIS ROUTINE SCALES A VECTOR X.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       ------------------
C
C       DX : DIAGONAL SCALING MATRIX FOR X
C       N  : DIMENSION OF X
C
C       OUTPUT PARAMETERS :
C       ------------------
C
C       X  : SCALED VECTOR X
C
C**********************************************************************

       INTEGER I

       DO 10 I = 1,N
          X(I) = DX(I)*X(I)
 10    CONTINUE

       RETURN
       END

       SUBROUTINE TSSLCT(RESIDT,RESIDN,ITRMCD,FC,M,N,DN,DT,G,DF,NWTAKE)

       INTEGER M,N,ITRMCD
       DOUBLE PRECISION RESIDT,RESIDN
       DOUBLE PRECISION FC(M),DF(N),DN(N),DT(N),G(N)
       LOGICAL NWTAKE

C*********************************************************************
C THIS ROUTINE DECIDES WHICH DIRECTION TO CHOOSE: THE TENSOR OR THE
C STANDARD DIRECTION. THE STANDARD DIRECTION IS CHOSEN WHENEVER
C THE TENSOR DIRECTION IS NOT DESCENT OR THE TENSOR DIRECTION IS TO A
C MINIMIZER OF THE TENSOR MODEL AND DOESN'T PROVIDE ENOUGH DECREASE
C IN THE TENSOR MODEL, OR THE QUADRATIC SYSTEM OF Q EQUATIONS IN P
C UNKNOWNS CANNOT BE SOLVED BY UNCMIN WITHIN 150 ITERATIONS.
C*********************************************************************
C
C       INPUT PARAMETERS :
C       ------------------
C
C       RESIDT : TENSOR RESIDUAL
C       RESIDN : NEWTON RESIDUAL
C       ITRMCD : UNCMIN TERMINATION CODE
C       FC : FUNCTION VALUE AT XC
C       M,N: DIMENSIONS OF PROBLEM
C       DN : STANDARD STEP
C       DT : TENSOR STEP
C       G  : GRADIENT VALUE AT XC
C
C
C       OUTPUT PARAMETERS :
C       -----------------
C
C       DF : EITHER THE STANDARD OR TENSOR STEP ON EXIT
C       NWTAKE : BOOLEAN VALUE WITH THE FOLLOWING MEANINGS:
C               =.TRUE.  : STANDARD STEP IS TAKEN
C               =.FALSE. : TENSOR STEP IS TAKEN
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ....  DCOPY,DDOT,DNRM2
C
C*********************************************************************

       DOUBLE PRECISION ANRMFC,DTNORM,GNORM
       DOUBLE PRECISION TEMP,TEMP1,BETA,GAMA
       DOUBLE PRECISION TENTH,ONETT,HALF
       DOUBLE PRECISION DNRM2,DDOT
       DATA ONETT,TENTH,HALF/1.0D-4,1.0D-1,0.50D0/

       NWTAKE = .FALSE.
       ANRMFC = DNRM2(M,FC,1)
       DTNORM = DNRM2(N,DT,1)
       GNORM = DNRM2(N,G,1)
       TEMP = DDOT(N,DT,1,G,1)

       GAMA = HALF
       IF(M.GT.N) THEN
          BETA = TENTH
       ELSE
          BETA = ONETT
       ENDIF

       TEMP1 = -BETA*DTNORM*GNORM

       IF(RESIDT.GE.GAMA*(ANRMFC+RESIDN).OR.(TEMP.GT.TEMP1).OR.
     +    ITRMCD.EQ.4) THEN
          CALL DCOPY(N,DN,1,DF,1)
          NWTAKE = .TRUE.
       ELSE
          CALL DCOPY(N,DT,1,DF,1)
       ENDIF

       RETURN
       END

       SUBROUTINE TSSMIN(ANLS,FQ,ADT,AG,CONST1,CONST2,DLT,NR,M,N,
     +                   P,NWTAKE,IERR,EPSM,VN,VNP,VNS,SOL)

       DOUBLE PRECISION DLT,EPSM
       INTEGER NR,M,N,P,IERR
       DOUBLE PRECISION ADT(N),AG(N),VN(M),VNP(M)
       DOUBLE PRECISION VNS(M),ANLS(NR,P),FQ(M)
       DOUBLE PRECISION CONST1(P),CONST2(P)
       LOGICAL NWTAKE

C***********************************************************************
C THIS ROUTINE MINIMIZES THE TENSOR MODEL OVER THE SUBSPACE SPANNED BY
C THE TENSOR STEP AND THE STEEPEST DECENT DIRECTION. IF THE NEWTON STEP
C WERE CHOSEN, IT WILL MINIMIZE THE NEWTON MODEL OVER THE SUBSPACE
C SPANNED BY THE NEWTON STEP AND THE STEEPEST DESCENT DIRECTION.
C***********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       ANLS : TENSOR TERM MATRIX AT CURRENT ITERATE
C       FQ   : FUNCTION VALUE AT CURRENT ITERATE MULTIPLIED BY
C              ORTHOGONAL MATRICES
C       ADT  : JACOBIAN MATRIX TIMES DT (SEE SUBROUTINE TS2DTR)
C        AG  : JACOBIAN MATRIX TIMES GBAR (SEE SUBROUTINE TS2DTR)
C      CONST1: SHAT-TRANS * DT  (SEE SUBROUTINE TS2DTR)
C      CONST2: SHAT-TRANS * GBAR (SEE SUBROUTINE TS2DTR)
C      ALPHA : POINT AT WHICH DERIVATIVE IS EVALUATED
C        DLT : CURRENT TRUST RADIUS
C        NR  : LEADING DIMENSION OF ANLS
C         M,N: DIMENSIONS OF THE PROBLEM
C           P: COLUMN DIMENSION OF MATRIX ANLS
C       NWTAKE : LOGICAL VARIABLE WITH THE FOLLOWING MEANINGS :
C                NWTAKE = .TRUE.  : STANDARD STEP TAKEN
C                NWTAKE = .FALSE. : TENSOR STEP TAKEN
C       IERR   : RETURN CODE FROM QRP FACTORIZATION ROUTINE
C                IERR = 0 : NO SINGULARITY OF JACOBIAN DETECTED
C                IERR = 1 : SINGULARITY OF JACOBIAN DETECTED
C       EPSM   : MACHINE PRECISION
C       VN,VNP,VNS : WORKING VECTORS
C
C
C       OUTPUT PARAMETERS :
C       -----------------
C
C       SOL   : GLOBAL MINIMIZER OF THE ONE VARIABLE FUNCTION IN ALPHA
C               ||F + J*(ALPHA*DT + BETA*GBAR) + 0.5*A*(ALPHA*DT +
C               BETA*GBAR)**2||**2 WHERE BETA = SQRT(DLT**2 - ALPHA**2)
C
C       SUBPROGRAMS CALLED:
C
C       TENSOLVE      ...  TSFAFA,TSMFDA,TSLMIN
C
C**********************************************************************

       INTEGER INT
       DOUBLE PRECISION SOL,TOL,DL,S,SP,C,TSFAFA,A,TSMFDA
       DOUBLE PRECISION D,S1,B,Q,BC,OPTIM,AC,GLOPT,BLOOP,ELOOP,INCR
       DOUBLE PRECISION ZERO,OHUND,TENTH,TWO,THREE,TEN
       LOGICAL FIRST
       DATA ZERO,OHUND,TENTH,TWO,THREE,TEN/0.0D0,1.0D-2,1.0D-1,2.0D0,
     + 3.0D0,10.0D0/

       FIRST = .TRUE.
       TOL = EPSM**(TWO/THREE)
       INT = 40
       DL = TOL
       IF(DLT.LT.TOL) THEN
          DL = TOL*TENTH
       ELSEIF(DLT.GT.TOL*TEN) THEN
          DL = TOL*TEN
       ENDIF
       IF(DLT.LE.OHUND) THEN
          INT = 10
       ENDIF

c find global minimizer of FALPHA

       BLOOP = -DLT+DL
       ELOOP = DLT*(INT-2)/INT
       INCR = TWO*(DLT-DL)/INT
       S = BLOOP
 10    CONTINUE

       SP = S
       S1 = S+INCR

c evaluate FALPHA(SP) and the derivative of FALPHA at SP

       IF(FIRST) THEN
          C = TSFAFA(ANLS,FQ,ADT,AG,CONST1,CONST2,SP,DLT,
     +               NR,M,N,P,NWTAKE,IERR,VN)
          A = TSMFDA(ANLS,ADT,AG,CONST1,CONST2,SP,DLT,
     +               NR,M,N,P,NWTAKE,IERR,VN,VNP)
       ELSE
         C = D
         A = B
       ENDIF

c evaluate FALPHA(S1) and the derivative of FALPHA at S1

       D = TSFAFA(ANLS,FQ,ADT,AG,CONST1,CONST2,S1,DLT,
     +            NR,M,N,P,NWTAKE,IERR,VN)
       B = TSMFDA(ANLS,ADT,AG,CONST1,CONST2,S1,DLT,
     +            NR,M,N,P,NWTAKE,IERR,VN,VNP)

c minimize FALPHA in the subinterval [SP,S1]

       IF((A.LE.ZERO).AND.(B.GE.ZERO)) THEN
          IF(C.GT.D) THEN
             Q = D
             BC = B
             CALL TSLMIN(S1,SP,BC,Q,ANLS,FQ,ADT,AG,CONST1,CONST2,
     +                   DLT,NR,M,N,P,NWTAKE,IERR,TOL,VN,VNP,
     +                   VNS,OPTIM)
          ELSE
             Q = C
             AC = A
             CALL TSLMIN(SP,S1,AC,Q,ANLS,FQ,ADT,AG,CONST1,CONST2,
     +                   DLT,NR,M,N,P,NWTAKE,IERR,TOL,VN,VNP,
     +                   VNS,OPTIM)
          ENDIF
        ELSEIF((A.LE.ZERO).AND.(B.LE.ZERO)) THEN
          IF(C.LE.D) THEN
             Q = C
             AC = A
             CALL TSLMIN(SP,S1,AC,Q,ANLS,FQ,ADT,AG,CONST1,CONST2,
     +                   DLT,NR,M,N,P,NWTAKE,IERR,TOL,VN,VNP,
     +                   VNS,OPTIM)
          ELSE
             OPTIM = S1
             Q = D
          ENDIF
       ELSEIF((A.GE.ZERO).AND.(B.GE.ZERO)) THEN
          IF(C.GE.D) THEN
             Q = D
             BC = B
             CALL TSLMIN(S1,SP,BC,Q,ANLS,FQ,ADT,AG,CONST1,CONST2,
     +                   DLT,NR,M,N,P,NWTAKE,IERR,TOL,VN,VNP,
     +                   VNS,OPTIM)
          ELSE
             OPTIM = SP
             Q = C
          ENDIF
       ENDIF

c update the global minimizer of FALPHA so far

       IF(FIRST) THEN
         IF(A.GT.ZERO .AND. B.LT.ZERO) THEN
            GLOPT = C
            SOL = SP
            IF(C.GT.D) THEN
               GLOPT = D
               SOL = S1
            ENDIF
         ELSE
            GLOPT = Q
            SOL = OPTIM
         ENDIF
         FIRST = .FALSE.
       ELSEIF(GLOPT.GE.Q) THEN
         GLOPT = Q
         SOL = OPTIM
       ENDIF

       S = S + INCR
       IF(S .LE. ELOOP) GO TO 10

       RETURN
       END

       SUBROUTINE TSSMRD(VECT,RESNEW,X,MU,IERR,M,N)

       INTEGER M,N,IERR
       DOUBLE PRECISION RESNEW,MU
       DOUBLE PRECISION VECT(M),X(N)

C**********************************************************************
C THIS ROUTINE COMPUTES THE RESIDUAL OF THE STANDARD MODEL.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       VECT : RIGHT HAND SIDE VECTOR OF THE NEWTON/GAUSS-NEWTON
C              EQUATIONS AFTER BEING MULTIPLIED BY ORTHOGONAL MATRICES
C              (SEE SUBROUTINE TSCPSS)
C       X    : STANDARD STEP COMPUTED BY THE SUBROUTINE TSCPSS
C       MU   : A SMALL PERTURBATION USED IN COMPUTING THE STANDARD
C              STEP WHEN THE JACOBIAN IS SINGULAR
C       IERR : RETURN CODE WITH THE FOLLOWING MEANINGS :
C              IERR = 0 : NO SINGULARITY OF JACOBIAN DETECTED
C              IERR = 1 : OTHERWISE
C       M,N  : DIMENSION OF PROBLEM
C
C       OUTPUT PARAMETERS :
C       ------------------
C
C       RESNEW : RESIDUAL OF THE STANDARD MODEL
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DNRM2
C
C**********************************************************************

       DOUBLE PRECISION TEMP,PROD,ZERO
       DOUBLE PRECISION DNRM2
       INTRINSIC SQRT
       DATA ZERO/0.0D0/

       IF(IERR .EQ .0) THEN
          IF(M .EQ. N) THEN
             RESNEW = ZERO
          ELSE
             RESNEW = DNRM2(M-N,VECT(N+1),1)
          ENDIF
       ELSE
          TEMP = DNRM2(M,VECT(N+1),1)**2
          PROD = MU * DNRM2(N,X,1)**2
          RESNEW = SQRT(TEMP-PROD)
       ENDIF

       RETURN
       END

       SUBROUTINE TSSQP1(AJA,ANLS,S,F,MAXM,MAXN,M,N,Q,ROOT,RESTNS)

       INTEGER MAXM,MAXN,M,N,Q
       DOUBLE PRECISION ROOT,RESTNS
       DOUBLE PRECISION AJA(MAXM,N),S(MAXN,*),ANLS(MAXM,*),F(M)

C**********************************************************************
C THIS ROUTINE SOLVES THE LOWER M-N+Q QUADRATIC EQUATIONS IN P UNKNOWNS
C OF THE TENSOR MODEL WHEN Q > 1 AND P = 1.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       AJA  : JACOBIAN MATRIX AT CURRENT ITERATE
C       ANLS : TENSOR TERM MATRIX AT CURRENT ITERATE
C       S    : MATRIX OF PAST LINEARLY INDEPENDENT DIRECTIONS
C       F    : FUNCTION VALUE AT CURRENT ITERATE MULTIPLIED BY AN
C              ORTHOGONAL MATRIX
C       MAXM : LEADING DIMENSION OF AJA AND ANLS
C       MAXN : LEADING DIMENSION OF S
C       M,N  : ROW AND COLUMN DIMENSIONS OF AJA
C       Q    : NUMERICAL RANK OF JACOBIAN :
C              Q > P : JACOBIAN IS SINGULAR
C              Q = P : OTHERWISE
C
C
C       OUTPUT PARAMETERS :
C       -----------------
C
C       ROOT   : SOLUTION TO THE SYSTEM
C       RESTNS : TENSOR RESIDUAL
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DDOT,DNRM2
C
C **********************************************************************

       INTEGER I
       DOUBLE PRECISION TEMP,A,B,C,D,RES1,RES2,RES3,RES,S1,S2,S3
       DOUBLE PRECISION T,T0,T1,T2,T3,PI,A1,A2,A3,THETA
       DOUBLE PRECISION ZERO,QUART,HALF,ONE,TWO,THREE,FOUR,NINE
       DOUBLE PRECISION TSEVEN,FFOUR,ONETRD
       DOUBLE PRECISION DDOT,DNRM2
       INTRINSIC ABS,SQRT
       PARAMETER (PI = 3.141592653589793D0)
       DATA ZERO,QUART,HALF,ONE,TWO,THREE,FOUR,NINE,TSEVEN,FFOUR/0.0D0,
     + 0.250D0,0.50D0,1.0D0,2.0D0,3.0D0,4.0D0,9.0D0,27.0D0,54.0D0/

c compute the coefficients of a third degree polynomial

       ONETRD = ONE/THREE
       A = ZERO
       B = ZERO
       C = ZERO

       TEMP = DNRM2(M-N+Q,F(N-Q+1),1)**2
       D = TWO*DDOT(M-N+Q,AJA(N-Q+1,N),1,F(N-Q+1),1)
       T0 = S(N+2,1)**2
       T1 = T0**2
       DO 10 I = N-Q+1,M
          T2 = AJA(I,N)
          T3 = ANLS(I,1) * T0
          C = C + TWO * (T2**2 + F(I) * T3)
          B = B + THREE * T2 * T3
          A = A + ANLS(I,1)**2 * T1
 10    CONTINUE

c compute the roots of the third degree polynomial

       IF(A.EQ.ZERO) THEN
         IF(B.NE.ZERO) THEN
            T0 = SQRT(C**2-FOUR*B*D)
            T1 = TWO*B
            S1 = (-C+T0)/T1
            S2 = (-C-T0)/T1
            RES1 = ABS(TEMP+D*S1+HALF*C*(S1**2)+ONETRD*B*(S1**3))
            RES2 = ABS(TEMP+D*S2+HALF*C*(S2**2)+ONETRD*B*(S2**3))
            IF(RES1 .GT. RES2) THEN
               ROOT =  S2
               RES  =  RES2
            ELSE
               ROOT =  S1
               RES  =  RES1
            ENDIF
            RESTNS  =  SQRT(RES)
            RETURN
         ELSEIF(C.NE.ZERO) THEN
            ROOT = -D/C
            RES = ABS(TEMP+D*ROOT+HALF*C*(ROOT**2))
            RESTNS = SQRT(RES)
            RETURN
         ELSE
            ROOT = ZERO
            RESTNS = SQRT(TEMP)
            RETURN
         ENDIF
       ELSEIF(D.EQ.ZERO) THEN
         ROOT = ZERO
         RESTNS = SQRT(TEMP)
         RETURN
       ENDIF
       T3 = D

       A1 = B/A
       A2 = C/A
       A3 = D/A
       T0 = (THREE*A2-A1**2)/NINE
       T1 = (NINE*A1*A2-TSEVEN*A3-TWO*A1**3)/FFOUR
       D = T0**3 + T1**2

       IF(D.GT.0) THEN
          T2 = T1+SQRT(D)
          T = T1-SQRT(D)
          IF(T.LT.0) THEN
             T = -(-T)**ONETRD
          ELSE
             T = T**ONETRD
          ENDIF
          IF(T2.LT.0)THEN
             T2 = -(-T2)**ONETRD
          ELSE
             T2 = T2**ONETRD
          ENDIF
          S1 = T2+T-A1/THREE
          S3 = S1
          S2 = S1
       ELSE
          T = T1/SQRT(-T0**3)
          THETA = ACOS(T)
          THETA = THETA/THREE
          T = TWO*SQRT(-T0)
          S1 = T*COS(THETA)-A1/THREE
          S2 = T*COS(THETA+PI*TWO/THREE)-A1/THREE
          S3 = T*COS(THETA+PI*FOUR/THREE)-A1/THREE
       ENDIF

c compute the tensor residual for each root

       RES1 = ABS(TEMP+T3*S1+HALF*C*(S1**2)+ONETRD*B*(S1**3)+
     +        QUART*A*(S1**4))
       RES2 = ABS(TEMP+T3*S2+HALF*C*(S2**2)+ONETRD*B*(S2**3)+
     +        QUART*A*(S2**4))
       RES3 = ABS(TEMP+T3*S3+HALF*C*(S3**2)+ONETRD*B*(S3**3)+
     +        QUART*A*(S3**4))

c select the root that produces the smallest tensor residual

       RES  =  RES1
       ROOT  =  S1
       IF(RES .GT. RES2) THEN
           RES  =  RES2
           ROOT  =  S2
        ENDIF
        IF(RES .GT. RES3) THEN
           RES  =  RES3
           ROOT  =  S3
        ENDIF
        RESTNS  =  SQRT(RES)

       RETURN
       END

       SUBROUTINE TSSSTP(AJA,FN,M,N,NR,EPSM,IGLOBL,AL2NRM,Y,W,
     +                   DN,FQ,PIVOT,PBAR,IERR)

       INTEGER NR,M,N,IERR,IGLOBL
       INTEGER PIVOT(N),PBAR(N)
       DOUBLE PRECISION EPSM,AJA(NR,N),AL2NRM(M),FN(M)
       DOUBLE PRECISION DN(N),Y(N),W(M),FQ(M)

C**********************************************************************
C THIS ROUTINE FINDS THE STANDARD STEP WHEN THE ITERATION NUMBER IS
C EQUAL TO 1 OR THE INPUT PARAMETER "METHOD" IS SET TO 0.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       AJA   : JACOBIAN MATRIX AT CURRENT ITERATE
C       FN    : FUNCTION VALUE AT CURRENT ITERATE
C       M,N   : DIMENSIONS OF PROBLEM
C       NR    : LEADING DIMENSION OF AJA
C       EPSM  : MACHINE EPSILON
C       IGLOBL: GLOBAL STRATEGY USED :
C                = 0 : LINE SEARCH USED
C                = 1 : 2-DIMENSIONAL TRUST REGION USED
C       AL2NRM,Y,W : WORKING VECTORS
C
C
C       OUTPUT PARAMETERS :
C       ------------------
C
C       DN : STANDARD STEP
C       FQ : FUNCTION VALUE AT CURRENT ITERATE MULTIPLIED BY
C              ORTHOGONAL MATRICES (THIS IS USED IN THE CASE WHERE
C              THE GLOBAL STRATEGY IS THE 2-DIMENSIONAL
C              TRUST REGION)
C       PIVOT,PBAR : PIVOT VECTORS
C       IERR : RETURNED CODE WITH THE FOLLOWING MEANING :
C              IERR  =  1 : SINGULARITY OF JACOBIAN DETECTED (ZERO1
C                           IS USED TO KEEP TRACK OF THE FIRST COLUMN
C                           WHERE SINGULARITY IS DETECTED)
C              IERR  =  0 : OTHERWISE
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DCOPY,DSCAL
C       TENSOLVE      ...  TSDLOD,TSQRFC,TSQMUV,TSBSLV,TSPRMV,TSCPMU
C
C**********************************************************************

       INTEGER ZERO1,ZEROTM,I,J
       DOUBLE PRECISION MU,ZERO,ONE
       DATA ZERO,ONE/0.0D0,1.0D0/

c perform a QR factorization of AJA

       CALL TSQRFC(AJA,NR,N,M,1,N,IERR,EPSM,AL2NRM,PIVOT,ZERO1)

       CALL DSCAL(M,-ONE,FN,1)

       IF(IERR.EQ.0) THEN
          IF(M.EQ.N) THEN
             ZEROTM = ZERO1-1
          ELSE
             ZEROTM = ZERO1
          ENDIF

c multiply (-FN) by the orthogonal matrix resulting from the QR
c decomposition of AJA

          CALL TSQMUV(AJA,FN,W,NR,M,1,ZEROTM,.FALSE.)

c solve AJA*DN  =  W

          CALL TSBSLV(AJA,NR,M,N,W,Y)
          CALL TSPRMV(DN,Y,PIVOT,N,0)

          IF(IGLOBL.EQ.1) THEN
             IERR = 0
             CALL DCOPY(M,W,1,FQ,1)
             CALL DSCAL(M,-ONE,FQ,1)
          ENDIF
          RETURN
       ELSE

c AJA is singular

          CALL TSQMUV(AJA,FN,W,NR,M,1,ZERO1,.FALSE.)

c solve ( AJA-trans AJA + MU I ) DN = -AJA-trans FN

c put the diagonal elements stored in row m+2 of AJA into their
c propre positions and zero out the unwanted portions of AJA

          DO 10 J = 1, ZERO1-1
             AJA(J,J) = AJA(M+2,J)
             CALL TSDLOD (M+N-J,ZERO,AJA(J+1,J),1)
 10       CONTINUE

          DO 20 J = ZERO1, N
             CALL TSDLOD (M+N-ZERO1+1,ZERO,AJA(ZERO1,J),1)
 20       CONTINUE

c compute a perturbation MU

          CALL TSCPMU(AJA,NR,N,EPSM,MU)

c form the augmented Jacobian matrix by adding an nxn diag(mu) in
c the bottom of AJA

          DO 70 I = M+1,M+N
             AJA(I,I-M) = MU
 70       CONTINUE

c factorize the transformed matrix AJA from 1 to n and compute
c the standard step DN

          CALL TSQRFC(AJA,NR,N,M+N,1,N,IERR,EPSM,AL2NRM,PBAR,ZERO1)
          CALL TSQMUV(AJA,W,FQ,NR,M+N,1,N+1,.FALSE.)
          CALL TSBSLV(AJA,NR,M+N,N,FQ,DN)
          CALL TSPRMV(Y,DN,PBAR,N,0)
          CALL TSPRMV(DN,Y,PIVOT,N,0)
       ENDIF

       IF(IGLOBL.EQ.1) THEN
          IERR = 1
          CALL DSCAL(M+N,-ONE,FQ,1)
       ENDIF

       END

        SUBROUTINE TSSTMX(S,X,NR,N,P,WRK1,WRK2)
        INTEGER NR,N,P
        DOUBLE PRECISION X(P),S(NR,P),WRK1(P),WRK2(P)

C*********************************************************************
C THIS ROUTINE COMPUTES SHAT-TRANS * X, WHERE SHAT IS AN UPSIDE DOWN
C TRIANGULAR MATRIX RESULTING FROM A QL FACTORIZATION OF A MATRIX
C A AND X IS A VECTOR.
C*********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       SHAT  : UPSIDE DOWN TRIANGULAR MATRIX RESULTING FROM A QL
C               FACTORIZATION
C       X     : VECTOR TO BE MULTIPLIED BY SHAT
C       NR    : LEADING DIMENSION OF SHAT
C       N     : ROW DIMENSION OF THE MATRIX A
C       P     : COLUMN DIMENSION OF SHAT
C       WRK   : WORKSPACE
C
C
C       OUTPUT PARAMETERS :
C       -----------------
C
C       WRK2  :  SHAT-TRANS * X
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DCOPY,DDOT
C       TENSOLVE      ...  TSDLOD
C
C*********************************************************************

        INTEGER COL
        DOUBLE PRECISION ZERO
        DOUBLE PRECISION DDOT
        DATA ZERO/0.0D0/

        CALL TSDLOD(P,ZERO,WRK1,1)

        WRK2(1) = S(N+2,1) * X(P)
        IF(P .GT. 1) THEN
           WRK1(P) = S(N,2)
           WRK1(P-1) = S(N+2,2)
           WRK2(2) = DDOT(P,WRK1,1,X,1)
           DO 10 COL = 3, P
              CALL DCOPY(COL-1,S(N-COL+2,COL),1,WRK1(P-COL+2),1)
              WRK1(P-COL+1) = S(N+2,COL)
              WRK2(COL) = DDOT(P,WRK1,1,X,1)
 10        CONTINUE
        ENDIF

        RETURN
        END

      SUBROUTINE TSTRUD(M,N,X,F,G,SC,NWTAKE,STEPMX,STEPTL,DLT,MXTAKE,
     +                  DXN,DFN,FVEC,SCRES,IRETCD,XPLSP,FPLSP,FPREV,
     +                  XPLS,FP,FPLS)

      INTEGER M,N,IRETCD
      DOUBLE PRECISION F,STEPMX,STEPTL,DLT,SCRES,FPLSP,FPLS
      DOUBLE PRECISION X(N),XPLS(N),G(N)
      DOUBLE PRECISION SC(N),XPLSP(N),FPREV(M),FP(M)
      DOUBLE PRECISION DXN(N),DFN(M)
      LOGICAL NWTAKE,MXTAKE
      EXTERNAL FVEC

C***********************************************************************
C THIS ROUTINE DECIDES WHETHER TO ACCEPT XPLS=X+SC AS THE NEXT ITERATE
C AND UPDATES THE TRUST REGION DLT.
C***********************************************************************
C
C
C
C PARAMETERS
C ----------
C M,N          --> DIMENSIONS OF PROBLEM
C X(N)         --> OLD ITERATE X[K-1]
C F            --> 0.50D0 * || FC ||**2
C G(N)         --> GRADIENT AT OLD ITERATE, G(X), OR APPROXIMATE
C SC(N)        --> CURRENT STEP
C NWTAKE       --> BOOLEAN, =.TRUE. IF INPUT STEP TAKEN
C STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE
C STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES
C                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM
C DLT         <--> TRUST REGION RADIUS
C MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED
C DXN         --->DIAGONAL SCALING MATRIX FOR X
C DFN         --->DIAGONAL SCALING MATRIX FOR F
C FVEC        --->SUBROUTINE TO EVALUATE FUNCTION
C
C IRETCD      <--> RETURN CODE
C                    =0 XPLS ACCEPTED AS NEXT ITERATE;
C                       DLT TRUST REGION FOR NEXT ITERATION.
C                    =1 XPLS UNSATISFACTORY BUT ACCEPTED AS NEXT ITERATE
C                       BECAUSE XPLS-X .LT. SMALLEST ALLOWABLE
C                       STEP LENGTH.
C                    =2 F(XPLS) TOO LARGE.  CONTINUE CURRENT ITERATION
C                       WITH NEW REDUCED DLT.
C                    =3 F(XPLS) SUFFICIENTLY SMALL, BUT QUADRATIC MODEL
C                       PREDICTS F(XPLS) SUFFICIENTLY WELL TO CONTINUE
C                       CURRENT ITERATION WITH NEW DOUBLED DLT.
C XPLSP(N)    <--> WORKSPACE [VALUE NEEDS TO BE RETAINED BETWEEN
C                  SUCCESSIVE CALLS OF K-TH GLOBAL STEP]
C FPLSP       <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C FPREV       ---> WORKING VECTOR
C XPLS(N)     <--  NEW ITERATE X[K]
C FP          <--  FUNCTION VALUE AT NEXT ITERATE
C FPLS        <--  FUNCTION VALUE AT NEW ITERATE, F(XPLS)
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DCOPY,DDOT,DNRM2
C       TENSOLVE      ...  TSFSCL
C
C**********************************************************************

      INTEGER I
      DOUBLE PRECISION STEPLN,DLTFP,SLOPE,DLTF,SLP,PQ,RLN,DLTMP
      DOUBLE PRECISION DNRM2,DDOT
      INTRINSIC ABS,MAX
      DOUBLE PRECISION ZERO,TENTH,HALF,ZNN,ONE,TWO
      DATA ZERO,TENTH,HALF,ZNN,ONE,TWO/0.0D0,0.10D0,0.50D0,
     + 0.99D0,1.0D0,2.0D0/

      MXTAKE = .FALSE.
      DO 90 I = 1,N
        XPLS(I) = X(I)+SC(I)
  90  CONTINUE
      STEPLN = DNRM2(N,SC,1)
      CALL TSFSCL(XPLS,DXN,DFN,M,N,FVEC,FP)
      FPLS = HALF*DNRM2(M,FP,1)**2
      DLTF = FPLS-F
      SLOPE = DDOT(N,G,1,SC,1)
      SLP = HALF*SCRES**2-F

c next statement added for case of compilers which do not optimize
c evaluation of next "IF" statement (in which case fplsp could be
c undefined)
c
      IF(IRETCD.EQ.4) FPLSP=0.0
C$      WRITE(*,961) IRETCD,FPLS,FPLSP,DLTF,SLP
      IF(IRETCD.NE.3 .OR. (FPLS.LT.FPLSP .AND. DLTF.LE. 1.E-4*SLP))
     +     GO TO 130
C     IF(IRETCD.EQ.3 .AND. (FPLS.GE.FPLSP .OR. DLTF.GT. 1.E-4*SLP))
C     THEN
C
C       reset XPLS to XPLSP and terminate global step
C
        IRETCD = 0
        CALL DCOPY(N,XPLSP,1,XPLS,1)
        FPLS = FPLSP
        DLT = HALF*DLT
        CALL DCOPY(M,FPREV,1,FP,1)
C$        WRITE(*,951)
        GO TO 230
C     ELSE
C
C       FPLS too large
C
  130   IF(DLTF.LE. 1.E-4*SLP) GO TO 170
C       IF(DLTF.GT. 1.E-4*SLP)
C       THEN
C$          WRITE(*,952)
          PQ = ONE
          RLN = ZERO
          DO 140 I = 1,N
            RLN = MAX(RLN,ABS(SC(I))/MAX(ABS(XPLS(I)),ONE/PQ))
  140     CONTINUE
C$          WRITE(*,962) RLN
          IF(RLN.GE.STEPTL) GO TO 150
C         IF(RLN.LT.STEPTL)
C         THEN
C
C           cannot find satisfactory XPLS sufficiently distinct from X
C
            IRETCD = 1
C$            WRITE(*,954)
            GO TO 230
C         ELSE
C
C           reduce trust region and continue global step
C
  150       IRETCD = 2
            DLTMP = -SLOPE*STEPLN/(TWO*(DLTF-SLOPE))
C$            WRITE(*,963) DLTMP
            IF(DLTMP.GE. TENTH*DLT) GO TO 155
C           IF(DLTMP.LT. TENTH*DLT)
C           THEN
              DLT = TENTH*DLT
              GO TO 160
C           ELSE
  155         IF(DLTMP.GT.HALF*DLT) THEN
                  DLT = HALF*DLT
              ELSE
                  DLT = DLTMP
              ENDIF
C           ENDIF
  160       CONTINUE
C$            WRITE(*,955)
            GO TO 230
C         ENDIF
C       ELSE
C
C         FPLS sufficiently small
C
  170     CONTINUE
C$          WRITE(*,958)
          DLTFP = HALF*SCRES**2-F
C$           WRITE(*,964) DLTFP,NWTAKE
          IF(IRETCD.EQ.2 .OR. ((ABS(DLTFP-DLTF).GT. TENTH*ABS(DLTF))
     +         .AND. (DLTFP.GT.SLOPE)).OR.NWTAKE
     +         .OR. (DLT.GT. ZNN*STEPMX)) GO TO 210
C         IF(IRETCD.NE.2 .AND. (ABS(DLTFP-DLTF) .LE. .1*ABS(DLTF))
C    +         .AND. (.NOT.NWTAKE) .AND. (DLT.LE. .99*STEPMX))
C         THEN
C
C           double trust region and continue global step
C
            IRETCD = 3
            CALL DCOPY(N,XPLS,1,XPLSP,1)
            FPLSP = FPLS
            DLT = MIN(TWO*DLT,STEPMX)
            CALL DCOPY(M,FP,1,FPREV,1)
C$          WRITE(*,959)
            GO TO 230
C         ELSE
C
C           accept XPLS as next iterate.  Choose new trust region.
C
  210       CONTINUE
C$            WRITE(*,960)
            IRETCD = 0
            IF(DLT.GT. ZNN*STEPMX) MXTAKE = .TRUE.
            IF(DLTF.LT. TENTH*DLTFP) GO TO 220
C           IF(DLTF.GE. TENTH*DLTFP)
C           THEN
C
C             Decrease trust region for next iteration
C
              DLT = HALF*DLT
              GO TO 230
C           ELSE
C             Check whether to increase trust region for next iteration
C
  220         IF(DLTF.LE. .75*DLTFP) DLT = MIN(TWO*DLT,STEPMX)
C           ENDIF
C         ENDIF
C       ENDIF
C     ENDIF
  230 CONTINUE
C$      WRITE(*,953)
C$      WRITE(*,956) IRETCD,MXTAKE,DLT,FPLS
C$      WRITE(*,957)
C$      WRITE(*,965) (XPLS(I),I = 1,N)
      RETURN
C
C$  951 FORMAT(' TSTRUD    RESET XPLS TO XPLSP. TERMINATION GLOBAL STEP')
C$  952 FORMAT(' TSTRUD    FPLS TOO LARGE')
C$  953 FORMAT(' TSTRUD    VALUES AFTER CALL TO TREGUP')
C$  954 FORMAT(' TSTRUD    CANNOT FIND SATISFACTORY XPLS DISTINCT FROM',
C$     +       ' X.  TERMINATE GLOBAL STEP')
C$  955 FORMAT(' TSTRUD    REDUCE TRUST REGION. CONTINUE GLOBAL STEP')
C$  956 FORMAT(' TSTRUD       IRETCD=',I3/
C$     +       ' TSTRUD       MXTAKE=',L1/
C$     +       ' TSTRUD       DLT   =',E20.13/
C$     +       ' TSTRUD       FPLS  =',E20.13)
C$  957 FORMAT(' TSTRUD       NEW ITERATE (XPLS)')
C$  958 FORMAT(' TSTRUD    FPLS SUFFICIENTLY SMALL')
C$  959 FORMAT(' TSTRUD    DOUBLE TRUST REGION.  CONTINUE GLOBAL STEP')
C$  960 FORMAT(' TSTRUD    ACCEPT XPLS AS NEW ITERATE.  CHOOSE NEW',
C$     +       ' TRUST REGION.  TERMINATE GLOBAL STEP')
C$  961 FORMAT(' TSTRUD    IRETCD=',I5/
C$     +       ' TSTRUD    FPLS  =',E20.13/
C$     +       ' TSTRUD    FPLSP =',E20.13/
C$     +       ' TSTRUD    DLTF  =',E20.13/
C$     +       ' TSTRUD    SLP   =',E20.13)
C$  962 FORMAT(' TSTRUD    RLN   =',E20.13)
C$  963 FORMAT(' TSTRUD    DLTMP =',E20.13)
C$  964 FORMAT(' TSTRUD    DLTFP =',E20.13/
C$     +       ' TSTRUD    NWTAKE=',L1)
C$  965 FORMAT(' TSTRUD       ',5(E20.13,3X))
      END

      SUBROUTINE TSUDQV(SHAT,WRK1,NR,N,P,CONST1)

      INTEGER NR,N,P
      DOUBLE PRECISION CONST1(P),SHAT(NR,P),WRK1(N)

C**********************************************************************
C THIS ROUTINE COMPUTES SHAT-TRANS * WRK1, WHERE SHAT IS AN UPSIDE
C DOWN TRIANGULAR MATRIX RESULTING FROM A QL FACTORIZATION OF A
C MATRIX A AND WRK1 IS A VECTOR OF LENGTH N.
C**********************************************************************
C
C  INPUT PARAMETERS
C  ----------------
C
C      SHAT : UPSIDE DOWN TRIANGULAR MATRIX RESULTING FROM A QL
C             FACTORIZATION
C      WRK1 : VECTOR TO BE MULTIPLIED BY SHAT
C      NR   : LEADING DIMENSION OF SHAT
C      N    : DIMENSION OF MATRIX A
C      P    : COLUMN DIMENSION OF SHAT
C
C  OUTPUT PARAMETERS
C  -----------------
C
C      CONST1 : SHAT * WRK1
C
C      SUBPROGRAMS CALLED:
C
C      LEVEL 1 BLAS  ...  DDOT
C
C **********************************************************************

      INTEGER J
      DOUBLE PRECISION DDOT

      CONST1(1) = SHAT(N+2,1) * WRK1(N)
      IF(P .GT. 1) THEN
         CONST1(2) = SHAT(N,2) * WRK1(N) + SHAT(N+2,2) * WRK1(N-1)
      ENDIF
      DO 20 J = 3,P
         CONST1(J) = DDOT(J-1,SHAT(N-J+2,J),1,WRK1(N-J+2),1)
     +               + SHAT(N+2,J) * WRK1(N-J+1)
 20   CONTINUE

      RETURN
      END

       SUBROUTINE TSUNSF(F,DF,M)

       INTEGER M
       DOUBLE PRECISION F(M),DF(M)

C*********************************************************************
C THIS ROUTINE UNSCALES A FUNCTION VALUE F.
C*********************************************************************
C
C       INPUT PARAMETERS :
C       ------------------
C
C       DF : DIAGONAL SCALING MATRIX FOR F
C       M  : DIMENSION OF F
C
C       INPUT-OUTPUT PARAMETERS :
C       ------------------
C
C       F  : SCALED FUNCTION VALUE ON ENTRY AND UNSCALED FUNCTION
C            VALUE ON EXIT
C
C**********************************************************************

       INTEGER I

       DO 10 I = 1,M
          F(I) = F(I)/DF(I)
 10    CONTINUE

       RETURN
       END

       SUBROUTINE TSUNSX(X,DX,N)

       INTEGER N
       DOUBLE PRECISION  X(N),DX(N)

C**********************************************************************
C THIS ROUTINE UNSCALES A VECTOR X.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       ------------------
C
C       DX : DIAGONAL SCALING MATRIX FOR X
C       N  : DIMENSION OF X
C
C       OUTPUT PARAMETERS :
C       ------------------
C
C       X  : UNSCALED VECTOR X
C
C**********************************************************************

       INTEGER I

       DO 10 I = 1,N
          X(I) = X(I)/DX(I)
 10    CONTINUE

       RETURN
       END

      SUBROUTINE TSUPSF(FC,XC,XP,SQRN,ITN,MAXM,MAXN,M,N,STEP,S,FV)

      INTEGER MAXM,MAXN,M,N,ITN,SQRN
      DOUBLE PRECISION S(MAXN,*),FV(MAXM,*)
      DOUBLE PRECISION FC(M),XC(N),XP(N),STEP(N)

C**********************************************************************
C THIS ROUTINE UPDATE THE MATRIX S OF PAST DIRECTIONS AND THE MATRIX
C FV OF FUNCTION VALUES.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       ----------------
C
C       FC  : FUNCTION VALUE AT CURRENT ITERATE
C       XC  : CURRENT ITERATE X[K-1]
C       XP  : NEW ITERATE X[K]
C       SQRN: MAXIMUM COLUMN DIMENSION OF S AND FV
C       ITN : ITERATION NUMBER
C       MAXM: LEADING DIMENSION OF FV
C       MAXN: LEADING DIMENSION OF S
C       M   : ROW DIMENSION OF MATRIX FV
C       N   : ROW DIMENSION OF MATRIX S
C       STEP: WORKING VECTOR
C
C
C       INPUT-OUTPUT PARAMETERS :
C        -----------------------
C
C       S   :  MATRIX OF PAST DIRECTIONS (XK - XC)
C       FV  :  MATRIX OF PAST FUNCTIONS VALUES
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DCOPY
C
C**********************************************************************

       INTEGER I,J,L

c update FV

       IF(SQRN.LT.ITN) THEN
          L = SQRN
       ELSE
          L = ITN
       ENDIF
       DO 10 J = L-1,1,-1
          CALL DCOPY(M,FV(1,J),1,FV(1,J+1),1)
 10    CONTINUE

       CALL DCOPY(M,FC,1,FV(1,1),1)

c update S

       DO 30 I = 1,N
          STEP(I) = XC(I)-XP(I)
 30    CONTINUE

       DO 50 J = L-1,1,-1
          DO 40 I = 1,N
             S(I,J+1) = S(I,J) + STEP(I)
 40       CONTINUE
 50    CONTINUE
       CALL DCOPY(N,STEP,1,S(1,1),1)

       RETURN
       END

       SUBROUTINE TSUTMD(AJA,D,NR,M,N,V)

       INTEGER NR,M,N
       DOUBLE PRECISION AJA(NR,N),D(N),V(N)

C**********************************************************************
C THIS ROUTINE MULTIPLIES AN UPPER TRIANGULAR MATRIX (AS STORED IN
C STEWART) BY A VECTOR D.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       AJA  : JACOBIAN AT CURRENT ITERATE
C       D    : VECTOR TO BE MULTIPLIED BY AJA
C       NR   : LEADING DIMENSION OF AJA
C       M,N  : DIMENSIONS OF PROBLEM
C
C       OUTPUT PARAMETERS :
C       -----------------
C
C       V  : VECTOR AJA * D ON EXIT
C
C       SUBPROGRAMS CALLED:
C
C       LEVEL 1 BLAS  ...  DAXPY
C
C**********************************************************************

       INTEGER J

       V(1) = AJA(M+2,1) * D(1) + AJA(1,2) * D(2)
       V(2) = AJA(M+2,2) * D(2)
       DO 20 J = 3, N
          V(J) = AJA(M+2,J) * D(J)
          CALL DAXPY(J-1,D(J),AJA(1,J),1,V,1)
 20    CONTINUE

       RETURN
       END

C****************************** uncmin.f *********************

      SUBROUTINE BAKSLV(NR,N,A,X,B)
C
C PURPOSE
C -------
C SOLVE  AX=B  WHERE A IS UPPER TRIANGULAR MATRIX.
C NOTE THAT A IS INPUT AS A LOWER TRIANGULAR MATRIX AND
C THAT THIS ROUTINE TAKES ITS TRANSPOSE IMPLICITLY.
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C A(N,N)       --> LOWER TRIANGULAR MATRIX (PRESERVED)
C X(N)        <--  SOLUTION VECTOR
C B(N)         --> RIGHT-HAND SIDE VECTOR
C
C NOTE
C ----
C IF B IS NO LONGER REQUIRED BY CALLING ROUTINE,
C THEN VECTORS B AND X MAY SHARE THE SAME STORAGE.
C
      INTEGER NR,N,I,IP1,J
      DOUBLE PRECISION SUM
      DOUBLE PRECISION A(NR,N),X(N),B(N)
C
C SOLVE (L-TRANSPOSE)X=B. (BACK SOLVE)
C
      I=N
      X(I)=B(I)/A(I,I)
      IF(N.EQ.1) RETURN
   30 IP1=I
      I=I-1
      SUM=0.
      DO 40 J=IP1,N
        SUM=SUM+A(J,I)*X(J)
   40 CONTINUE
      X(I)=(B(I)-SUM)/A(I,I)
      IF(I.GT.1) GO TO 30
      RETURN
      END
      SUBROUTINE CHOLDC(NR,N,A,DIAGMX,TOL,ADDMAX)
C
C PURPOSE
C -------
C FIND THE PERTURBED L(L-TRANSPOSE) [WRITTEN LL+] DECOMPOSITION
C OF A+D, WHERE D IS A NON-NEGATIVE DIAGONAL MATRIX ADDED TO A IF
C NECESSARY TO ALLOW THE CHOLESKY DECOMPOSITION TO CONTINUE.
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C A(N,N)      <--> ON ENTRY: MATRIX FOR WHICH TO FIND PERTURBED
C                       CHOLESKY DECOMPOSITION
C                  ON EXIT:  CONTAINS L OF LL+ DECOMPOSITION
C                  IN LOWER TRIANGULAR PART AND DIAGONAL OF "A"
C DIAGMX       --> MAXIMUM DIAGONAL ELEMENT OF "A"
C TOL          --> TOLERANCE
C ADDMAX      <--  MAXIMUM AMOUNT IMPLICITLY ADDED TO DIAGONAL OF "A"
C                  IN FORMING THE CHOLESKY DECOMPOSITION OF A+D
C INTERNAL VARIABLES
C ------------------
C AMINL    SMALLEST ELEMENT ALLOWED ON DIAGONAL OF L
C AMNLSQ   =AMINL**2
C OFFMAX   MAXIMUM OFF-DIAGONAL ELEMENT IN COLUMN OF A
C
C
C DESCRIPTION
C -----------
C THE NORMAL CHOLESKY DECOMPOSITION IS PERFORMED.  HOWEVER, IF AT ANY
C POINT THE ALGORITHM WOULD ATTEMPT TO SET L(I,I)=SQRT(TEMP)
C WITH TEMP < TOL*DIAGMX, THEN L(I,I) IS SET TO SQRT(TOL*DIAGMX)
C INSTEAD.  THIS IS EQUIVALENT TO ADDING TOL*DIAGMX-TEMP TO A(I,I)
C
C
      INTEGER NR,N,J,JM1,K,JP1,I
      DOUBLE PRECISION DIAGMX,TOL,ADDMAX,AMINL,SUM,TEMP,AMNLSQ,OFFMAX
      DOUBLE PRECISION A(NR,N)
C
      ADDMAX=0.
      AMINL=SQRT(DIAGMX*TOL)
      AMNLSQ=AMINL*AMINL
C
C FORM COLUMN J OF L
C
      DO 100 J=1,N
C FIND DIAGONAL ELEMENTS OF L
        SUM=0.
        IF(J.EQ.1) GO TO 20
        JM1=J-1
        DO 10 K=1,JM1
          SUM=SUM + A(J,K)*A(J,K)
   10   CONTINUE
   20   TEMP=A(J,J)-SUM
        IF(TEMP.LT.AMNLSQ) GO TO 30
C       IF(TEMP.GE.AMINL**2)
C       THEN
          A(J,J)=SQRT(TEMP)
          GO TO 40
C       ELSE
C
C FIND MAXIMUM OFF-DIAGONAL ELEMENT IN COLUMN
   30     OFFMAX=0.
          IF(J.EQ.N) GO TO 37
          JP1=J+1
          DO 35 I=JP1,N
            IF(ABS(A(I,J)).GT.OFFMAX) OFFMAX=ABS(A(I,J))
   35     CONTINUE
   37     IF(OFFMAX.LE.AMNLSQ) OFFMAX=AMNLSQ
C
C ADD TO DIAGONAL ELEMENT  TO ALLOW CHOLESKY DECOMPOSITION TO CONTINUE
          A(J,J)=SQRT(OFFMAX)
          ADDMAX=MAX(ADDMAX,OFFMAX-TEMP)
C       ENDIF
C
C FIND I,J ELEMENT OF LOWER TRIANGULAR MATRIX
   40   IF(J.EQ.N) GO TO 100
        JP1=J+1
        DO 70 I=JP1,N
          SUM=0.0
          IF(J.EQ.1) GO TO 60
          JM1=J-1
          DO 50 K=1,JM1
            SUM=SUM+A(I,K)*A(J,K)
   50     CONTINUE
   60     A(I,J)=(A(I,J)-SUM)/A(J,J)
   70   CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE CHLHSN(NR,N,A,EPSM,SX,UDIAG)
C
C PURPOSE
C -------
C FIND THE L(L-TRANSPOSE) [WRITTEN LL+] DECOMPOSITION OF THE PERTURBED
C MODEL HESSIAN MATRIX A+MU*I(WHERE MU\0 AND I IS THE IDENTITY MATRIX)
C WHICH IS SAFELY POSITIVE DEFINITE.  IF A IS SAFELY POSITIVE DEFINITE
C UPON ENTRY, THEN MU=0.
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C A(N,N)      <--> ON ENTRY; "A" IS MODEL HESSIAN (ONLY LOWER
C                  TRIANGULAR PART AND DIAGONAL STORED)
C                  ON EXIT:  A CONTAINS L OF LL+ DECOMPOSITION OF
C                  PERTURBED MODEL HESSIAN IN LOWER TRIANGULAR
C                  PART AND DIAGONAL AND CONTAINS HESSIAN IN UPPER
C                  TRIANGULAR PART AND UDIAG
C EPSM         --> MACHINE EPSILON
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C UDIAG(N)    <--  ON EXIT: CONTAINS DIAGONAL OF HESSIAN
C
C INTERNAL VARIABLES
C ------------------
C TOL              TOLERANCE
C DIAGMN           MINIMUM ELEMENT ON DIAGONAL OF A
C DIAGMX           MAXIMUM ELEMENT ON DIAGONAL OF A
C OFFMAX           MAXIMUM OFF-DIAGONAL ELEMENT OF A
C OFFROW           SUM OF OFF-DIAGONAL ELEMENTS IN A ROW OF A
C EVMIN            MINIMUM EIGENVALUE OF A
C EVMAX            MAXIMUM EIGENVALUE OF A
C
C DESCRIPTION
C -----------
C 1. IF "A" HAS ANY NEGATIVE DIAGONAL ELEMENTS, THEN CHOOSE MU>0
C SUCH THAT THE DIAGONAL OF A:=A+MU*I IS ALL POSITIVE
C WITH THE RATIO OF ITS SMALLEST TO LARGEST ELEMENT ON THE
C ORDER OF SQRT(EPSM).
C
C 2. "A" UNDERGOES A PERTURBED CHOLESKY DECOMPOSITION WHICH
C RESULTS IN AN LL+ DECOMPOSITION OF A+D, WHERE D IS A
C NON-NEGATIVE DIAGONAL MATRIX WHICH IS IMPLICITLY ADDED TO
C "A" DURING THE DECOMPOSITION IF "A" IS NOT POSITIVE DEFINITE.
C "A" IS RETAINED AND NOT CHANGED DURING THIS PROCESS BY
C COPYING L INTO THE UPPER TRIANGULAR PART OF "A" AND THE
C DIAGONAL INTO UDIAG.  THEN THE CHOLESKY DECOMPOSITION ROUTINE
C IS CALLED.  ON RETURN, ADDMAX CONTAINS MAXIMUM ELEMENT OF D.
C
C 3. IF ADDMAX=0, "A" WAS POSITIVE DEFINITE GOING INTO STEP 2
C AND RETURN IS MADE TO CALLING PROGRAM.  OTHERWISE,
C THE MINIMUM NUMBER SDD WHICH MUST BE ADDED TO THE
C DIAGONAL OF A TO MAKE IT SAFELY STRICTLY DIAGONALLY DOMINANT
C IS CALCULATED.  SINCE A+ADDMAX*I AND A+SDD*I ARE SAFELY
C POSITIVE DEFINITE, CHOOSE MU=MIN(ADDMAX,SDD) AND DECOMPOSE
C A+MU*I TO OBTAIN L.
C
      INTEGER NR,N,I,J,IM1,JP1,IP1,JM1
      DOUBLE PRECISION EPSM,TOL,DIAGMX,DIAGMN,POSMAX,AMU,OFFMAX
      DOUBLE PRECISION ADDMAX,EVMIN,EVMAX,OFFROW,SDD,ZERO
      DOUBLE PRECISION A(NR,N),SX(N),UDIAG(N)
      DATA ZERO/0.0D0/
C
C SCALE HESSIAN
C PRE- AND POST- MULTIPLY "A" BY INV(SX)
C
      DO 20 J=1,N
        DO 10 I=J,N
          A(I,J)=A(I,J)/(SX(I)*SX(J))
   10   CONTINUE
   20 CONTINUE
C
C STEP1
C -----
C NOTE:  IF A DIFFERENT TOLERANCE IS DESIRED THROUGHOUT THIS
C ALGORITHM, CHANGE TOLERANCE HERE:
      TOL=SQRT(EPSM)
C
      DIAGMX=A(1,1)
      DIAGMN=A(1,1)
      IF(N.EQ.1) GO TO 35
      DO 30 I=2,N
        IF(A(I,I).LT.DIAGMN) DIAGMN=A(I,I)
        IF(A(I,I).GT.DIAGMX) DIAGMX=A(I,I)
   30 CONTINUE
   35 POSMAX=MAX(DIAGMX,0.0D0)
C
C DIAGMN .LE. 0
C
      IF(DIAGMN.GT.POSMAX*TOL) GO TO 100
C     IF(DIAGMN.LE.POSMAX*TOL)
C     THEN
        AMU=TOL*(POSMAX-DIAGMN)-DIAGMN
        IF(AMU.NE.0.) GO TO 60
C       IF(AMU.EQ.0.)
C       THEN
C
C FIND LARGEST OFF-DIAGONAL ELEMENT OF A
          OFFMAX=0.
          IF(N.EQ.1) GO TO 50
          DO 45 I=2,N
            IM1=I-1
            DO 40 J=1,IM1
              IF(ABS(A(I,J)).GT.OFFMAX) OFFMAX=ABS(A(I,J))
   40       CONTINUE
   45     CONTINUE
   50     AMU=OFFMAX
          IF(AMU.NE.0.) GO TO 55
C         IF(AMU.EQ.0.)
C         THEN
            AMU=1.0
            GO TO 60
C         ELSE
   55       AMU=AMU*(1.0+TOL)
C         ENDIF
C       ENDIF
C
C A=A + MU*I
C
   60   DO 65 I=1,N
          A(I,I)=A(I,I)+AMU
   65   CONTINUE
        DIAGMX=DIAGMX+AMU
C     ENDIF
C
C STEP2
C -----
C COPY LOWER TRIANGULAR PART OF "A" TO UPPER TRIANGULAR PART
C AND DIAGONAL OF "A" TO UDIAG
C
  100 CONTINUE
      DO 110 J=1,N
        UDIAG(J)=A(J,J)
        IF(J.EQ.N) GO TO 110
        JP1=J+1
        DO 105 I=JP1,N
          A(J,I)=A(I,J)
  105   CONTINUE
  110 CONTINUE
C
      CALL CHOLDC(NR,N,A,DIAGMX,TOL,ADDMAX)
C
C
C STEP3
C -----
C IF ADDMAX=0, "A" WAS POSITIVE DEFINITE GOING INTO STEP 2,
C THE LL+ DECOMPOSITION HAS BEEN DONE, AND WE RETURN.
C OTHERWISE, ADDMAX>0.  PERTURB "A" SO THAT IT IS SAFELY
C DIAGONALLY DOMINANT AND FIND LL+ DECOMPOSITION
C
      IF(ADDMAX.LE.0.) GO TO 170
C     IF(ADDMAX.GT.0.)
C     THEN
C
C RESTORE ORIGINAL "A" (LOWER TRIANGULAR PART AND DIAGONAL)
C
        DO 120 J=1,N
          A(J,J)=UDIAG(J)
          IF(J.EQ.N) GO TO 120
          JP1=J+1
          DO 115 I=JP1,N
            A(I,J)=A(J,I)
  115     CONTINUE
  120   CONTINUE
C
C FIND SDD SUCH THAT A+SDD*I IS SAFELY POSITIVE DEFINITE
C NOTE:  EVMIN<0 SINCE A IS NOT POSITIVE DEFINITE;
C
        EVMIN=0.
        EVMAX=A(1,1)
        DO 150 I=1,N
          OFFROW=0.
          IF(I.EQ.1) GO TO 135
          IM1=I-1
          DO 130 J=1,IM1
            OFFROW=OFFROW+ABS(A(I,J))
  130     CONTINUE
  135     IF(I.EQ.N) GO TO 145
          IP1=I+1
          DO 140 J=IP1,N
            OFFROW=OFFROW+ABS(A(J,I))
  140     CONTINUE
  145     EVMIN=MIN(EVMIN,A(I,I)-OFFROW)
          EVMAX=MAX(EVMAX,A(I,I)+OFFROW)
  150   CONTINUE
        SDD=TOL*(EVMAX-EVMIN)-EVMIN
C
C PERTURB "A" AND DECOMPOSE AGAIN
C
        AMU=MIN(SDD,ADDMAX)
        DO 160 I=1,N
          A(I,I)=A(I,I)+AMU
          UDIAG(I)=A(I,I)
  160   CONTINUE
C
C "A" NOW GUARANTEED SAFELY POSITIVE DEFINITE
C
        CALL CHOLDC(NR,N,A,ZERO,TOL,ADDMAX)
C     ENDIF
C
C UNSCALE HESSIAN AND CHOLESKY DECOMPOSITION MATRIX
C
  170 DO 190 J=1,N
        DO 175 I=J,N
          A(I,J)=SX(I)*A(I,J)
  175   CONTINUE
        IF(J.EQ.1) GO TO 185
        JM1=J-1
        DO 180 I=1,JM1
          A(I,J)=SX(I)*SX(J)*A(I,J)
  180   CONTINUE
  185   UDIAG(J)=UDIAG(J)*SX(J)*SX(J)
  190 CONTINUE
      RETURN
      END
      SUBROUTINE DFAUT(N,TYPSIZ,FSCALE,METHOD,IEXP,MSG,NDIGIT,
     +     ITNLIM,IAGFLG,IAHFLG,IPR,DLT,GRADTL,STEPMX,STEPTL)
C
C PURPOSE
C -------
C SET DEFAULT VALUES FOR EACH INPUT VARIABLE TO
C MINIMIZATION ALGORITHM.
C
C PARAMETERS
C ----------
C N            --> DIMENSION OF PROBLEM
C TYPSIZ(N)   <--  TYPICAL SIZE FOR EACH COMPONENT OF X
C FSCALE      <--  ESTIMATE OF SCALE OF MINIMIZATION FUNCTION
C METHOD      <--  ALGORITHM TO USE TO SOLVE MINIMIZATION PROBLEM
C IEXP        <--  =0 IF MINIMIZATION FUNCTION NOT EXPENSIVE TO EVALUATE
C MSG         <--  MESSAGE TO INHIBIT CERTAIN AUTOMATIC CHECKS + OUTPUT
C NDIGIT      <--  NUMBER OF GOOD DIGITS IN MINIMIZATION FUNCTION
C ITNLIM      <--  MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
C IAGFLG      <--  =0 IF ANALYTIC GRADIENT NOT SUPPLIED
C IAHFLG      <--  =0 IF ANALYTIC HESSIAN NOT SUPPLIED
C IPR         <--  DEVICE TO WHICH TO SEND OUTPUT
C DLT         <--  TRUST REGION RADIUS
C GRADTL      <--  TOLERANCE AT WHICH GRADIENT CONSIDERED CLOSE ENOUGH
C                  TO ZERO TO TERMINATE ALGORITHM
C STEPMX      <--  VALUE OF ZERO TO TRIP DEFAULT MAXIMUM IN OPTCHK
C STEPTL      <--  TOLERANCE AT WHICH SUCCESSIVE ITERATES CONSIDERED
C                  CLOSE ENOUGH TO TERMINATE ALGORITHM
C
      INTEGER N,METHOD,IEXP,MSG,NDIGIT,ITNLIM,IAGFLG,IAHFLG,IPR,I
      DOUBLE PRECISION FSCALE,DLT,GRADTL,STEPMX,STEPTL,EPSM,DPMEPS
      DOUBLE PRECISION TYPSIZ(N),ZERO,ONE,THREE
      DATA ZERO,ONE,THREE/0.0D0,1.0D0,3.0D0/
C
C SET TYPICAL SIZE OF X AND MINIMIZATION FUNCTION
      DO 10 I=1,N
        TYPSIZ(I)=ONE
   10 CONTINUE
      FSCALE=ONE
C
C SET TOLERANCES
      DLT=-ONE
      EPSM=DPMEPS()
      GRADTL=EPSM**(ONE/THREE)
      STEPMX=ZERO
      STEPTL=SQRT(EPSM)
C
C SET FLAGS
      METHOD=1
      IEXP=1
      MSG=0
      NDIGIT=-1
      ITNLIM=150
      IAGFLG=0
      IAHFLG=0
      IPR=6
C
      RETURN
      END
      SUBROUTINE DOGDRV(NR,N,X,F,G,A,P,XPLS,FPLS,FCN,SX,STEPMX,
     + STEPTL,DLT,IRETCD,MXTAKE,SC,WRK1,WRK2,WRK3,IPR,
     + AJA,ANLS,SHAT,VECT1,VECT2,VECT3,VECT4,VECT5,VECT6,
     + NRM,NRN,MM,NN,IQ)
C
C PURPOSE
C -------
C FIND A NEXT NEWTON ITERATE (XPLS) BY THE DOUBLE DOGLEG METHOD
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C X(N)         --> OLD ITERATE X[K-1]
C F            --> FUNCTION VALUE AT OLD ITERATE, F(X)
C G(N)         --> GRADIENT  AT OLD ITERATE, G(X), OR APPROXIMATE
C A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN
C                  IN LOWER TRIANGULAR PART AND DIAGONAL
C P(N)         --> NEWTON STEP
C XPLS(N)     <--  NEW ITERATE X[K]
C FPLS        <--  FUNCTION VALUE AT NEW ITERATE, F(XPLS)
C FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE
C STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES
C                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM
C DLT         <--> TRUST REGION RADIUS
C                  [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C IRETCD      <--  RETURN CODE
C                    =0 SATISFACTORY XPLS FOUND
C                    =1 FAILED TO FIND SATISFACTORY XPLS SUFFICIENTLY
C                       DISTINCT FROM X
C MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED
C SC(N)        --> WORKSPACE [CURRENT STEP]
C WRK1(N)      --> WORKSPACE (AND PLACE HOLDING ARGUMENT TO TREGUP)
C WRK2(N)      --> WORKSPACE
C WRK3(N)      --> WORKSPACE
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C
      INTEGER N,IPR,NRM,NRN,MM,NN,IQ,I,NR,IRETCD
      DOUBLE PRECISION F,FPLS,STEPMX,STEPTL,DLT,TMP,RNWTLN,CLN
      DOUBLE PRECISION ETA,FPLSP
      DOUBLE PRECISION X(N),XPLS(N),G(N),P(N)
      DOUBLE PRECISION SX(N)
      DOUBLE PRECISION SC(N),WRK1(N),WRK2(N),WRK3(N)
      DOUBLE PRECISION A(NR,N)
      DOUBLE PRECISION AJA(NRM,NN),ANLS(NRM,N),SHAT(NRN,N)
      DOUBLE PRECISION VECT1(MM),VECT2(MM),VECT3(MM)
      DOUBLE PRECISION VECT4(MM),VECT5(MM),VECT6(MM)
      LOGICAL FSTDOG,NWTAKE,MXTAKE
      EXTERNAL FCN
C
      IRETCD=4
      FSTDOG=.TRUE.
      TMP=0.
      DO 5 I=1,N
        TMP=TMP+(SX(I)*P(I))**2
    5 CONTINUE
      RNWTLN=SQRT(TMP)
C
  100 CONTINUE
C
C FIND NEW STEP BY DOUBLE DOGLEG ALGORITHM
      CALL DOGSTP(NR,N,G,A,P,SX,RNWTLN,DLT,NWTAKE,FSTDOG,
     +     WRK1,WRK2,CLN,ETA,SC,IPR,STEPMX)
C
C CHECK NEW POINT AND UPDATE TRUST REGION
      CALL TREGUP(NR,N,X,F,G,A,FCN,SC,SX,NWTAKE,STEPMX,STEPTL,DLT,
     + IRETCD,WRK3,FPLSP,XPLS,FPLS,MXTAKE,IPR,2,WRK1,
     + AJA,ANLS,SHAT,VECT1,VECT2,VECT3,VECT4,VECT5,VECT6,
     + NRM,NRN,MM,NN,IQ)
      IF(IRETCD.LE.1) RETURN
      GO TO 100
      END
      SUBROUTINE DOGSTP(NR,N,G,A,P,SX,RNWTLN,DLT,NWTAKE,FSTDOG,
     +     SSD,V,CLN,ETA,SC,IPR,STEPMX)
C
C PURPOSE
C -------
C FIND NEW STEP BY DOUBLE DOGLEG ALGORITHM
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C G(N)         --> GRADIENT AT CURRENT ITERATE, G(X)
C A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN IN
C                  LOWER PART AND DIAGONAL
C P(N)         --> NEWTON STEP
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C RNWTLN       --> NEWTON STEP LENGTH
C DLT         <--> TRUST REGION RADIUS
C NWTAKE      <--> BOOLEAN, =.TRUE. IF NEWTON STEP TAKEN
C FSTDOG      <--> BOOLEAN, =.TRUE. IF ON FIRST LEG OF DOGLEG
C SSD(N)      <--> WORKSPACE [CAUCHY STEP TO THE MINIMUM OF THE
C                  QUADRATIC MODEL IN THE SCALED STEEPEST DESCENT
C                  DIRECTION] [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C V(N)        <--> WORKSPACE  [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C CLN         <--> CAUCHY LENGTH
C                  [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C ETA              [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C SC(N)       <--  CURRENT STEP
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE
C
C INTERNAL VARIABLES
C ------------------
C CLN              LENGTH OF CAUCHY STEP
C
      INTEGER NR,N,IPR,I,J
      DOUBLE PRECISION STEPMX,DLT,RNWTLN,CLN,ETA,ALPHA,BETA,TMP
      DOUBLE PRECISION DOT1,DOT2,ALAM
      DOUBLE PRECISION G(N),P(N)
      DOUBLE PRECISION SX(N)
      DOUBLE PRECISION SC(N),SSD(N),V(N)
      DOUBLE PRECISION A(NR,N)
      DOUBLE PRECISION DDOT
      LOGICAL NWTAKE,FSTDOG
C
C CAN WE TAKE NEWTON STEP
C
      IF(RNWTLN.GT.DLT) GO TO 100
C     IF(RNWTLN.LE.DLT)
C     THEN
        NWTAKE=.TRUE.
        DO 10 I=1,N
          SC(I)=P(I)
   10   CONTINUE
        DLT=RNWTLN
        GO TO 700
C     ELSE
C
C NEWTON STEP TOO LONG
C CAUCHY STEP IS ON DOUBLE DOGLEG CURVE
C
  100   NWTAKE=.FALSE.
        IF(.NOT.FSTDOG) GO TO 200
C       IF(FSTDOG)
C       THEN
C
C         CALCULATE DOUBLE DOGLEG CURVE (SSD)
          FSTDOG=.FALSE.
          ALPHA=0.
          DO 105 I = 1, N
            SSD(I)  = G(I)/SX(I)
  105     CONTINUE
          DO 110 I=1,N
            ALPHA=ALPHA + SSD(I)*SSD(I)
  110     CONTINUE
          BETA=0.
          DO 130 I=1,N
            TMP=0.
            DO 120 J=I,N
              TMP=TMP + (A(J,I)/SX(J))*SSD(J)
  120       CONTINUE
            BETA=BETA+TMP*TMP
  130     CONTINUE
          DO 140 I=1,N
            SSD(I)=-(ALPHA/BETA)*SSD(I)
  140     CONTINUE
          CLN=ALPHA*SQRT(ALPHA)/BETA
          ETA=.2 + (.8*ALPHA*ALPHA)/(-BETA*DDOT(N,G,1,P,1))
          DO 150 I=1,N
            V(I)=ETA*SX(I)*P(I) - SSD(I)
  150     CONTINUE
          IF (DLT .EQ. (-1.0)) DLT = MIN(CLN, STEPMX)
C       ENDIF
  200   IF(ETA*RNWTLN.GT.DLT) GO TO 220
C       IF(ETA*RNWTLN .LE. DLT)
C       THEN
C
C         TAKE PARTIAL STEP IN NEWTON DIRECTION
C
          DO 210 I=1,N
            SC(I)=(DLT/RNWTLN)*P(I)
  210     CONTINUE
          GO TO 700
C       ELSE
  220     IF(CLN.LT.DLT) GO TO 240
C         IF(CLN.GE.DLT)
C         THEN
C           TAKE STEP IN STEEPEST DESCENT DIRECTION
C
            DO 230 I=1,N
              SC(I)=(DLT/CLN)*SSD(I)/SX(I)
  230       CONTINUE
            GO TO 700
C         ELSE
C
C           CALCULATE CONVEX COMBINATION OF SSD AND ETA*P
C           WHICH HAS SCALED LENGTH DLT
C
  240       DOT1=DDOT(N,V,1,SSD,1)
            DOT2=DDOT(N,V,1,V,1)
            ALAM=(-DOT1+SQRT((DOT1*DOT1)-DOT2*(CLN*CLN-DLT*DLT)))/DOT2
            DO 250 I=1,N
              SC(I)=(SSD(I) + ALAM*V(I))/SX(I)
  250       CONTINUE
C         ENDIF
C       ENDIF
C     ENDIF
  700 CONTINUE
      RETURN
      END
      SUBROUTINE D1FCN(N,X,G,AJA,ANLS,SHAT,VECT1,VECT2,VECT3,
     +                 VECT4,VECT5,VECT6,NRM,NRN,MM,NN,IQ)
C
C PURPOSE
C -------
C DUMMY ROUTINE TO PREVENT UNSATISFIED EXTERNAL DIAGNOSTIC
C WHEN SPECIFIC ANALYTIC GRADIENT FUNCTION NOT SUPPLIED.
C
      INTEGER N,NRM,NRN,MM,NN,IQ
      DOUBLE PRECISION X(N),G(N)
      DOUBLE PRECISION AJA(NRM,NN),ANLS(NRM,N),SHAT(NRN,N)
      DOUBLE PRECISION VECT1(MM),VECT2(MM),VECT3(MM)
      DOUBLE PRECISION VECT4(MM),VECT5(MM),VECT6(MM)
      STOP
      END
      SUBROUTINE D2FCN(NR,N,X,H)
C
C PURPOSE
C -------
C DUMMY ROUTINE TO PREVENT UNSATISFIED EXTERNAL DIAGNOSTIC
C WHEN SPECIFIC ANALYTIC HESSIAN FUNCTION NOT SUPPLIED.
C
      INTEGER NR,N
      DOUBLE PRECISION X(N),H(NR,N)
      STOP
      END
       SUBROUTINE FORSLV (NR,N,A,X,B)
C
C PURPOSE
C --------
C SOLVE  AX=B  WHERE A  IS LOWER TRIANGULAR  MATRIX
C
C PARAMETERS
C ---------
C
C NR            -----> ROW DIMENSION OF MATRIX
C N             -----> DIMENSION OF PROBLEM
C A(N,N)        -----> LOWER TRIANGULAR MATRIX (PRESERVED)
C X(N)          <----  SOLUTION VECTOR
C B(N)           ----> RIGHT-HAND SIDE VECTOR
C
C NOTE
C-----
C THEN VECTORS B AND X MAY SHARE THE SAME STORAGE
C
      INTEGER NR,N,I,J,IM1
      DOUBLE PRECISION SUM
      DOUBLE PRECISION A(NR,N),X(N),B(N)
C
C SOLVE LX=B.  (FORWARD  SOLVE)
C
      X(1)=B(1)/A(1,1)
      DO 20 I=2,N
        SUM=0.0
        IM1=I-1
        DO 10 J=1,IM1
          SUM=SUM+A(I,J)*X(J)
   10   CONTINUE
        X(I)=(B(I)-SUM)/A(I,I)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE FSTOCD (N, X, FCN, SX, RNOISE, G)
C PURPOSE
C -------
C FIND CENTRAL DIFFERENCE APPROXIMATION G TO THE FIRST DERIVATIVE
C (GRADIENT) OF THE FUNCTION DEFINED BY FCN AT THE POINT X.
C
C PARAMETERS
C ----------
C N            --> DIMENSION OF PROBLEM
C X            --> POINT AT WHICH GRADIENT IS TO BE APPROXIMATED.
C FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION.
C SX           --> DIAGONAL SCALING MATRIX FOR X.
C RNOISE       --> RELATIVE NOISE IN FCN [F(X)].
C G           <--  CENTRAL DIFFERENCE APPROXIMATION TO GRADIENT.
C
C
      INTEGER N,I
      DOUBLE PRECISION RNOISE,THIRD,XTEMPI,FPLUS,FMINUS,STEPI
      DOUBLE PRECISION X(N)
      DOUBLE PRECISION SX(N)
      DOUBLE PRECISION G(N)
      EXTERNAL FCN
C
C FIND I TH  STEPSIZE, EVALUATE TWO NEIGHBORS IN DIRECTION OF I TH
C UNIT VECTOR, AND EVALUATE I TH  COMPONENT OF GRADIENT.
C
      THIRD = 1.0/3.0
      DO 10 I = 1, N
         STEPI = RNOISE**THIRD * MAX(ABS(X(I)), 1.0/SX(I))
         XTEMPI = X(I)
         X(I) = XTEMPI + STEPI
         CALL FCN (N, X, FPLUS)
         X(I) = XTEMPI - STEPI
         CALL FCN (N, X, FMINUS)
         X(I) = XTEMPI
         G(I) = (FPLUS - FMINUS)/(2.0*STEPI)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FSTOFD(NR,M,N,XPLS,FCN,FPLS,A,SX,RNOISE,FHAT,ICASE,
     + AJA,ANLS,SHAT,VECT1,VECT2,VECT3,VECT4,VECT5,VECT6,NRM,NRN,
     + MM,NN,IQ)
C PURPOSE
C -------
C FIND FIRST ORDER FORWARD FINITE DIFFERENCE APPROXIMATION "A" TO THE
C FIRST DERIVATIVE OF THE FUNCTION DEFINED BY THE SUBPROGRAM "FNAME"
C EVALUATED AT THE NEW ITERATE "XPLS".
C
C
C FOR OPTIMIZATION USE THIS ROUTINE TO ESTIMATE:
C 1) THE FIRST DERIVATIVE (GRADIENT) OF THE OPTIMIZATION FUNCTION "FCN
C    ANALYTIC USER ROUTINE HAS BEEN SUPPLIED;
C 2) THE SECOND DERIVATIVE (HESSIAN) OF THE OPTIMIZATION FUNCTION
C    IF NO ANALYTIC USER ROUTINE HAS BEEN SUPPLIED FOR THE HESSIAN BUT
C    ONE HAS BEEN SUPPLIED FOR THE GRADIENT ("FCN") AND IF THE
C    OPTIMIZATION FUNCTION IS INEXPENSIVE TO EVALUATE
C
C NOTE
C ----
C _M=1 (OPTIMIZATION) ALGORITHM ESTIMATES THE GRADIENT OF THE FUNCTION
C      (FCN).   FCN(X) # F: R(N)-->R(1)
C _M=N (SYSTEMS) ALGORITHM ESTIMATES THE JACOBIAN OF THE FUNCTION
C      FCN(X) # F: R(N)-->R(N).
C _M=N (OPTIMIZATION) ALGORITHM ESTIMATES THE HESSIAN OF THE OPTIMIZATIO
C      FUNCTION, WHERE THE HESSIAN IS THE FIRST DERIVATIVE OF "FCN"
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C M            --> NUMBER OF ROWS IN A
C N            --> NUMBER OF COLUMNS IN A; DIMENSION OF PROBLEM
C XPLS(N)      --> NEW ITERATE:  X[K]
C FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION
C FPLS(M)      --> _M=1 (OPTIMIZATION) FUNCTION VALUE AT NEW ITERATE:
C                       FCN(XPLS)
C                  _M=N (OPTIMIZATION) VALUE OF FIRST DERIVATIVE
C                       (GRADIENT) GIVEN BY USER FUNCTION FCN
C                  _M=N (SYSTEMS)  FUNCTION VALUE OF ASSOCIATED
C                       MINIMIZATION FUNCTION
C A(NR,N)     <--  FINITE DIFFERENCE APPROXIMATION (SEE NOTE).  ONLY
C                  LOWER TRIANGULAR MATRIX AND DIAGONAL ARE RETURNED
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C RNOISE       --> RELATIVE NOISE IN FCN [F(X)]
C FHAT(M)      --> WORKSPACE
C ICASE        --> =1 OPTIMIZATION (GRADIENT)
C                  =2 SYSTEMS
C                  =3 OPTIMIZATION (HESSIAN)
C
C INTERNAL VARIABLES
C ------------------
C STEPSZ - STEPSIZE IN THE J-TH VARIABLE DIRECTION
C
      INTEGER NR,M,N,ICASE,NRM,NRN,MM,NN,IQ,I,J,NM1,JP1
      DOUBLE PRECISION RNOISE,STEPSZ,XTMPJ,SQRTR,RSTEPSZ,HALF,ONE
      DOUBLE PRECISION XPLS(N),FPLS(M)
      DOUBLE PRECISION FHAT(M)
      DOUBLE PRECISION SX(N)
      DOUBLE PRECISION A(NR,N)
      DOUBLE PRECISION AJA(NRM,NN),ANLS(NRM,N),SHAT(NRN,N)
      DOUBLE PRECISION VECT1(MM),VECT2(MM),VECT3(MM)
      DOUBLE PRECISION VECT4(MM),VECT5(MM),VECT6(MM)
      DATA HALF,ONE/0.50D0,1.0D0/
      EXTERNAL FCN
C
C FIND J-TH COLUMN OF A
C EACH COLUMN IS DERIVATIVE OF F(FCN) WITH RESPECT TO XPLS(J)
C
      SQRTR = SQRT(RNOISE)
      DO 30 J=1,N
        STEPSZ=SQRTR*MAX(ABS(XPLS(J)),ONE/SX(J))
        XTMPJ=XPLS(J)
        XPLS(J)=XTMPJ+STEPSZ
        CALL FCN(N,XPLS,FHAT,AJA,ANLS,SHAT,VECT1,VECT2,VECT3,
     +  VECT4,VECT5,VECT6,NRM,NRN,MM,NN,IQ)
        XPLS(J)=XTMPJ
        RSTEPSZ = ONE/STEPSZ
        DO 20 I=1,M
          A(I,J)=(FHAT(I)-FPLS(I))*RSTEPSZ
   20   CONTINUE
   30 CONTINUE
      IF(ICASE.NE.3) RETURN
C
C IF COMPUTING HESSIAN, A MUST BE SYMMETRIC
C
      IF(N.EQ.1) RETURN
      NM1=N-1
      DO 50 J=1,NM1
        JP1=J+1
        DO 40 I=JP1,M
          A(I,J)=(A(I,J)+A(J,I))*HALF
   40   CONTINUE
   50 CONTINUE
      RETURN
      END
      SUBROUTINE HOOKDR(NR,N,X,F,G,A,UDIAG,P,XPLS,FPLS,FCN,SX,STEPMX,
     +     STEPTL,DLT,IRETCD,MXTAKE,AMU,DLTP,PHI,PHIP0,
     +     SC,XPLSP,WRK0,EPSM,ITNCNT,IPR)
C
C PURPOSE
C -------
C FIND A NEXT NEWTON ITERATE (XPLS) BY THE MORE-HEBDON METHOD
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C X(N)         --> OLD ITERATE X[K-1]
C F            --> FUNCTION VALUE AT OLD ITERATE, F(X)
C G(N)         --> GRADIENT AT OLD ITERATE, G(X), OR APPROXIMATE
C A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN IN LOWER
C                  TRIANGULAR PART AND DIAGONAL.
C                  HESSIAN IN UPPER TRIANGULAR PART AND UDIAG.
C UDIAG(N)     --> DIAGONAL OF HESSIAN IN A(.,.)
C P(N)         --> NEWTON STEP
C XPLS(N)     <--  NEW ITERATE X[K]
C FPLS        <--  FUNCTION VALUE AT NEW ITERATE, F(XPLS)
C FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE
C STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES
C                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM
C DLT         <--> TRUST REGION RADIUS
C IRETCD      <--  RETURN CODE
C                    =0 SATISFACTORY XPLS FOUND
C                    =1 FAILED TO FIND SATISFACTORY XPLS SUFFICIENTLY
C                       DISTINCT FROM X
C MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED
C AMU         <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C DLTP        <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C PHI         <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C PHIP0       <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C SC(N)        --> WORKSPACE
C XPLSP(N)     --> WORKSPACE
C WRK0(N)      --> WORKSPACE
C EPSM         --> MACHINE EPSILON
C ITNCNT       --> ITERATION COUNT
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C
      INTEGER NR,N,IRETCD,ITNCNT,IPR,I,J
      DOUBLE PRECISION FPLS,STEPMX,STEPTL,DLT,AMU,DLTP,PHI
      DOUBLE PRECISION PHIP0,EPSM,TMP,RNWTLN,ALPHA,BETA,F,ZERO,ONE
      DOUBLE PRECISION X(N),G(N),P(N),XPLS(N),SX(N)
      DOUBLE PRECISION A(NR,N),UDIAG(N)
      DOUBLE PRECISION SC(N),XPLSP(N),WRK0(N)
      LOGICAL MXTAKE,NWTAKE
      LOGICAL FSTIME
      DATA ZERO,ONE/0.0D0,1.0D0/
      EXTERNAL FCN
C
      IRETCD=4
      FSTIME=.TRUE.
      TMP=ZERO
      DO 5 I=1,N
        TMP=TMP+(SX(I)*P(I))**2
    5 CONTINUE
      RNWTLN=SQRT(TMP)
      IF(ITNCNT.EQ.1) THEN
        AMU=ZERO
C
C       IF FIRST ITERATION AND TRUST REGION NOT PROVIDED BY USER,
C       COMPUTE INITIAL TRUST REGION.
C
        IF(DLT.EQ. -ONE) THEN
          ALPHA=ZERO
          DO 10 I=1,N
            ALPHA=ALPHA+(G(I)/SX(I))**2
   10     CONTINUE
          BETA=ZERO
          DO 30 I=1,N
            TMP=ZERO
            DO 20 J=I,N
              TMP=TMP + (A(J,I)*G(J))/(SX(J)*SX(J))
   20       CONTINUE
            BETA=BETA+TMP*TMP
   30     CONTINUE
          DLT=ALPHA*SQRT(ALPHA)/BETA
          DLT = MIN(DLT, STEPMX)
        ENDIF
      ENDIF
C
  100 CONTINUE
C
C FIND NEW STEP BY MORE-HEBDON ALGORITHM
      CALL HOOKST(NR,N,G,A,UDIAG,P,SX,RNWTLN,DLT,AMU,
     +     DLTP,PHI,PHIP0,FSTIME,SC,NWTAKE,WRK0,EPSM,IPR)
      DLTP=DLT
C
C CHECK NEW POINT AND UPDATE TRUST REGION
C     CALL TREGUP(NR,N,X,F,G,A,FCN,SC,SX,NWTAKE,STEPMX,STEPTL,
C    +         DLT,IRETCD,XPLSP,FPLSP,XPLS,FPLS,MXTAKE,IPR,3,UDIAG)
      IF(IRETCD.LE.1) RETURN
      GO TO 100
      END
      SUBROUTINE HOOKST(NR,N,G,A,UDIAG,P,SX,RNWTLN,DLT,AMU,
     +     DLTP,PHI,PHIP0,FSTIME,SC,NWTAKE,WRK0,EPSM,IPR)
C
C PURPOSE
C -------
C FIND NEW STEP BY MORE-HEBDON ALGORITHM
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C G(N)         --> GRADIENT AT CURRENT ITERATE, G(X)
C A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN IN
C                  LOWER TRIANGULAR PART AND DIAGONAL.
C                  HESSIAN OR APPROX IN UPPER TRIANGULAR PART
C UDIAG(N)     --> DIAGONAL OF HESSIAN IN A(.,.)
C P(N)         --> NEWTON STEP
C SX(N)        --> DIAGONAL SCALING MATRIX FOR N
C RNWTLN       --> NEWTON STEP LENGTH
C DLT         <--> TRUST REGION RADIUS
C AMU         <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C DLTP         --> TRUST REGION RADIUS AT LAST EXIT FROM THIS ROUTINE
C PHI         <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C PHIP0       <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C FSTIME      <--> BOOLEAN. =.TRUE. IF FIRST ENTRY TO THIS ROUTINE
C                  DURING K-TH ITERATION
C SC(N)       <--  CURRENT STEP
C NWTAKE      <--  BOOLEAN, =.TRUE. IF NEWTON STEP TAKEN
C WRK0         --> WORKSPACE
C EPSM         --> MACHINE EPSILON
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C
      INTEGER NR,N,IPR,I,JP1,J
      DOUBLE PRECISION RNWTLN,DLT,AMU,DLTP,PHI,PHIP0,HI,ALO,PHIP
      DOUBLE PRECISION AMUUP,ADDMAX,STEPLN,AMULO,EPSM,ZERO
      DOUBLE PRECISION G(N),P(N),SX(N),SC(N),WRK0(N)
      DOUBLE PRECISION A(NR,N),UDIAG(N)
      DOUBLE PRECISION DNRM2
      LOGICAL NWTAKE,DONE
      LOGICAL FSTIME
      DATA ZERO/0.0D0/
C
C HI AND ALO ARE CONSTANTS USED IN THIS ROUTINE.
C CHANGE HERE IF OTHER VALUES ARE TO BE SUBSTITUTED.
      HI=1.5
      ALO=.75
C -----
      IF(RNWTLN.LE.HI*DLT) THEN
C
C       TAKE NEWTON STEP
C
        NWTAKE=.TRUE.
        DO 10 I=1,N
          SC(I)=P(I)
   10   CONTINUE
        DLT=MIN(DLT,RNWTLN)
        AMU=0.
        RETURN
      ENDIF
C
C     NEWTON STEP NOT TAKEN
C
      NWTAKE=.FALSE.
C
C SET PHIP TO 1.0 FOR COMPILATION. THIS SUBROUTINE IS NOT CURRENTLY
C USED BY TENSOLVE.
C
      PHIP = 1.0
      IF(AMU.GT.ZERO) THEN
        AMU=AMU - (PHI+DLTP) *((DLTP-DLT)+PHI)/(DLT*PHIP)
      ENDIF
      PHI=RNWTLN-DLT
      IF(FSTIME) THEN
        DO 25 I=1,N
          WRK0(I)=SX(I)*SX(I)*P(I)
   25   CONTINUE
C
C         SOLVE L*Y = (SX**2)*P
C
        CALL FORSLV(NR,N,A,WRK0,WRK0)
        PHIP0=-DNRM2(N,WRK0,1)**2/RNWTLN
        FSTIME=.FALSE.
      ENDIF
      PHIP=PHIP0
      AMULO=-PHI/PHIP
      AMUUP=0.0
      DO 30 I=1,N
        AMUUP=AMUUP+(G(I)*G(I))/(SX(I)*SX(I))
   30 CONTINUE
      AMUUP=SQRT(AMUUP)/DLT
      DONE=.FALSE.
C
C       TEST VALUE OF AMU; GENERATE NEXT AMU IF NECESSARY
C
  100 CONTINUE
      IF(DONE) RETURN
      IF(AMU.GE.AMULO .AND. AMU.LE.AMUUP) GO TO 110
C     IF(AMU.LT.AMULO .OR.  AMU.GT.AMUUP)
C     THEN
        AMU=MAX(SQRT(AMULO*AMUUP),AMUUP*1.0E-3)
C     ENDIF
  110 CONTINUE
C
C     COPY (H,UDIAG) TO L
C     WHERE H <-- H+AMU*(SX**2) [DO NOT ACTUALLY CHANGE (H,UDIAG)]
      DO 130 J=1,N
        A(J,J)=UDIAG(J) + AMU*SX(J)*SX(J)
        IF(J.EQ.N) GO TO 130
        JP1=J+1
        DO 120 I=JP1,N
          A(I,J)=A(J,I)
  120   CONTINUE
  130 CONTINUE
C
C     FACTOR H=L(L+)
C
      CALL CHOLDC(NR,N,A,ZERO,SQRT(EPSM),ADDMAX)
C
C     SOLVE H*P = L(L+)*SC = -G
C
      DO 140 I=1,N
        WRK0(I)=-G(I)
  140 CONTINUE
      CALL LLTSLV(NR,N,A,SC,WRK0)
C
C     RESET H.  NOTE SINCE UDIAG HAS NOT BEEN DESTROYED WE NEED DO
C     NOTHING HERE.  H IS IN THE UPPER PART AND IN UDIAG, STILL INTACT
C
      STEPLN=0.
      DO 150 I=1,N
        STEPLN=STEPLN + SX(I)*SX(I)*SC(I)*SC(I)
  150 CONTINUE
      STEPLN=SQRT(STEPLN)
      PHI=STEPLN-DLT
      DO 160 I=1,N
        WRK0(I)=SX(I)*SX(I)*SC(I)
  160 CONTINUE
      CALL FORSLV(NR,N,A,WRK0,WRK0)
      PHIP=-DNRM2(N,WRK0,1)**2/STEPLN
      IF((ALO*DLT.GT.STEPLN .OR. STEPLN.GT.HI*DLT) .AND.
     +     (AMUUP-AMULO.GT.0.)) GO TO 170
C     IF((ALO*DLT.LE.STEPLN .AND. STEPLN.LE.HI*DLT) .OR.
C          (AMUUP-AMULO.LE.0.))
C     THEN
C
C       SC IS ACCEPTABLE HOOKSTEP
C
        DONE=.TRUE.
        GO TO 100
C     ELSE
C
C       SC NOT ACCEPTABLE HOOKSTEP.  SELECT NEW AMU
C
  170   CONTINUE
        AMULO=MAX(AMULO,AMU-(PHI/PHIP))
        IF(PHI.LT.0.) AMUUP=MIN(AMUUP,AMU)
        AMU=AMU-(STEPLN*PHI)/(DLT*PHIP)
        GO TO 100
C      ENDIF
C     ENDIF
      END
      SUBROUTINE HSNINT(NR,N,A,SX,METHOD)
C
C
C PURPOSE
C -------
C PROVIDE INITIAL HESSIAN WHEN USING SECANT UPDATES
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C A(N,N)      <--  INITIAL HESSIAN (LOWER TRIANGULAR MATRIX)
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C METHOD       --> ALGORITHM TO USE TO SOLVE MINIMIZATION PROBLEM
C                    =1,2 FACTORED SECANT METHOD USED
C                    =3   UNFACTORED SECANT METHOD USED
C
      INTEGER NR,N,METHOD,J,I,JP1
      DOUBLE PRECISION A(NR,N),SX(N)
C
      DO 100 J=1,N
        IF(METHOD.EQ.3) A(J,J)=SX(J)*SX(J)
        IF(METHOD.NE.3) A(J,J)=SX(J)
        IF(J.EQ.N) GO TO 100
        JP1=J+1
        DO 90 I=JP1,N
          A(I,J)=0.
   90   CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE LLTSLV(NR,N,A,X,B)
C
C PURPOSE
C -------
C SOLVE AX=B WHERE A HAS THE FORM L(L-TRANSPOSE)
C BUT ONLY THE LOWER TRIANGULAR PART, L, IS STORED.
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C A(N,N)       --> MATRIX OF FORM L(L-TRANSPOSE).
C                  ON RETURN A IS UNCHANGED.
C X(N)        <--  SOLUTION VECTOR
C B(N)         --> RIGHT-HAND SIDE VECTOR
C
C NOTE
C ----
C IF B IS NOT REQUIRED BY CALLING PROGRAM, THEN
C B AND X MAY SHARE THE SAME STORAGE.
C
      INTEGER NR,N
      DOUBLE PRECISION A(NR,N),X(N),B(N)
C
C FORWARD SOLVE, RESULT IN X
C
      CALL FORSLV(NR,N,A,X,B)
C
C BACK SOLVE, RESULT IN X
C
      CALL BAKSLV(NR,N,A,X,X)
      RETURN
      END
      SUBROUTINE OPTCHK(N,X,TYPSIZ,SX,FSCALE,GRADTL,ITNLIM,NDIGIT,EPSM,
     +     DLT,METHOD,IEXP,IAGFLG,IAHFLG,STEPMX,MSG,IPR)
C
C PURPOSE
C -------
C CHECK INPUT FOR REASONABLENESS
C
C PARAMETERS
C ----------
C N            --> DIMENSION OF PROBLEM
C X(N)         --> ON ENTRY, ESTIMATE TO ROOT OF FCN
C TYPSIZ(N)   <--> TYPICAL SIZE OF EACH COMPONENT OF X
C SX(N)       <--  DIAGONAL SCALING MATRIX FOR X
C FSCALE      <--> ESTIMATE OF SCALE OF OBJECTIVE FUNCTION FCN
C GRADTL       --> TOLERANCE AT WHICH GRADIENT CONSIDERED CLOSE
C                  ENOUGH TO ZERO TO TERMINATE ALGORITHM
C ITNLIM      <--> MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
C NDIGIT      <--> NUMBER OF GOOD DIGITS IN OPTIMIZATION FUNCTION FCN
C EPSM         --> MACHINE EPSILON
C DLT         <--> TRUST REGION RADIUS
C METHOD      <--> ALGORITHM INDICATOR
C IEXP        <--> EXPENSE FLAG
C IAGFLG      <--> =1 IF ANALYTIC GRADIENT SUPPLIED
C IAHFLG      <--> =1 IF ANALYTIC HESSIAN SUPPLIED
C STEPMX      <--> MAXIMUM STEP SIZE
C MSG         <--> MESSAGE AND ERROR CODE
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C
      INTEGER N,ITNLIM,NDIGIT,METHOD,IEXP,IAGFLG,I
      INTEGER IAHFLG,MSG,IPR
      DOUBLE PRECISION FSCALE,GRADTL,EPSM,DLT,STEPMX,STPSIZ
      DOUBLE PRECISION X(N),TYPSIZ(N),SX(N)
C
C COMPUTE SCALE MATRIX
C
      DO 10 I=1,N
        IF(TYPSIZ(I).EQ.0.) TYPSIZ(I)=1.0
        IF(TYPSIZ(I).LT.0.) TYPSIZ(I)=-TYPSIZ(I)
        SX(I)=1.0/TYPSIZ(I)
   10 CONTINUE
C
C CHECK MAXIMUM STEP SIZE
C
      STPSIZ = 0.0
      DO 15 I = 1, N
         STPSIZ = STPSIZ + X(I)*X(I)*SX(I)*SX(I)
   15 CONTINUE
      STPSIZ =DSQRT(STPSIZ)
      STEPMX = MAX(1.0D3*STPSIZ, 1.0D3)
C
C CHECK NUMBER OF DIGITS OF ACCURACY IN FUNCTION FCN
      NDIGIT=-DLOG10(EPSM)
      RETURN
      END
      SUBROUTINE OPTDRV(NR,N,X,FCN,D1FCN,D2FCN,TYPSIZ,FSCALE,
     +     METHOD,IEXP,MSG,NDIGIT,ITNLIM,IAGFLG,IAHFLG,IPR,
     +     DLT,GRADTL,STEPMX,STEPTL,XPLS,FPLS,GPLS,ITRMCD,
     +     A,UDIAG,G,P,SX,WRK0,WRK1,WRK2,WRK3,
     +     AJA,ANLS,SHAT,VECT1,VECT2,VECT3,VECT4,VECT5,VECT6,
     +     NRM,NRN,MM,NN,IQ)
C
C PURPOSE
C -------
C DRIVER FOR NON-LINEAR OPTIMIZATION PROBLEM
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C X(N)         --> ON ENTRY: ESTIMATE TO A ROOT OF FCN
C FCN          --> NAME OF SUBROUTINE TO EVALUATE OPTIMIZATION FUNCTION
C                  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE
C                            FCN: R(N) --> R(1)
C D1FCN        --> (OPTIONAL) NAME OF SUBROUTINE TO EVALUATE GRADIENT
C                  OF FCN.  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE
C D2FCN        --> (OPTIONAL) NAME OF SUBROUTINE TO EVALUATE HESSIAN OF
C                  OF FCN.  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE
C TYPSIZ(N)    --> TYPICAL SIZE FOR EACH COMPONENT OF X
C FSCALE       --> ESTIMATE OF SCALE OF OBJECTIVE FUNCTION
C METHOD       --> ALGORITHM TO USE TO SOLVE MINIMIZATION PROBLEM
C                    =1 LINE SEARCH
C                    =2 DOUBLE DOGLEG
C                    =3 MORE-HEBDON
C IEXP         --> =1 IF OPTIMIZATION FUNCTION FCN IS EXPENSIVE TO
C                  EVALUATE, =0 OTHERWISE.  IF SET THEN HESSIAN WILL
C                  BE EVALUATED BY SECANT UPDATE INSTEAD OF
C                  ANALYTICALLY OR BY FINITE DIFFERENCES
C MSG         <--> ON INPUT:  (.GT.0) MESSAGE TO INHIBIT CERTAIN
C                    AUTOMATIC CHECKS
C                  ON OUTPUT: (.LT.0) ERROR CODE; =0 NO ERROR
C NDIGIT       --> NUMBER OF GOOD DIGITS IN OPTIMIZATION FUNCTION FCN
C ITNLIM       --> MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
C IAGFLG       --> =1 IF ANALYTIC GRADIENT SUPPLIED
C IAHFLG       --> =1 IF ANALYTIC HESSIAN SUPPLIED
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C DLT          --> TRUST REGION RADIUS
C GRADTL       --> TOLERANCE AT WHICH GRADIENT CONSIDERED CLOSE
C                  ENOUGH TO ZERO TO TERMINATE ALGORITHM
C STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE
C STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES
C                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM
C XPLS(N)     <--> ON EXIT:  XPLS IS LOCAL MINIMUM
C FPLS        <--> ON EXIT:  FUNCTION VALUE AT SOLUTION, XPLS
C GPLS(N)     <--> ON EXIT:  GRADIENT AT SOLUTION XPLS
C ITRMCD      <--  TERMINATION CODE
C A(N,N)       --> WORKSPACE FOR HESSIAN (OR ESTIMATE)
C                  AND ITS CHOLESKY DECOMPOSITION
C UDIAG(N)     --> WORKSPACE [FOR DIAGONAL OF HESSIAN]
C G(N)         --> WORKSPACE (FOR GRADIENT AT CURRENT ITERATE)
C P(N)         --> WORKSPACE FOR STEP
C SX(N)        --> WORKSPACE (FOR DIAGONAL SCALING MATRIX)
C WRK0(N)      --> WORKSPACE
C WRK1(N)      --> WORKSPACE
C WRK2(N)      --> WORKSPACE
C WRK3(N)      --> WORKSPACE
C
C
C INTERNAL VARIABLES
C ------------------
C ANALTL           TOLERANCE FOR COMPARISON OF ESTIMATED AND
C                  ANALYTICAL GRADIENTS AND HESSIANS
C EPSM             MACHINE EPSILON
C F                FUNCTION VALUE: FCN(X)
C ITNCNT           CURRENT ITERATION, K
C RNF              RELATIVE NOISE IN OPTIMIZATION FUNCTION FCN.
C                       NOISE=10.**(-NDIGIT)
C
      INTEGER NR,N,METHOD,IEXP,MSG,NDIGIT,ITNLIM,IAGFLG,IAHFLG,IPR
      INTEGER NRM,NRN,MM,NN,IQ,I,ITRMCD,ITNCNT,IRETCD,ICSCMX
      DOUBLE PRECISION FSCALE,DLT,GRADTL,STEPMX,STEPTL,F,FPLS
      DOUBLE PRECISION EPSM,DPMEPS,RNF,ANALTL,DLTSAV
      DOUBLE PRECISION AMUSAV,AMU,DLPSAV,DLTP,PHISAV,PHI,PHPSAV,PHIP0
      DOUBLE PRECISION X(N),XPLS(N),G(N),GPLS(N),P(N)
      DOUBLE PRECISION TYPSIZ(N),SX(N),WRK(1)
      DOUBLE PRECISION A(NR,N),UDIAG(N)
      DOUBLE PRECISION WRK0(N),WRK1(N),WRK2(N),WRK3(N)
      DOUBLE PRECISION AJA(NRM,NN),ANLS(NRM,N),SHAT(NRN,N)
      DOUBLE PRECISION VECT1(MM),VECT2(MM),VECT3(MM)
      DOUBLE PRECISION VECT4(MM),VECT5(MM),VECT6(MM)
      LOGICAL MXTAKE
      EXTERNAL FCN,D1FCN,D2FCN
C
C INITIALIZATION
C --------------
      DO 10 I=1,N
        P(I)=0.
   10 CONTINUE
      ITNCNT=0
      IRETCD=-1
      EPSM=DPMEPS()
      CALL OPTCHK(N,X,TYPSIZ,SX,FSCALE,GRADTL,ITNLIM,NDIGIT,EPSM,
     +     DLT,METHOD,IEXP,IAGFLG,IAHFLG,STEPMX,MSG,IPR)
      IF(MSG.LT.0) RETURN
      RNF=MAX(10.0D0**(-NDIGIT),EPSM)
      ANALTL=MAX(1.D-2,SQRT(RNF))
C
C EVALUATE FCN(X)
C
      CALL FCN(N,X,F,AJA,ANLS,SHAT,VECT1,VECT2,VECT3,VECT4,
     +  VECT5,VECT6,NRM,NRN,MM,NN,IQ)
C
C EVALUATE ANALYTIC OR FINITE DIFFERENCE GRADIENT AND CHECK ANALYTIC
C GRADIENT, IF REQUESTED.
C
      IF (IAGFLG .EQ. 1) THEN
         CALL D1FCN (N, X, G,AJA,ANLS,SHAT,VECT1,VECT2,VECT3,
     +               VECT4,VECT5,VECT6,NRM,NRN,MM,NN,IQ)
      ELSE
        CALL FSTOFD (1, 1, N, X, FCN, F, G, SX, RNF, WRK, 1,
     +  AJA,ANLS,SHAT,VECT1,VECT2,VECT3,VECT4,VECT5,VECT6,
     +  NRM,NRN,MM,NN,IQ)
      ENDIF
C
      CALL OPTSTP(N,X,F,G,WRK1,ITNCNT,ICSCMX,
     +            ITRMCD,GRADTL,STEPTL,SX,FSCALE,ITNLIM,IRETCD,MXTAKE,
     +            IPR,MSG)
      IF(ITRMCD.NE.0) GO TO 700
C
      IF(IEXP.NE.1) GO TO 80
C
C IF OPTIMIZATION FUNCTION EXPENSIVE TO EVALUATE (IEXP=1), THEN
C HESSIAN WILL BE OBTAINED BY SECANT UPDATES.  GET INITIAL HESSIAN.
C
      CALL HSNINT(NR,N,A,SX,METHOD)
      GO TO 90
   80 CONTINUE
C
C EVALUATE ANALYTIC OR FINITE DIFFERENCE HESSIAN AND CHECK ANALYTIC
C HESSIAN IF REQUESTED (ONLY IF USER-SUPPLIED ANALYTIC HESSIAN
C ROUTINE D2FCN FILLS ONLY LOWER TRIANGULAR PART AND DIAGONAL OF A).
C
      IF (IAHFLG .EQ. 0) THEN
         IF (IAGFLG .EQ. 1) CALL FSTOFD (NR, N, N, X, D1FCN, G, A, SX,
     1      RNF, WRK1, 3,AJA,ANLS,SHAT,VECT1,VECT2,VECT3,VECT4,VECT5,
     1      VECT6,NRM,NRN,MM,NN,IQ)
         IF (IAGFLG .NE. 1) CALL SNDOFD (NR, N, X, FCN, F, A, SX, RNF,
     1      WRK1, WRK2)
C
      ENDIF
C
   90 CONTINUE
C
C
C ITERATION
C ---------
  100 ITNCNT=ITNCNT+1
C
C FIND PERTURBED LOCAL MODEL HESSIAN AND ITS LL+ DECOMPOSITION
C (SKIP THIS STEP IF LINE SEARCH OR DOGSTEP TECHNIQUES BEING USED WITH
C SECANT UPDATES.  CHOLESKY DECOMPOSITION L ALREADY OBTAINED FROM
C SECFAC.)
C
      IF(IEXP.EQ.1 .AND. METHOD.NE.3) GO TO 105
  103   CALL CHLHSN(NR,N,A,EPSM,SX,UDIAG)
  105 CONTINUE
C
C SOLVE FOR NEWTON STEP:  AP=-G
C
      DO 110 I=1,N
        WRK1(I)=-G(I)
  110 CONTINUE
      CALL LLTSLV(NR,N,A,P,WRK1)
C
C DECIDE WHETHER TO ACCEPT NEWTON STEP  XPLS=X + P
C OR TO CHOOSE XPLS BY A GLOBAL STRATEGY.
C
      IF (IAGFLG .NE. 0 .OR. METHOD .EQ. 1) GO TO 111
      DLTSAV = DLT
      IF (METHOD .EQ. 2) GO TO 111
      AMUSAV = AMU
      DLPSAV = DLTP
      PHISAV = PHI
      PHPSAV = PHIP0
  111 CONTINUE
      IF(METHOD.EQ.2)
     +     CALL DOGDRV(NR,N,X,F,G,A,P,XPLS,FPLS,FCN,SX,STEPMX,
     +     STEPTL,DLT,IRETCD,MXTAKE,WRK0,WRK1,WRK2,WRK3,IPR,
     +     AJA,ANLS,SHAT,VECT1,VECT2,VECT3,VECT4,VECT5,VECT6,
     +     NRM,NRN,MM,NN,IQ)
      IF(METHOD.EQ.3)
     +     CALL HOOKDR(NR,N,X,F,G,A,UDIAG,P,XPLS,FPLS,FCN,SX,STEPMX,
     +     STEPTL,DLT,IRETCD,MXTAKE,AMU,DLTP,PHI,PHIP0,WRK0,
     +     WRK1,WRK2,EPSM,ITNCNT,IPR)
C
C IF COULD NOT FIND SATISFACTORY STEP AND FORWARD DIFFERENCE
C GRADIENT WAS USED, RETRY USING CENTRAL DIFFERENCE GRADIENT.
C
      IF (IRETCD .NE. 1 .OR. IAGFLG .NE. 0) GO TO 112
C     IF (IRETCD .EQ. 1 .AND. IAGFLG .EQ. 0)
C     THEN
C
C        SET IAGFLG FOR CENTRAL DIFFERENCES
C
         IAGFLG = -1
         CALL FSTOCD (N, X, FCN, SX, RNF, G)
         IF (METHOD .EQ. 1) GO TO 105
         DLT = DLTSAV
         IF (METHOD .EQ. 2) GO TO 105
         AMU = AMUSAV
         DLTP = DLPSAV
         PHI = PHISAV
         PHIP0 = PHPSAV
         GO TO 103
C
C CALCULATE STEP FOR OUTPUT
C
  112 CONTINUE
      DO 114 I = 1, N
         P(I) = XPLS(I) - X(I)
  114 CONTINUE
C
C CALCULATE GRADIENT AT XPLS
C
      IF (IAGFLG .EQ. (-1)) GO TO 116
      IF (IAGFLG .EQ. 0) GO TO 118
C
C ANALYTIC GRADIENT
C
      CALL D1FCN (N, XPLS, GPLS,
     +  AJA,ANLS,SHAT,VECT1,VECT2,VECT3,VECT4,VECT5,VECT6,
     +  NRM,NRN,MM,NN,IQ)
      GO TO 120
C
C CENTRAL DIFFERENCE GRADIENT
C
  116 CALL FSTOCD (N, XPLS, FCN, SX, RNF, GPLS)
      GO TO 120
C
C FORWARD DIFFERENCE GRADIENT
C
  118 CALL FSTOFD (1, 1, N, XPLS, FCN, FPLS, GPLS, SX, RNF, WRK, 1,
     +  AJA,ANLS,SHAT,VECT1,VECT2,VECT3,VECT4,VECT5,VECT6,
     +  NRM,NRN,MM,NN,IQ)
  120 CONTINUE
C
C CHECK WHETHER STOPPING CRITERIA SATISFIED
C
      CALL OPTSTP(N,XPLS,FPLS,GPLS,X,ITNCNT,ICSCMX,
     +            ITRMCD,GRADTL,STEPTL,SX,FSCALE,ITNLIM,IRETCD,MXTAKE,
     +            IPR,MSG)
      IF(ITRMCD.NE.0) GO TO 690
C
C EVALUATE HESSIAN AT XPLS
C
      IF(IEXP.NE.0) GO TO 150
      IF(IAHFLG.EQ.1) GO TO 140
      IF(IAGFLG.EQ.1)
     +     CALL FSTOFD(NR,N,N,XPLS,D1FCN,GPLS,A,SX,RNF,WRK1,3,
     +     AJA,ANLS,SHAT,VECT1,VECT2,VECT3,VECT4,VECT5,VECT6,
     +     NRM,NRN,MM,NN,IQ)
      IF(IAGFLG.NE.1) CALL SNDOFD(NR,N,XPLS,FCN,FPLS,A,SX,
     +                            RNF,WRK1,WRK2)
      GO TO 150
  140 CALL D2FCN(NR,N,XPLS,A)
  150 CONTINUE
C
C X <-- XPLS  AND  G <-- GPLS  AND  F <-- FPLS
C
      F=FPLS
      DO 160 I=1,N
        X(I)=XPLS(I)
        G(I)=GPLS(I)
  160 CONTINUE
      GO TO 100
C
C TERMINATION
C -----------
C RESET XPLS,FPLS,GPLS,  IF PREVIOUS ITERATE SOLUTION
C
  690 IF(ITRMCD.NE.3) GO TO 710
  700 CONTINUE
      FPLS=F
      DO 705 I=1,N
        XPLS(I)=X(I)
        GPLS(I)=G(I)
  705 CONTINUE
  710 CONTINUE
      MSG=0
      RETURN
      END
      SUBROUTINE OPTIF9(NR,N,X,FCN,D1FCN,D2FCN,TYPSIZ,FSCALE,METHOD,
     +     IEXP,MSG,NDIGIT,ITNLIM,IAGFLG,IAHFLG,IPR,DLT,GRADTL,
     +     STEPMX,STEPTL,XPLS,FPLS,GPLS,ITRMCD,A,WRK,AJA,ANLS,SHAT,
     +     VECT1,VECT2,VECT3,VECT4,VECT5,VECT6,NRM,NRN,MM,NN,IQ)
C
C PURPOSE
C -------
C PROVIDE COMPLETE INTERFACE TO MINIMIZATION PACKAGE.
C USER HAS FULL CONTROL OVER OPTIONS.
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C X(N)         --> ON ENTRY: ESTIMATE TO A ROOT OF FCN
C FCN          --> NAME OF SUBROUTINE TO EVALUATE OPTIMIZATION FUNCTION
C                  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE
C                            FCN: R(N) --> R(1)
C D1FCN        --> (OPTIONAL) NAME OF SUBROUTINE TO EVALUATE GRADIENT
C                  OF FCN.  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE
C D2FCN        --> (OPTIONAL) NAME OF SUBROUTINE TO EVALUATE HESSIAN OF
C                  OF FCN.  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE
C TYPSIZ(N)    --> TYPICAL SIZE FOR EACH COMPONENT OF X
C FSCALE       --> ESTIMATE OF SCALE OF OBJECTIVE FUNCTION
C METHOD       --> ALGORITHM TO USE TO SOLVE MINIMIZATION PROBLEM
C                    =1 LINE SEARCH
C                    =2 DOUBLE DOGLEG
C                    =3 MORE-HEBDON
C IEXP         --> =1 IF OPTIMIZATION FUNCTION FCN IS EXPENSIVE TO
C                  EVALUATE, =0 OTHERWISE.  IF SET THEN HESSIAN WILL
C                  BE EVALUATED BY SECANT UPDATE INSTEAD OF
C                  ANALYTICALLY OR BY FINITE DIFFERENCES
C MSG         <--> ON INPUT:  (.GT.0) MESSAGE TO INHIBIT CERTAIN
C                    AUTOMATIC CHECKS
C                  ON OUTPUT: (.LT.0) ERROR CODE; =0 NO ERROR
C NDIGIT       --> NUMBER OF GOOD DIGITS IN OPTIMIZATION FUNCTION FCN
C ITNLIM       --> MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
C IAGFLG       --> =1 IF ANALYTIC GRADIENT SUPPLIED
C IAHFLG       --> =1 IF ANALYTIC HESSIAN SUPPLIED
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C DLT          --> TRUST REGION RADIUS
C GRADTL       --> TOLERANCE AT WHICH GRADIENT CONSIDERED CLOSE
C                  ENOUGH TO ZERO TO TERMINATE ALGORITHM
C STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE
C STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES
C                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM
C XPLS(N)     <--> ON EXIT:  XPLS IS LOCAL MINIMUM
C FPLS        <--> ON EXIT:  FUNCTION VALUE AT SOLUTION, XPLS
C GPLS(N)     <--> ON EXIT:  GRADIENT AT SOLUTION XPLS
C ITRMCD      <--  TERMINATION CODE
C A(N,N)       --> WORKSPACE FOR HESSIAN (OR ESTIMATE)
C                  AND ITS CHOLESKY DECOMPOSITION
C WRK(N,8)     --> WORKSPACE
C
      INTEGER NR,N,METHOD,IEXP,MSG,NDIGIT,ITNLIM,IAGFLG,IAHFLG,IPR
      INTEGER NRM,NRN,MM,NN,IQ,ITRMCD
      DOUBLE PRECISION FSCALE,DLT,GRADTL,STEPMX,STEPTL,FPLS
      DOUBLE PRECISION X(N),XPLS(N),GPLS(N),TYPSIZ(N)
      DOUBLE PRECISION A(NR,N),WRK(NR,8)
      DOUBLE PRECISION AJA(NRM,NN),ANLS(NRM,N),SHAT(NRN,N)
      DOUBLE PRECISION VECT1(MM),VECT2(MM),VECT3(MM)
      DOUBLE PRECISION VECT4(MM),VECT5(MM),VECT6(MM)
      EXTERNAL FCN,D1FCN,D2FCN
C
C EQUIVALENCE WRK(N,1) = UDIAG(N)
C             WRK(N,2) = G(N)
C             WRK(N,3) = P(N)
C             WRK(N,4) = SX(N)
C             WRK(N,5) = WRK0(N)
C             WRK(N,6) = WRK1(N)
C             WRK(N,7) = WRK2(N)
C             WRK(N,8) = WRK3(N)
C
      CALL OPTDRV(NR,N,X,FCN,D1FCN,D2FCN,TYPSIZ,FSCALE,
     +     METHOD,IEXP,MSG,NDIGIT,ITNLIM,IAGFLG,IAHFLG,IPR,
     +     DLT,GRADTL,STEPMX,STEPTL,XPLS,FPLS,GPLS,ITRMCD,
     +     A,WRK(1,1),WRK(1,2),WRK(1,3),WRK(1,4),WRK(1,5),
     +     WRK(1,6),WRK(1,7),WRK(1,8),
     +     AJA,ANLS,SHAT,VECT1,VECT2,VECT3,VECT4,VECT5,VECT6,
     +     NRM,NRN,MM,NN,IQ)
      RETURN
      END
      SUBROUTINE OPTSTP(N,XPLS,FPLS,GPLS,X,ITNCNT,ICSCMX,
     +      ITRMCD,GRADTL,STEPTL,SX,FSCALE,ITNLIM,IRETCD,
     +      MXTAKE,IPR,MSG)
C
C UNCONSTRAINED MINIMIZATION STOPPING CRITERIA
C --------------------------------------------
C FIND WHETHER THE ALGORITHM SHOULD TERMINATE, DUE TO ANY
C OF THE FOLLOWING:
C 1) PROBLEM SOLVED WITHIN USER TOLERANCE
C 2) CONVERGENCE WITHIN USER TOLERANCE
C 3) ITERATION LIMIT REACHED
C 4) DIVERGENCE OR TOO RESTRICTIVE MAXIMUM STEP (STEPMX) SUSPECTED
C
C
C PARAMETERS
C ----------
C N            --> DIMENSION OF PROBLEM
C XPLS(N)      --> NEW ITERATE X[K]
C FPLS         --> FUNCTION VALUE AT NEW ITERATE F(XPLS)
C GPLS(N)      --> GRADIENT AT NEW ITERATE, G(XPLS), OR APPROXIMATE
C X(N)         --> OLD ITERATE X[K-1]
C ITNCNT       --> CURRENT ITERATION K
C ICSCMX      <--> NUMBER CONSECUTIVE STEPS .GE. STEPMX
C                  [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C ITRMCD      <--  TERMINATION CODE
C GRADTL       --> TOLERANCE AT WHICH RELATIVE GRADIENT CONSIDERED CLOSE
C                  ENOUGH TO ZERO TO TERMINATE ALGORITHM
C STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES
C                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C FSCALE       --> ESTIMATE OF SCALE OF OBJECTIVE FUNCTION
C ITNLIM       --> MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
C IRETCD       --> RETURN CODE
C MXTAKE       --> BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C MSG          --> IF MSG INCLUDES A TERM 8, SUPPRESS OUTPUT
C
C
      INTEGER N,MSG,ITNLIM,IPR,I
      INTEGER JTRMCD,ITNCNT,IRETCD,ITRMCD,ICSCMX
      DOUBLE PRECISION FSCALE,GRADTL,STEPTL,FPLS,D,RGX
      DOUBLE PRECISION RELGRD,RELSTP,RSX
      DOUBLE PRECISION SX(N)
      DOUBLE PRECISION XPLS(N),GPLS(N),X(N)
      LOGICAL MXTAKE
C
      ITRMCD=0
C
C LAST GLOBAL STEP FAILED TO LOCATE A POINT LOWER THAN X
      IF(IRETCD.NE.1) GO TO 50
C     IF(IRETCD.EQ.1)
C     THEN
        JTRMCD=3
        GO TO 600
C     ENDIF
   50 CONTINUE
C
C FIND DIRECTION IN WHICH RELATIVE GRADIENT MAXIMUM.
C CHECK WHETHER WITHIN TOLERANCE
C
      D=MAX(ABS(FPLS),FSCALE)
      RGX=0.0
      DO 100 I=1,N
        RELGRD=ABS(GPLS(I))*MAX(ABS(XPLS(I)),1./SX(I))/D
        RGX=MAX(RGX,RELGRD)
  100 CONTINUE
      JTRMCD=1
      IF(RGX.LE.GRADTL) GO TO 600
C
      IF(ITNCNT.EQ.0) RETURN
C
C FIND DIRECTION IN WHICH RELATIVE STEPSIZE MAXIMUM
C CHECK WHETHER WITHIN TOLERANCE.
C
      RSX=0.0
      DO 120 I=1,N
        RELSTP=ABS(XPLS(I)-X(I))/MAX(ABS(XPLS(I)),1./SX(I))
        RSX=MAX(RSX,RELSTP)
  120 CONTINUE
      JTRMCD=2
      IF(RSX.LE.STEPTL) GO TO 600
C
C CHECK ITERATION LIMIT
C
      JTRMCD=4
      IF(ITNCNT.GE.ITNLIM) GO TO 600
C
C CHECK NUMBER OF CONSECUTIVE STEPS \ STEPMX
C
      IF(MXTAKE) GO TO 140
C     IF(.NOT.MXTAKE)
C     THEN
        ICSCMX=0
        RETURN
C     ELSE
  140   CONTINUE
        ICSCMX=ICSCMX+1
        IF(ICSCMX.LT.5) RETURN
        JTRMCD=5
C     ENDIF
C
  600 ITRMCD=JTRMCD
C
      RETURN
      END
      SUBROUTINE SNDOFD(NR,N,XPLS,FCN,FPLS,A,SX,RNOISE,STEPSZ,ANBR)
C PURPOSE
C -------
C FIND SECOND ORDER FORWARD FINITE DIFFERENCE APPROXIMATION "A"
C TO THE SECOND DERIVATIVE (HESSIAN) OF THE FUNCTION DEFINED BY THE SUBP
C "FCN" EVALUATED AT THE NEW ITERATE "XPLS"
C
C FOR OPTIMIZATION USE THIS ROUTINE TO ESTIMATE
C 1) THE SECOND DERIVATIVE (HESSIAN) OF THE OPTIMIZATION FUNCTION
C    IF NO ANALYTICAL USER FUNCTION HAS BEEN SUPPLIED FOR EITHER
C    THE GRADIENT OR THE HESSIAN AND IF THE OPTIMIZATION FUNCTION
C    "FCN" IS INEXPENSIVE TO EVALUATE.
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C XPLS(N)      --> NEW ITERATE:   X[K]
C FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION
C FPLS         --> FUNCTION VALUE AT NEW ITERATE, F(XPLS)
C A(N,N)      <--  FINITE DIFFERENCE APPROXIMATION TO HESSIAN
C                  ONLY LOWER TRIANGULAR MATRIX AND DIAGONAL
C                  ARE RETURNED
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C RNOISE       --> RELATIVE NOISE IN FNAME [F(X)]
C STEPSZ(N)    --> WORKSPACE (STEPSIZE IN I-TH COMPONENT DIRECTION)
C ANBR(N)      --> WORKSPACE (NEIGHBOR IN I-TH DIRECTION)
C
C
      INTEGER NR,N,I,J,IP1
      DOUBLE PRECISION FPLS,RNOISE,OV3,XTMPI,XTMPJ,FHAT
      DOUBLE PRECISION XPLS(N)
      DOUBLE PRECISION SX(N)
      DOUBLE PRECISION STEPSZ(N),ANBR(N)
      DOUBLE PRECISION A(NR,N)
      EXTERNAL FCN
C
C FIND I-TH STEPSIZE AND EVALUATE NEIGHBOR IN DIRECTION
C OF I-TH UNIT VECTOR.
C
      OV3 = 1.0/3.0
      DO 10 I=1,N
        STEPSZ(I)=RNOISE**OV3 * MAX(ABS(XPLS(I)),1./SX(I))
        XTMPI=XPLS(I)
        XPLS(I)=XTMPI+STEPSZ(I)
        CALL FCN(N,XPLS,ANBR(I))
        XPLS(I)=XTMPI
   10 CONTINUE
C
C CALCULATE COLUMN I OF A
C
      DO 30 I=1,N
        XTMPI=XPLS(I)
        XPLS(I)=XTMPI+2.0*STEPSZ(I)
        CALL FCN(N,XPLS,FHAT)
        A(I,I)=((FPLS-ANBR(I))+(FHAT-ANBR(I)))/(STEPSZ(I)*STEPSZ(I))
C
C CALCULATE SUB-DIAGONAL ELEMENTS OF COLUMN
        IF(I.EQ.N) GO TO 25
        XPLS(I)=XTMPI+STEPSZ(I)
        IP1=I+1
        DO 20 J=IP1,N
          XTMPJ=XPLS(J)
          XPLS(J)=XTMPJ+STEPSZ(J)
          CALL FCN(N,XPLS,FHAT)
          A(J,I)=((FPLS-ANBR(I))+(FHAT-ANBR(J)))/(STEPSZ(I)*STEPSZ(J))
          XPLS(J)=XTMPJ
   20   CONTINUE
   25   XPLS(I)=XTMPI
   30 CONTINUE
      RETURN
      END
      SUBROUTINE TREGUP(NR,N,X,F,G,A,FCN,SC,SX,NWTAKE,STEPMX,STEPTL,
     + DLT,IRETCD,XPLSP,FPLSP,XPLS,FPLS,MXTAKE,IPR,METHOD,UDIAG,
     + AJA,ANLS,SHAT,VECT1,VECT2,VECT3,VECT4,VECT5,VECT6,NRM,NRN,MM,
     + NN,IQ)
C
C PURPOSE
C -------
C DECIDE WHETHER TO ACCEPT XPLS=X+SC AS THE NEXT ITERATE AND UPDATE THE
C TRUST REGION DLT.
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C X(N)         --> OLD ITERATE X[K-1]
C F            --> FUNCTION VALUE AT OLD ITERATE, F(X)
C G(N)         --> GRADIENT AT OLD ITERATE, G(X), OR APPROXIMATE
C A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN IN
C                  LOWER TRIANGULAR PART AND DIAGONAL.
C                  HESSIAN OR APPROX IN UPPER TRIANGULAR PART
C FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION
C SC(N)        --> CURRENT STEP
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C NWTAKE       --> BOOLEAN, =.TRUE. IF NEWTON STEP TAKEN
C STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE
C STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES
C                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM
C DLT         <--> TRUST REGION RADIUS
C IRETCD      <--> RETURN CODE
C                    =0 XPLS ACCEPTED AS NEXT ITERATE;
C                       DLT TRUST REGION FOR NEXT ITERATION.
C                    =1 XPLS UNSATISFACTORY BUT ACCEPTED AS NEXT ITERATE
C                       BECAUSE XPLS-X .LT. SMALLEST ALLOWABLE
C                       STEP LENGTH.
C                    =2 F(XPLS) TOO LARGE.  CONTINUE CURRENT ITERATION
C                       WITH NEW REDUCED DLT.
C                    =3 F(XPLS) SUFFICIENTLY SMALL, BUT QUADRATIC MODEL
C                       PREDICTS F(XPLS) SUFFICIENTLY WELL TO CONTINUE
C                       CURRENT ITERATION WITH NEW DOUBLED DLT.
C XPLSP(N)    <--> WORKSPACE [VALUE NEEDS TO BE RETAINED BETWEEN
C                  SUCCESSIVE CALLS OF K-TH GLOBAL STEP]
C FPLSP       <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C XPLS(N)     <--  NEW ITERATE X[K]
C FPLS        <--  FUNCTION VALUE AT NEW ITERATE, F(XPLS)
C MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C METHOD       --> ALGORITHM TO USE TO SOLVE MINIMIZATION PROBLEM
C                    =1 LINE SEARCH
C                    =2 DOUBLE DOGLEG
C                    =3 MORE-HEBDON
C UDIAG(N)     --> DIAGONAL OF HESSIAN IN A(.,.)
C
      INTEGER NR,N,IRETCD,IPR,METHOD,NRM,NRN,MM,NN,IQ,I,J,IP1
      DOUBLE PRECISION STEPMX,STEPTL,DLT,FPLSP,FPLS,SLP,RLN
      DOUBLE PRECISION DLTMP,DLTFP,TEMP,DLTF,F
      DOUBLE PRECISION X(N),XPLS(N),G(N)
      DOUBLE PRECISION SX(N),SC(N),XPLSP(N)
      DOUBLE PRECISION A(NR,N)
      DOUBLE PRECISION  UDIAG(N)
      DOUBLE PRECISION AJA(NRM,NN),ANLS(NRM,N),SHAT(NRN,N),VECT1(MM)
      DOUBLE PRECISION VECT2(MM),VECT3(MM),VECT4(MM),VECT5(MM)
      DOUBLE PRECISION VECT6(MM)
      DOUBLE PRECISION DDOT
      LOGICAL NWTAKE,MXTAKE
      EXTERNAL FCN
C
      MXTAKE=.FALSE.
      DO 100 I=1,N
        XPLS(I)=X(I)+SC(I)
  100 CONTINUE
      CALL FCN(N,XPLS,FPLS,AJA,ANLS,SHAT,VECT1,VECT2,VECT3,
     +   VECT4,VECT5,VECT6,NRM,NRN,MM,NN,IQ)
      DLTF=FPLS-F
      SLP=DDOT(N,G,1,SC,1)
C
C NEXT STATEMENT ADDED FOR CASE OF COMPILERS WHICH DO NOT OPTIMIZE
C EVALUATION OF NEXT "IF" STATEMENT (IN WHICH CASE FPLSP COULD BE
C UNDEFINED).
      IF(IRETCD.EQ.4) FPLSP=0.0
      IF(IRETCD.NE.3 .OR. (FPLS.LT.FPLSP .AND. DLTF.LE. 1.E-4*SLP))
     +                                                     GO TO 130
C     IF(IRETCD.EQ.3 .AND. (FPLS.GE.FPLSP .OR. DLTF.GT. 1.E-4*SLP))
C     THEN
C
C       RESET XPLS TO XPLSP AND TERMINATE GLOBAL STEP
C
        IRETCD=0
        DO 110 I=1,N
          XPLS(I)=XPLSP(I)
  110   CONTINUE
        FPLS=FPLSP
        DLT=.5*DLT
        GO TO 230
C     ELSE
C
C       FPLS TOO LARGE
C
  130   IF(DLTF.LE. 1.E-4*SLP) GO TO 170
C       IF(DLTF.GT. 1.E-4*SLP)
C       THEN
          RLN=0.
          DO 140 I=1,N
            RLN=MAX(RLN,ABS(SC(I))/MAX(ABS(XPLS(I)),1./SX(I)))
  140     CONTINUE
          IF(RLN.GE.STEPTL) GO TO 150
C         IF(RLN.LT.STEPTL)
C         THEN
C
C           CANNOT FIND SATISFACTORY XPLS SUFFICIENTLY DISTINCT FROM X
C
            IRETCD=1
            GO TO 230
C         ELSE
C
C           REDUCE TRUST REGION AND CONTINUE GLOBAL STEP
C
  150       IRETCD=2
            DLTMP=-SLP*DLT/(2.*(DLTF-SLP))
            IF(DLTMP.GE. .1*DLT) GO TO 155
C           IF(DLTMP.LT. .1*DLT)
C           THEN
              DLT=.1*DLT
              GO TO 160
C           ELSE
  155         DLT=DLTMP
C           ENDIF
  160       CONTINUE
            GO TO 230
C         ENDIF
C       ELSE
C
C         FPLS SUFFICIENTLY SMALL
C
  170     CONTINUE
          DLTFP=0.
          IF (METHOD .EQ. 3) GO TO 180
C
C         IF (METHOD .EQ. 2)
C         THEN
C
          DO 177 I = 1, N
             TEMP = 0.0
             DO 173 J = I, N
                TEMP = TEMP + (A(J, I)*SC(J))
  173        CONTINUE
             DLTFP = DLTFP + TEMP*TEMP
  177     CONTINUE
          GO TO 190
C
C         ELSE
C
  180     DO 187 I = 1, N
             DLTFP = DLTFP + UDIAG(I)*SC(I)*SC(I)
             IF (I .EQ. N) GO TO 187
             TEMP = 0
             IP1 = I + 1
             DO 183 J = IP1, N
                TEMP = TEMP + A(I, J)*SC(I)*SC(J)
  183        CONTINUE
             DLTFP = DLTFP + 2.0*TEMP
  187     CONTINUE
C
C         END IF
C
  190     DLTFP = SLP + DLTFP/2.0
          IF(IRETCD.EQ.2 .OR. (ABS(DLTFP-DLTF).GT. .1*ABS(DLTF))
     +         .OR. NWTAKE .OR. (DLT.GT. .99*STEPMX)) GO TO 210
C         IF(IRETCD.NE.2 .AND. (ABS(DLTFP-DLTF) .LE. .1*ABS(DLTF))
C    +         .AND. (.NOT.NWTAKE) .AND. (DLT.LE. .99*STEPMX))
C         THEN
C
C           DOUBLE TRUST REGION AND CONTINUE GLOBAL STEP
C
            IRETCD=3
            DO 200 I=1,N
              XPLSP(I)=XPLS(I)
  200       CONTINUE
            FPLSP=FPLS
            DLT=MIN(2.0D0*DLT,STEPMX)
            GO TO 230
C         ELSE
C
C           ACCEPT XPLS AS NEXT ITERATE.  CHOOSE NEW TRUST REGION.
C
  210       CONTINUE
            IRETCD=0
            IF(DLT.GT. .99*STEPMX) MXTAKE=.TRUE.
            IF(DLTF.LT. .1*DLTFP) GO TO 220
C           IF(DLTF.GE. .1*DLTFP)
C           THEN
C
C             DECREASE TRUST REGION FOR NEXT ITERATION
C
              DLT=.5*DLT
              GO TO 230
C           ELSE
C
C             CHECK WHETHER TO INCREASE TRUST REGION FOR NEXT ITERATION
C
  220         IF(DLTF.LE. .75*DLTFP) DLT=MIN(2.*DLT,STEPMX)
C           ENDIF
C         ENDIF
C       ENDIF
C     ENDIF
  230 CONTINUE
      RETURN
      END


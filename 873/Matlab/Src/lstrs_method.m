%
%   File:              lstrs_method.m
%
%   Description:       LSTRS: A Matrix-Free Method for the
%                      Large-Scale Trust-Region Subproblem
%               
%                      min   1/2 {x'Hx + g'x}
%                      s.t.  ||x|| <= Delta
%
%                      Based on: M. Rojas, S.A. Santos and D.C.Sorensen.
%                                A New Matrix-Free Algorithm for the Large-Scale
%                                Trust-Region Subproblem, SIAM J. Optim., 
%                                11(3):611-646, 2000.
%
%   Authors:           Marielba Rojas
%                        mr@imm.dtu.dk
%                      Sandra Santos
%                      Danny Sorensen
%
%   Version:           1.2
%
%   System:            MATLAB 6.0 or higher 
%
%   Date:              15 March 2007
%
%   Functions called:  adjust_alpha, b_epairs, convergence, func2str, 
%                      init_lo_bounds, init_up_bounds, inter_point, isa,
%                      isempty, length, min, norm, output, strcmp, 
%                      upd_alpha_safe, upd_deltaU, upd_param0, upd_paramk
%
%   Call: [x,lambda,info,moreinfo] = ...
%         lstrs_method(H,g,Delta,epsilon,eigensolver,lopts,Hpar,eigensolverpar);
%

function [x,lambda,info,moreinfo] = ...
         lstrs_method(H,g,Delta,epsilon,eigensolver,lopts,Hpar,eigensolverpar)

%
% Declare and initialize global variable to count matrix-vector products
% Routines that use it: lstrs_method, matvec, output
%
global mvp_lstrs;
mvp_lstrs = 0;

%
%  Extract some values from lopts
%
     message_level = lopts.message_level;
     maxiter       = lopts.maxiter;

%
%  Useful quantities and fixed part of bordered matrix
%
     norg      = norm(g);

     Deltasqp1 = Delta^2 + 1;

     Balpha.H    = H;
     Balpha.Hpar = Hpar;
     Balpha.dim  = length(g) + 1;
     Balpha.g    = g;
     Balpha.bord = 1;

     clear H Hpar g


%
%  Initial upper bounds for alpha and delta1
%
%      and
%
%  Initial guess for alpha
%

     [alphaU,deltaU] = init_up_bounds(Balpha,norg,Delta,lopts.deltaU,message_level);
     

     if (isa(lopts.alpha,'char'))
        if (strcmp(lopts.alpha,'min'))
           alpha  = min([0,alphaU]);
        else
           alpha  = deltaU;
        end
     else
        alpha = lopts.alpha;
     end

     Balpha.alpha = alpha;

%
%  Chooses the initial eigensolver
%
     if (~isa(eigensolver,'char'))
        sesolver = func2str(eigensolver);
     else
        sesolver = eigensolver;
     end

     if (strcmp(sesolver,'tcheigs_lstrs_gateway'))
        initial_eigensolver = @eigs_lstrs_gateway;
     else
        initial_eigensolver = eigensolver;
     end

%
%  Computes iterate 0 with initial guess for alpha
%
     [nconv,epair1,epair2,v1] = ...
      b_epairs(Balpha,initial_eigensolver,eigensolverpar,lopts);
     eigensolverpar.v0 = v1;

     if (nconv == 0)
        disp('Initialization of LSTRS: No eigenvalue converged!'); 
        disp('You might want to try the following:')
        disp('- Ask for a lower accuracy in the eigenpairs');
        disp('- Increase the number of Lanczos vectors');
	error(' ');
     end

     adj = 1;

%
%  Initial lower bound for alpha
%
     alphaL = init_lo_bounds(epair1,norg,Delta,message_level);

%
%  Initial value for z, approximate vector in S1
%
     z    = 0;

     iter = 0;

%
%  Updates upper bound for delta1
%
     deltaU = upd_deltaU(epair1,Balpha.g,deltaU,message_level);

%
%  Adjusts the parameter alpha0
%

     [alpha,alphaU,nconv,epair1,epair2,z,deltaU,adjk,iterate,v1] =       ...
      adjust_alpha(nconv,epair1,epair2,Balpha,norg,alphaL,alphaU,deltaU, ...
                   epsilon,z,eigensolver,eigensolverpar,lopts);

     Balpha.alpha = alpha;

     adj = adj + adjk;

     [convergence_cond,xtilde,lambdatilde,stopcond] =             ...
      convergence(Balpha,nconv,epair1,epair2,Delta,Deltasqp1,     ...
                  alphaL,alphaU,iter,maxiter,iterate,epsilon,message_level);

     if (~ convergence_cond)
%
%  Chooses interpolation point
%
        ipointb = inter_point(epair1,epair2,iterate,Balpha.g,message_level);

%
%  Updates bounds for alpha (safeguarding interval)
%
        [alphaL,alphaU] = ...
         upd_alpha_safe(alpha,alphaL,alphaU,iterate,ipointb,Delta,message_level);

        alpha = upd_param0(ipointb,Delta,alpha,alphaL,alphaU,deltaU,message_level);

        Balpha.alpha = alpha;

%
%  (Possibly) adjusts accuracy of eigenpairs
%

       if (isfield(eigensolverpar,'tol') & ~isempty(lopts.maxeigentol))
          eigensolverpar.tol = ...
          adjust_eigentol(eigensolverpar.tol,iterate,epair1,epair2,Delta, ...
                          lopts.maxeigentol,message_level);
       end

%
%  Computes iterate 1 with alpha from one-point interpolation
%

       [nconv,epair1,epair2,v1] = ...
        b_epairs(Balpha,eigensolver,eigensolverpar,lopts);
       eigensolverpar.v0 = v1;

       adj  = adj + 1;

       if (nconv > 0)
          deltaU = upd_deltaU(epair1,Balpha.g,deltaU,message_level);
       end

       [alpha,alphaU,nconv,epair1,epair2,z,deltaU,adjk,iterate,v1] =        ...
        adjust_alpha(nconv,epair1,epair2,Balpha,norg,alphaL,alphaU,deltaU,  ...
                     epsilon,z,eigensolver,eigensolverpar,lopts);

       Balpha.alpha = alpha;

       adj  = adj + adjk;
       iter = 1;

       [convergence_cond,xtilde,lambdatilde,stopcond] =             ...
        convergence(Balpha,nconv,epair1,epair2,Delta,Deltasqp1,     ...
                    alphaL,alphaU,iter,maxiter,iterate,epsilon,message_level);


       while (~ convergence_cond)
%
%  Saves previous interpolation point; chooses new interpolation point
%
          ipointa = ipointb;

          ipointb = inter_point(epair1,epair2,iterate,Balpha.g,message_level);
   
          [alphaL,alphaU] = upd_alpha_safe(alpha,alphaL,alphaU,iterate, ...
                                           ipointb,Delta,message_level);

          alpha = upd_paramk(ipointa,ipointb,Delta,alphaL,alphaU, ...
                             deltaU,message_level);

          Balpha.alpha = alpha;

          if (isfield(eigensolverpar,'tol') & ~isempty(lopts.maxeigentol))
             eigensolverpar.tol = ...
             adjust_eigentol(eigensolverpar.tol,iterate,epair1,epair2,Delta, ...
                             lopts.maxeigentol,message_level);
          end

          [nconv,epair1,epair2,v1] = ...
           b_epairs(Balpha,eigensolver,eigensolverpar,lopts);
          eigensolverpar.v0 = v1;

          adj = adj + 1;

          if (nconv > 0)
             deltaU = upd_deltaU(epair1,Balpha.g,deltaU,message_level);
          end

          [alpha,alphaU,nconv,epair1,epair2,z,deltaU,adjk,iterate,v1] =       ...
           adjust_alpha(nconv,epair1,epair2,Balpha,norg,alphaL,alphaU,deltaU, ...
                        epsilon,z,eigensolver,eigensolverpar,lopts);

          Balpha.alpha = alpha;

          adj  = adj + adjk;
          iter = iter + 1;

          [convergence_cond,xtilde,lambdatilde,stopcond] =                ...
           convergence(Balpha,nconv,epair1,epair2,Delta,Deltasqp1,        ...
                       alphaL,alphaU,iter,maxiter,iterate,epsilon,message_level);

       end

    end

    iter = iter + 1;

    [x,lambda,info,moreinfo] =                                     ...
     output(epair1,epair2,iterate,lambdatilde,xtilde,Delta,       ...
            Balpha,norg,z,iter,adj,stopcond,lopts);


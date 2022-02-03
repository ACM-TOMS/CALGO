%
%   File:              lstrs.m
%
%   Description:       Interface for lstrs_method which implements LSTRS:
%                      A New Matrix-Free Algorithm for the
%                      Large-Scale Trust-Region Subproblem
%
%                      min   1/2 x'Hx + g'x
%                      s.t.  ||x|| <= Delta
%
%                      H real, square, symmetric matrix; g real vector; 
%                      Delta positive scalar
%
%                      Based on: M. Rojas, S.A. Santos and D.C. Sorensen.
%                                A New Matrix-Free Algorithm for the Large-Scale
%                                Trust-Region Subproblem, SIAM J. Optim.,
%                                11(3):611-646, 2000.
%
%                      The following special case of this problem is used
%                      for Regularization of discrete forms of ill-posed problems:
%
%                      min   1/2 ||Ax - b||^2
%                      s.t.  ||x|| <= Delta
%
%                      A real, rectangular matrix; b real vector; 
%                      Delta positive scalar
%
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
%   Functions called:  clf, disp, error, exist, figure, func2str, isa,
%                      isempty, isstruct, legend, length, lstrs_method, lower,
%                      plot, sprintf, strcmp, title, warning
%
%   General Call: [X,LAMBDA,INFO,MOREINFO] =  ...
%                 lstrs(H,G,DELTA,EPSILON,EIGENSOLVER,LOPTS,HPAR,EIGENSOLVERPAR);
%
%   Output Parameters:
%   X:         real, nx1 vector, the solution to the trust-region subproblem
%
%   LAMBDA:    nonpositive scalar, the Lagrange multiplier
%
%   INFO:      the result of the computation 
%              0: Boundary solution
%              1: Interior solution
%              2: Quasi-Optimal solution
%             -1: An interior solution was detected for a regularization problem,
%                 but the linear system was not solved. The last iterate is returned.
%             -2: The safeguarding interval cannot be further reduced.
%                 The last iterate is returned.
%             -3: The maximum number of iterations was reached. The last iterate
%                 is returned.
%             -4: No iterate could be computed: x is empty. 
%
%   MOREINFO: a structure with the following fields
%             EXITCOND: an array of strings describing all the stopping criteria
%                       satisfied. Possible values:
%                       BS:  Boundary Solution
%                       IS:  Interior Solution
%                       QO:  Quasi-Optimal Solution
%                       AI:  Small safeguarding interval 
%                       IT:  Maximum number of iterations reached 
%             MVP:      number of matrix-vector products 
%             ITER:     number of LSTRS iterations
%             SOLVES:   number of calls to the eigensolver 
%             KKT:      The value of ||(H-\lambdaI)x + g||/||g||
%                       or -1 if no iterate could be computed
%             ALPHA:    The last value of the scalar parameter alpha
%
%   Required Input Parameters:
%   
%   H:               real, nxn, symmetric matrix, or
%                    string or function-handle specifying a
%                    matrix-vector multiplication routine 
%
%   See 1. below for calling specifications in case H is a routine.
%
%   G:               real, nx1 vector
%
%   DELTA:           positive scalar (trust-region radius)
%
%   Optional Input Parameters:
%
%
%   EPSILON:         structure containing tolerances for the method LSTRS
%                    (default values in parentheses):
%     epsilon.Delta: accuracy of the norm of the solution with respect
%                    to Delta (1e-4)
%     epsilon.HC:    accuracy in the Hard Case (1e-4)
%     epsilon.Int:   threshold for Interior solutions (1e-8)
%     epsilon.alpha: threshold for safeguarding interval for
%                    parameter alpha (1e-8)
%     epsilon.nu:    threshold for small components (1e-2)
%
%   EIGENSOLVER:     string or function handle specifying the eigensolver routine.
%
%   Current choices:
%   'eig_gateway', 'eigs_lstrs_gateway' (ARPACK),
%   'tcheigs_lstrs_gateway'(ARPACK+Tchebyshev Spectral Transformation),
%   user-provided
%
%   Default: 'eigs_lstrs_gateway'
%
%   See 2. below for calling specifications.
%
%
%   LOPTS:           a structure specifying choices for optional parameters,
%                    with fields:
%
%   lopts.maxiter:        maximum number of LSTRS iterations allowed. Default: 50.
%
%   lopts.message_level:  possible values are
%                         0 - no messages
%                         1 - a message per iteration plus a final summary
%                         2 - more detailed messages
%                         Default: 1
%
%   lopts.name:           string containing the name of the problem to solve
%
%   lopts.plot:           string indicating if a plot of the solution is desired.
%                         Possible values: a string starting with 'y' or 'Y' (plot),
%                         other string (no plot). Default: no plot.
%
%   lopts.correction:     string indicating if, in the hard case, a correction term
%                         in the direction of an appropriate eigenvector should be
%                         added. Adding the correction term might not be desirable
%                         in some cases, such as in regularization problems.
%                         Possible values: a string starting with 'y' or 'Y' (add),
%                         other string (do not add). Default: add term.
%
%   lopts.interior:       string indicating if, when the existence of an interior
%                         solution is detected, such solution should be computed.
%                         Note that, in regularization problems, an interior
%                         solution is usually of no interest. Possible values:
%                         a string starting with 'y' or 'Y' (compute), other string
%                         (do not compute). Default: compute interior solution.
%
%   lopts.intsoltol:      a scalar indicating the tolerance for the residual, i.e.
%                         the stopping criterion for the linear system solver
%                         when computing an interior solution. Default: epsilon.Delta
%
%   lopts.deltaU:         a string indicating how to initialize deltaU, an upper bound
%                         for the smaller eigenvalue of H, or a scalar with the initial
%                         value. Possible values: 'rayleigh', 'mindiag', any scalar.
%                         Default: 'rayleigh'.
%
%   lopts.alpha:          a string indicating how to initialize alpha, or a scalar
%                         with the initial value. Possible values: 'min', 'deltaU',
%                         any scalar. Default: 'min'.
% 
%   lopts.maxeigentol:    a scalar indicating the maximum relative accuracy in the
%                         eigenpairs, in case the user wants to adjust this accuracy
%                         at each iteration. It could also be a structure containing
%                         the maximum relative accuracy in the eigenpairs (maxeigentol)
%                         and the accuracy in the current iterate (itermaxacc). Two different
%                         adjustment strategies are implemented in 'adjust_eigentol.m'
%                   
%                         Default: '[]' (no adjustment, i.e. fixed eigenpair accuracy).
%
%   lopts.heuristics:     scalar indicating if a heuristics should be used to compute
%                         iterates. This option can only be used in combination with
%                         the eigensolver eigs_lstrs_gateway.
%                         Possible values: any scalar. Default: 0.
%
%
%   HPAR:            structure containing parameters for H
%
%   EIGENSOLVERPAR:  structure containing parameters for 
%                    eigensolver routine EIGENSOLVER. See 2. below.
%
%   NOTE: All optional parameters can be given the value `[]'
%         in which case a default value is used, if appropriate.
%         The order in which parameters of the type struct appear
%         indicates which one takes on the default value, according to:
%         epsilon, eigensolver, lopts, Hpar, eigensolverpar.
%
%
%   1. Calling Specifications for H (routine):
%
%   If H is a matrix-vector multiplication routine that computes
%      W = <Hessian matrix>*V 
%   H is called as:
%      W = H(V,HPAR) where HPAR is a structure containing additional
%                    parameters for H
%
%   2. Calling Specifications for EIGENSOLVER:
%
%      [NCONV,LAMBDA1,Y1,LAMBDA2,Y2,V1] = ...
%      EIGENSOLVER(H,alpha,g,Hpar,eigensolverpar,message_level);
%
%   Output Parameters:
%   NCONV:          number of converged eigenvalues
%   LAMBDA1, Y1:    smallest eigenvalue and corresponding eigenvector
%   LAMBDA2, Y2:    another eigenvalue and corresponding eigenvector
%                   (in general, the second smallest eigenvalue or close to it)
%   V1:             First column of Lanczos basis matrix (or first eigenvector)
%
%   Input Parameters:
%   H:              Hessian matrix or matrix-vector multiplication routine
%   ALPHA:          key scalar parameter in LSTRS
%   G:              gradient vector in objective function 
%   HPAR:           parameters for matrix-vector multiplication routine H
%                   See 1. above
%   EIGENSOLVERPAR: structure containing parameters for eigensolver routine
%   MESSAGE_LEVEL:  as above
%  

function [x,lambda,info,moreinfo] = lstrs(H,g,Delta,varargin)

%
% Next lines of code correspond to the interface: 
% processing the input arguments and setting some default values
%

%
% Check type of required parameters
%

if (~isa(H,'double') & ~isa(H,'char') & ~isa(H,'function_handle'))
   error('LSTRS: Type of H must be double, string or function-handle');
elseif (~isa(g,'double'))
   error('LSTRS: Type of g must be double');
elseif (length(g) < 2)
   error('LSTRS: g must be a vector of dimension at least 2');
elseif (~isa(Delta,'double'))
   error('LSTRS: Type of Delta must be double');
elseif (length(Delta) ~= 1)
   error('LSTRS: Delta must be a scalar');
end

%
% Check that the Trust-Region Subproblem is well-defined and nontrivial
%
if (~any(g))
   if (Delta > 0)
      disp(sprintf('Gradient g is zero.'));
      disp(sprintf('Solution is an eigenvector corresponding to'));
      disp(sprintf('the smallest eigenvalue of H, normalized so that ||x||=Delta.\n'));
      x=[]; lambda = []; info = []; moreinfo = [];
      return
   end
else
   n = length(g);
end

if (Delta < 0) 
   error('Delta must be greater than zero');
elseif (Delta == 0) 
   disp(sprintf('Delta = 0. Solution is x = 0.\n'));
   x=0; lambda = 0; info = 0; moreinfo = [];
   return
end

epsflag  = 0; eigensolverflag = 0; loptsflag = 0;
Hparflag = 0; eparflag        = 0;

if (nargin > 8)
   error('LSTRS: Too many input arguments!')
elseif (nargin < 3)
   error('LSTRS: Not enough input arguments!')
else
   if (length(varargin) > 5)
      error('LSTRS: Too many Optional input arguments!')
   else
      for i=1:length(varargin)
          if (isempty(varargin{i}))
             if (~epsflag)
                epsilon = [];
                epsflag = 1;
             elseif (~eigensolverflag)
               eigensolver     = 'eigs_lstrs_gateway';
               eigensolverflag = 1;
             elseif (~loptsflag)
                lopts     = [];
                loptsflag = 1;
             elseif (~Hparflag)
                Hpar     = [];
                Hparflag = 1;
             elseif (~eparflag)
                eigensolverpar = [];
                eparflag       = 1;
             end
          elseif (isa(varargin{i},'function_handle') & ~eigensolverflag)
             eigensolver     = varargin{i};
             eigensolverflag = 1;
          elseif (isa(varargin{i},'char') & ~eigensolverflag)
             eigensolver     = varargin{i};
             eigensolverflag = 1;
          elseif (isstruct(varargin{i}))
             if (~epsflag)
                epsilon = varargin{i};
                epsflag = 1;
             elseif (~loptsflag)
                lopts     = varargin{i};
                loptsflag = 1;
             elseif (~Hparflag)
                Hpar     = varargin{i};
                Hparflag = 1;
             else
                eigensolverpar = varargin{i};
                eparflag       = 1;
             end
          else
             error('LSTRS: Invalid Type for Input Argument');
          end
      end
   end
end

if (~epsflag)         epsilon        = [];                   end
if (~eigensolverflag) eigensolver    = 'eigs_lstrs_gateway'; end
if (~loptsflag)       lopts          = [];                   end
if (~Hparflag)        Hpar           = [];                   end
if (~eparflag)        eigensolverpar = [];                   end

%
% Set default values for epsilon
%
if (isempty(epsilon))
   epsilon.Delta = 1e-4;
   epsilon.Int   = 1e-8;
   epsilon.HC    = 1e-4;
   epsilon.nu    = 1e-2;
   epsilon.alpha = 1e-8;
else
   if (~isfield(epsilon,'Delta') | isempty(epsilon.Delta))
      epsilon.Delta = 1e-4; 
   end
   if (~isfield(epsilon,'Int') | isempty(epsilon.Int))
      epsilon.Int   = 1e-8;
   end
   if (~isfield(epsilon,'HC') | isempty(epsilon.HC))
      epsilon.HC    = 1e-4;
   end
   if (~isfield(epsilon,'nu') | isempty(epsilon.nu))
      epsilon.nu    = 1e-2;
   end
   if (~isfield(epsilon,'alpha') | isempty(epsilon.alpha))
      epsilon.alpha = 1e-8;
   end
end


%
% Set default values for eigensolver
%
 
if (isa(eigensolver,'char'))
   sesolver = eigensolver;
else
   sesolver = func2str(eigensolver);
end

iseig_gateway    = strcmp(sesolver,'eig_gateway');
isarpack_gateway = strcmp(sesolver,'eigs_lstrs_gateway') | ...
                   strcmp(sesolver,'tcheigs_lstrs_gateway');

if (isempty(eigensolverpar))
   if (iseig_gateway)
        eigensolverpar = [];
   elseif (isarpack_gateway)
        eigensolverpar.tol   = 1e-2;
        eigensolverpar.maxit = 13;
        eigensolverpar.k     = 2;
        if (n < 7) 
            eigensolverpar.p = n+1;
        else
            eigensolverpar.p = 7;
        end
        eigensolverpar.issym = 1;
        eigensolverpar.disp  = 0;
   end
elseif (isarpack_gateway)
      if (~isfield(eigensolverpar,'tol'))   eigensolverpar.tol   = 1e-2; end
      if (~isfield(eigensolverpar,'maxit')) eigensolverpar.maxit = 13;   end
      if (isfield(eigensolverpar,'k'))
         if (eigensolverpar.k < 2)
             eigensolverpar.k = 2;
         end
      else
         eigensolverpar.k = 2;
      end
      if (n < 7)
         eigensolverpar.p = n+1;
      elseif (~isfield(eigensolverpar,'p'))
         eigensolverpar.p = 7;
      end
      if (~isfield(eigensolverpar,'issym'))
         eigensolverpar.issym = 1;
      elseif (~eigensolverpar.issym)
         emessage = sprintf('LSTRS: opts.issym = 0 for eigs_lstrs !');
         emessage = strcat(emessage,' The Hessian matrix must be symmetric !');
         error(emessage)
      end
      if (~isfield(eigensolverpar,'disp'))  eigensolverpar.disp = 0;  end
end

if (iseig_gateway & ~isa(H,'double'))
   error('LSTRS: To use the eigensolver ''eig_gateway'', ''H'' must be a matrix !');
end 

if (~exist(sesolver))
   error(sprintf('LSTRS: Undefined eigensolver: ''%s''. Not in search path.',sesolver));
end

%
% Set default values for lopts
%
if (isempty(lopts))
   lopts.message_level = 1;
   lopts.heuristics    = 0;
   lopts.correction    = 'y';
   lopts.interior      = 'y';
   lopts.intsoltol     = epsilon.Delta;
   lopts.plot          = 'n';
   lopts.name          = 'no name available';
   lopts.maxiter       = 50;
   lopts.deltaU        = 'rayleigh';
   lopts.alpha         = 'min';
   lopts.maxeigentol   = [];
else
   if (isfield(lopts,'message_level'))
      if (~isa(lopts.message_level,'double') | (length(lopts.message_level) ~= 1))
         error('LSTRS: input parameter message_level must be a scalar');
      else
         if (lopts.message_level ~= 0 & lopts.message_level ~= 1 & ...
             lopts.message_level ~= 2)
            error('LSTRS: input parameter lopts.message_level must be 0, 1, or 2');
         end
      end
   else  lopts.message_level = 1; end

   if (isfield(lopts,'heuristics'))
      if (~isa(lopts.heuristics,'double') | (length(lopts.heuristics) ~= 1))
         error('LSTRS: input parameter lopts.heuristics must be a scalar');
      else
         if (~strcmp(sesolver,'eigs_lstrs_gateway'))
            error('LSTRS: Heuristics can only be used with eigs_lstrs_gateway');
         end
      end
   else
      lopts.heuristics = 0;
   end

   if (isfield(lopts,'correction'))
      if (~isa(lopts.correction,'char'))
         error('LSTRS: input parameter lopts.correction must be a string');
      end
   else  lopts.correction = 'y'; end

   if (isfield(lopts,'interior'))
      if (~isa(lopts.interior,'char'))
         error('LSTRS: input parameter lopts.interior must be a string');
      end
   else  lopts.interior = 'y'; end

   if (isfield(lopts,'intsoltol'))
      if (~isa(lopts.intsoltol,'double') | (length(lopts.intsoltol) ~= 1))
         error('LSTRS: input parameter lopts.intsoltol must be a scalar');
      end
   else  lopts.intsoltol = epsilon.Delta; end
         
   if (isfield(lopts,'plot'))
      if (~isa(lopts.plot,'char'))
         error('LSTRS: input parameter lopts.plot must be a string');
      end
   else lopts.plot = 'n'; end

   if (isfield(lopts,'name'))
      if (~isa(lopts.name,'char'))
         error('LSTRS: input parameter lopts.name must be a string');
      end
   else  lopts.name = 'no name available'; end

   if (isfield(lopts,'maxiter'))
      if (~isa(lopts.maxiter,'double') | ...
          (length(lopts.maxiter) > 1) | (lopts.maxiter < 1))
         error('LSTRS: input parameter lopts.maxiter must be a positive integer');
      end
   else  lopts.maxiter = 50; end

   if (isfield(lopts,'deltaU'))
      if (isa(lopts.deltaU,'char'))
         if (~strcmp(lopts.deltaU,'rayleigh') & ~strcmp(lopts.deltaU,'mindiag'))
            error('LSTRS: input parameter lopts.deltaU valid string values: ''rayleigh'', ''mindiag''');
         end
         if (strcmp(lopts.deltaU,'mindiag') & ~isa(H,'double'))
            error('LSTRS: lopts.deltaU must be different from ''mindiag'' if H is NOT an array');
         end
      elseif (~isa(lopts.deltaU,'double') | ...
             (isa(lopts.deltaU,'double') & (length(lopts.deltaU) > 1)))
         error('LSTRS: input parameter lopts.deltaU must be a string or a scalar');
      end
   else  lopts.deltaU = 'rayleigh'; end

   if (isfield(lopts,'alpha'))
      if (isa(lopts.alpha,'char'))
         if (~strcmp(lopts.alpha,'min') & ~strcmp(lopts.alpha,'deltaU'))
            error('LSTRS: input parameter lopts.alpha valid string values: ''min'', ''deltaU''');
         end
      elseif (~isa(lopts.alpha,'double') | ...
             (isa(lopts.alpha,'double') & (length(lopts.alpha) > 1)))
         error('LSTRS: input parameter lopts.alpha must be a string or a scalar');
      end
   else  lopts.alpha = 'min'; end

   if (isfield(lopts,'maxeigentol'))
      if (~isempty(lopts.maxeigentol) & ~isa(lopts.maxeigentol,'double') & ...
                                        ~isstruct(lopts.maxeigentol)     | ...
             (isa(lopts.maxeigentol,'double') & (length(lopts.maxeigentol) > 1)))
         error('LSTRS: input parameter lopts.maxeigentol must be a scalar, a structure or empty ([])');
      end
   else  lopts.maxeigentol = []; end

end

%
% Displays problem information
%

if (lopts.message_level > 0)
   disp(sprintf('\nProblem: %s.     Dimension: %d.    Delta: %e\n', ...
                lopts.name,length(g),Delta));
   disp(sprintf('Eigensolver: %s\n',sesolver));
end

%
% Actual call to LSTRS, the Trust-Region Solver
%

[x,lambda,info,moreinfo] = ...
 lstrs_method(H,g,Delta,epsilon,eigensolver,lopts,Hpar,eigensolverpar);


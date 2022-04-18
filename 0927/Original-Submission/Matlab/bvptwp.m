function sol=bvptwp(f,g,solinit,options)
% bvptwp solve boundary value problems for ODEs by using a
% deferred correction scheme based respectively on Mono-Implicit Runge-Kutta
% (MIRK) (options:solver = 'twpbvp_m', 'twpbvpc_m'), on LOBATTO formulas 
% (options: solver = 'twpbvp_l', 'twpbvpc_l').
% It allows also LOBATTO formulas with automatic continuation 
% (options: solver = 'acdc', 'acdcc').
%
% SOL = bvptwp(f, g, solinit) integrates a system of ordinary differential
% equations of the form y' = f(x,y) on the interval [a,b], subject to
% general two-point boundary conditions of the form g(y(a),y(b)) = 0.
% f and g are function handles. For a scalar X and a column vector Y, f(X,Y)
% must return a column vector representing f(x,y). For column vectors YA and YB,
% g(YA,YB) must return a column  vector representing g(y(a),y(b)).
% SOLINIT is a structure with fields
%         x -- ordered nodes of the initial mesh with
%              SOLINIT.x(1) = a, SOLINIT.x(end) = b
%         y -- initial guess for the solution with SOLINIT.y(:,i)
%              a guess for y(x(i)), the solution at the node SOLINIT.x(i)
%         fixpnt -- fixed points used in all meshes.
%
% The output SOL is a structure with
%         SOL.x  -- last mesh selected by bvptwp
%         SOL.y  -- approximation to y(x) at the mesh points of SOL.x
%         SOL.solver -- 'twpbvp_m' or 'twpbvp_l' or 'acdc' or 
%                       'twpbvpc_m' or 'twpbvpc_l' or 'acdcc'
%         SOL.lambda  -- the final value of the continuation 
%                            parameter used (only for acdc, acdcc)
%         SOL.iflbvp  -- 0, the code solved the problem
%                     -- -1, the code solved a problem with a different
%                            continuation parameter  (only for acdc, acdcc)
%                     -- 1, tha code failed (maximum number of meshpoints
%                           reached)
%                     -- 2, the code failed (maximum number of possible
%                           meshes reached, default 100, only for twpbvpc_m
%                           twpbvpc_l)
%                     -- 3, the code failed (maxmum number of continuation
%                           steps reached, only for acdc, acdcc)
%                     -- 4, the code failed (unknown error)
%         SOL.condpar -- information about the conditioning parameters 
%                        (only for twpbvpc_m, twpbvpc_l, acdcc)
%                         SOL.condpar.kappa : conditioning of the bvp 
%                                  in inf norm
%                         SOL.condpar.kappa1 : conditioning related to
%                              changes in the initial values (inf norm)
%                         SOL.condpar.kappa2 : conditioning related to
%                              the Green's function (inf norm)
%                         SOL.condpar.gamma1 : conditioning related to
%                              changes in the initial values (1 norm)
%                         SOL.condpar.sigma : stiffness parameter
%                         SOL.condpar.stabcond : 1 or 0 if the conditioning
%                              parameters stabilized
%         SOL.stats  -- Computational cost statistics and information about
%                       the maximum scaled error computed
%
% SOL = bvptwp(f, g, solinit, options) solves as above with default parameters
% replaced by values in OPTIONS, a structure created with the bvptwpset function.
% To reduce the run time greatly, use OPTIONS to supply a function for
% evaluating the Jacobian and/or Jvectorize.
%
%  The function bvpinit forms the guess structure in the most common situations:
%  SOLINIT = bvpinit(X,YINIT) forms the guess for an initial mesh X as described
%  for SOLINIT.x, and YINIT either a constant vector guess for the solution or
%  a function handle. If YINIT is a function handle then for a scalar X, YINIT(X)
%  must return a column vector, a guess for the solution at point x in [a,b].
%  
%  Fixed points in the row vector fixedpoints could be added to SOLINIT using
%        SOLINIT.fixpnt = fixedpoints
%
%   Example
%shock_bvptwp  The solution has a shock layer near x = 0
%   This is an example used in U. Ascher, R. Mattheij, and R. Russell,
%   Numerical Solution of Boundary Value Problems for Ordinary Differential
%   Equations, SIAM, Philadelphia, PA, 1995,  to illustrate the mesh
%   selection strategy of COLSYS.
%
%   For 0 < e << 1, the solution of
%
%       e*y'' + x*y' = -e*pi^2*cos(pi*x) - pi*x*sin(pi*x)
%
%   on the interval [-1,1] with boundary conditions y(-1) = -2 and y(1) = 0
%   has a rapid transition layer at x = 0.
%
%   For this problem,
%   analytical partial derivatives are easy to derive and the solver benefits
%   from using them.
%
%   By default, this example uses the twpbvpc_l.
%   Use syntax 
%   SHOCK_BVPTWP(solver) to solve this problem with the another solver
%     available solvers are: 'twpbvp_m', 'twpbvpc_m', 'twpbvp_l', 'twpbvpc_l',
%                            'acdc', 'acdcc'
%
%   For more examples see bvptwp, bvptwpset, bvptwpget, bvpMtest.
%   This code is based on the Fortran codes TWPBVP, TWPBVPC, TWPBVPLC, ACDC
%
%   For a full description of the input parameters see the report manual_bvptwp.pdf
%
%   
%
%       Authors:
%
%       Jeff R. Cash 
%            (Department of Mathematics, Imperial College,  London, England.)
%       Davy  Hollevoet 
%            (Vakgroep Toegepaste Wiskunde en Informatica, Universiteit Gent, Belgium.)
%       Francesca Mazzia  
%            (Dipartimento di Matematica, Universita' di Bari, Italy)
%       Abdelhameed Nagy Abdo
%            (Dipartimento di Matematica, Universit\`a di Bari, Italy)
%            (Dept. of Mathematics, Faculty of Sciences, Benha  University,Egypt)
%            
%


problem.NFUN = 0;   
problem.NBC = 0;
if nargin<3
    error('MATLAB:bvptwp:NotEnoughInputs','Need at least f, g and solinit');
end
if ~isa(f,'function_handle')
    error('MATLAB:bvptwp:FNoFunction','First argument should be a function handle');
end
if ~isa(g,'function_handle')
    error('MATLAB:bvptwp:DFNoFunction','Second argument should be a function handle');
end

if ~isstruct(solinit)
    error('MATLAB:bvptwp:SolinitNotStruct','Third argument should be a struct')
end
if ~isfield(solinit,'x')
    error('MATLAB:bvptwp:NoXInSolinit','Field ''x'' not present in initial guess');
end
if ~isfield(solinit,'y')
    error('MATLAB:bvptwp:NoYInSolinit','Field ''y'' not present in initial guess');
end

if nargin<4
    options=bvptwpset();
end

problem.linear=~isempty(bvptwpget(options,'Linear')) && strcmp(options.Linear,'on');

if ~isempty(bvptwpget(options,'Solver'))
    problem.solver=options.Solver;
    
else
    if problem.linear
      problem.solver='twpbvp_l';
    else
      problem.solver='twpbvp_m';
    end
end
if strcmp(problem.solver,'twpbvp_m') ||strcmp(problem.solver,'twpbvp_l') || strcmp(problem.solver,'acdc') 
     problem.conditioning=0;
elseif  strcmp(problem.solver,'twpbvpc_m')||strcmp(problem.solver,'twpbvpc_l') ||strcmp(problem.solver,'acdcc') 
    problem.conditioning=1;
else
      error('MATLAB:bvptwp:solver','Invalid solver %s', problem.solver);
end


if ~isempty(bvptwpget(options,'LambdaMin'))
    lambdamin= bvptwpget(options,'LambdaMin');
else
    lambdamin=1e-5;
end

if ~isempty(bvptwpget(options,'LambdaStart'))
    problem.lambda = bvptwpget(options,'LambdaStart');
else
    problem.lambda =  0.5;
end

if ~isempty(bvptwpget(options,'FJacobian')) && ~isa(options.FJacobian,'function_handle')
    error('MATLAB:bvptwp:FJacNoFunction','options.FJacobian is not a function handle');
end

if ~isempty(bvptwpget(options,'BCJacobian')) && ~isa(options.BCJacobian,'function_handle')
    error('MATLAB:bvptwp:BCJacNoFunction','options.BCJacobian is not a function handle');
end

if strcmp('on',bvptwpget(options,'Debug'))
    problem.debug=true;
else
    problem.debug=false;
end

if strcmp('on',bvptwpget(options,'Vectorized'))
    problem.f=f;
    problem.vectorized =1;
elseif strcmp(problem.solver,'acdc') ||strcmp(problem.solver,'acdcc')
    problem.tscalarf=f;
    problem.f=@fvectorizer_a;
    problem.vectorized =0;
else
    problem.tscalarf=f;
    problem.f=@fvectorizer;
    problem.vectorized =0;
end

if ~isempty(bvptwpget(options,'FJacobian'))
    if strcmp('on',bvptwpget(options,'JVectorized'))
        problem.df=options.FJacobian;
    elseif  strcmp(problem.solver,'acdc') ||strcmp(problem.solver,'acdcc')
        problem.tscalardf=options.FJacobian;
        problem.df=@dfvectorizer_a;
    else
        problem.tscalardf=options.FJacobian;
        problem.df=@dfvectorizer;
    end
else
    if  strcmp(problem.solver,'acdc') ||strcmp(problem.solver,'acdcc')
        problem.df=@dfnumerical_a;
    else
        problem.df=@dfnumerical;
    end
end

if ~isempty(bvptwpget(options,'BCJacobian'))
    problem.dg=options.BCJacobian;
else
    if strcmp(problem.solver,'acdc') ||strcmp(problem.solver,'acdcc')
        problem.dg=@dbcnumerical_a;
    else
        problem.dg=@dbcnumerical;
    end
end

problem.g=g;

problem.a=solinit.x(1,1);
problem.b=solinit.x(1,end);
problem.ncomp=size(solinit.y,1);
problem.y0=solinit.y;
problem.x0=solinit.x;


if isfield(solinit,'fixpnt')
    if isempty(solinit.fixpnt) || (isnumeric(solinit.fixpnt) && isvector(solinit.fixpnt) && size(solinit.fixpnt,1)==1)
        problem.fixpnt=sort(solinit.fixpnt);
        xn = sort( [solinit.x solinit.fixpnt]);
        xn = xn(find(diff(xn)));
        xn = [xn,problem.b];
        problem.y0=interpu(solinit.x,solinit.y,xn);
        problem.x0=xn;
    else
        error('MATLAB:bvptwp:SolinitInvalidFixedPoints','solinit.fixpnt should be a rowvector')
    end
else
    problem.fixpnt=[];
end

tol=bvptwpget(options,'RelTol',1e-3);
nmax=bvptwpget(options,'Nmax',floor(50000/problem.ncomp));


ltol=(1:problem.ncomp)';
if isscalar(tol)
    tol=repmat(tol,problem.ncomp,1);
else
    tol=tol(:);
    if size(tol,1)==problem.ncomp
        ltol(tol==0)=[];
        tol(tol==0)=[];
    else
        ltol=ltol(1:size(tol,1));
    end
end



liseries= bvptwpget(options,'MaxNumberOfMeshes',100);
liseries=fix(abs(liseries));

maxcon= bvptwpget(options,'MaxNumberOfContStep',150);
maxcon=fix(abs(maxcon));

odenumjacproxy=@checkodenumjac;


warning off MATLAB:nearlySingularMatrix


if strcmp(problem.solver,'twpbvpc_m')||strcmp(problem.solver,'twpbvp_m') 
    [y,t,err,success,maxmsh,stabcond,condpar,problem,iseries]=twpbvp_m(problem,problem.x0,nmax,ltol,tol,liseries);
    sol=struct('x',t,'y',y,'solver',problem.solver);
    sol.stats.iseries = sparse(iseries);
    sol.stats.maxmesh = max(iseries);
    if success == 1
       sol.iflbvp=~success;
    else
       if maxmsh 
         sol.iflbvp=1;
       else
         sol.iflbvp=success;
       end  
    end
elseif strcmp(problem.solver,'twpbvpc_l')||strcmp(problem.solver,'twpbvp_l') 
    [y,t,err,success,maxmsh,stabcond,condpar,problem,iseries]=twpbvp_l(problem,problem.x0,nmax,ltol,tol,liseries);
    sol=struct('x',t,'y',y,'solver',problem.solver);
    sol.stats.iseries = sparse(iseries);
    sol.stats.maxmesh = max(iseries);
    
    if success == 1
       sol.iflbvp=~success;
    else
       if maxmsh 
         sol.iflbvp=1;
       else
         sol.iflbvp=success;
       end  
    end
else
    problem.iflbvp = 4;
    problem.iprec = 0;
    problem.iback = 0;
    problem.ifinal = 0;
    [y,t,err,maxmsh,stabcond,condpar,problem,nc,nmaxmesh] = acdc(problem,problem.x0,nmax,problem.lambda,lambdamin,ltol,tol,maxcon);
    sol=struct('x',t,'y',y,'solver',problem.solver);
    sol.lambda= problem.lambda;
    sol.iflbvp=problem.iflbvp;
       if maxmsh 
         sol.iflbvp=1;
       elseif nc == maxcon && sol.iflbvp~=0
         sol.iflbvp=3;
       end  
end
if problem.conditioning==1;
    sol.condpar = condpar;
    sol.condpar.stabcond = stabcond;
end
if strcmp(problem.solver,'acdc')||strcmp(problem.solver,'acdcc')
    sol.stats.nc = nc;
    sol.stats.maxmesh=nmaxmesh;
    if (problem.iflbvp ~= 0)
        if nc == maxcon
            warning('MATLAB:bvptwp:maxnc','Epsmin changed, too many continuation steps');
        end
        if maxmsh
            error('MATLAB:bvptwp:maxmsh','Unable to meet the tolerance without using more than %d mesh points.',nmax);
        elseif problem.iflbvp==-1
            warning('MATLAB:bvptwp:epschanged','Epsmin changed, successfull termination');
        else
            error('MATLAB:bvptwp:failedunknown','Unknown error: bvptwp failed');
        end
    end    
else
    if ~(success==1)
        if maxmsh
            error('MATLAB:bvptwp:maxmsh','Unable to meet the tolerance without using more than %d mesh points.',nmax);
        elseif success==2
            error('MATLAB:bvptwp:maxnumberofmsh','Unable to meet the tolerance without using more than %d meshes.',liseries);
        else
            error('MATLAB:bvptwp:failedunknown','Unknown error: bvptwp failed');
        end
    end
end

printstats = bvptwpget(options,'Stats','off');
printstats = strcmp(printstats,'on');

sol.stats.nmeshpoints = length(t);
sol.stats.maxScaledError = err;
sol.stats.nODEevals =  problem.NFUN;
sol.stats.nBCevals = problem.NBC;
if printstats
    fprintf('-------------------------------------\n');
    fprintf('Solver %s.\n', problem.solver);
    fprintf('The solution was obtained on a mesh of %g points.\n',length(t));
    fprintf('The maximum scaled error  is %10.3e. \n',err);
    fprintf('There were %g calls to the ODE function. \n',problem.NFUN);
    fprintf('There were %g calls to the BC function. \n',problem.NBC);
end

    function fty=fvectorizer_a(t,y,lambda)
        sf=problem.tscalarf
        fty=zeros(size(y));
        for ti=1:size(t,2)
            fty(:,ti)=sf(t(1,ti),y(:,ti),lambda);
        end
    end

    function dfdy=dfvectorizer_a(t,y,lambda)
        sdf=problem.tscalardf;
        dfdy=zeros(problem.ncomp,problem.ncomp,size(t,2));
        for ti=1:size(t,2)
            dfdy(:,:,ti)=sdf(t(1,ti),y(:,ti),lambda);
        end
    end

    function dfdy=dfnumerical_a(t,y,lambda)
        
        dfdy=zeros(problem.ncomp,problem.ncomp,size(t,2));
        %todo: user-settable threshold or smth
        opts=struct('diffvar',2,'vectvars',[1,2],'thresh',1e-10,'fac',[]);
        
        for ti=1:size(t,2)
            dfdy(:,:,ti)=odenumjacproxy(problem.f,{t(1,ti),y(:,ti),lambda},problem.f(t(1,ti),y(:,ti),lambda),opts);
        end
    end

    function [dga,dgb]=dbcnumerical_a(ya,yb,lambda)
        %todo: user-settable threshold or smth
        opts=struct('diffvar',1,'vectvars',[],'thresh',1e-10,'fac',[]);
        
        bcargs={ya,yb,lambda};
        bcval=problem.g(bcargs{:});
        
        
        dga=odenumjacproxy(problem.g,bcargs,bcval,opts);
        opts.diffvar=2;
        dgb=odenumjacproxy(problem.g,bcargs,bcval,opts);
    end



    function fty=fvectorizer(t,y)
        sf=problem.tscalarf;
        fty=zeros(size(y));
        for ti=1:size(t,2)
            fty(:,ti)=sf(t(1,ti),y(:,ti));
        end
    end

    function dfdy=dfvectorizer(t,y)
        sdf=problem.tscalardf;
        dfdy=zeros(problem.ncomp,problem.ncomp,size(t,2));
        for ti=1:size(t,2)
            dfdy(:,:,ti)=sdf(t(1,ti),y(:,ti));
        end
    end

    function dfdy=dfnumerical(t,y)
        
        dfdy=zeros(problem.ncomp,problem.ncomp,size(t,2));
        %todo: user-settable threshold or smth
        opts=struct('diffvar',2,'vectvars',[1,2],'thresh',1e-10,'fac',[]);
        
        for ti=1:size(t,2)
            
            dfdy(:,:,ti)=odenumjacproxy(problem.f,{t(1,ti),y(:,ti)},problem.f(t(1,ti),y(:,ti)),opts);
        end
    end

    function [dga,dgb]=dbcnumerical(ya,yb)
        %todo: user-settable threshold or smth
        opts=struct('diffvar',1,'vectvars',[],'thresh',1e-10,'fac',[]);
        
        bcargs={ya,yb};
        bcval=problem.g(bcargs{:});
        
        
        dga=odenumjacproxy(problem.g,bcargs,bcval,opts);
        opts.diffvar=2;
        dgb=odenumjacproxy(problem.g,bcargs,bcval,opts);
    end

    function dappr=checkodenumjac(varargin)
        
        % We cannot use odenumjac since it is for private use only. We also
        % cannot include odenumjac in this package because we are not
        % allowed to redistribute it. So we just copy odenumjac.m to
        % somewhere we actually can use it. This is fragile ofcourse...
        
        if strcmp(getfield(functions(@odenumjac),'file'),'')
            
            % do not copy odenumjac to the current directory, but to the
            % directory of e.g. bvpsol, which is a steady name I suppose.
            dstpath=fileparts(getfield(functions(@bvpsol),'file'));
            copyfile(fullfile(matlabroot,'toolbox','matlab','funfun','private','odenumjac.m'),fullfile(dstpath,'odenumjac.m'));
            
            if strcmp(getfield(functions(@odenumjac),'file'),'')
                error('MATLAB:twpbvp:copyodenumjac','Numerical jacobians not available: could not make a local copy of odenumjac.m');
            end
        end
        
        odenumjacproxy=@odenumjac; % from now on, use the copied odenumjac
        dappr=odenumjac(varargin{:});
    end
end

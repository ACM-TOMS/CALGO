function [sol,time] = bvpMtest(problem,solver,tol,plotsol,opts,varargin)
%  bvpMtest solve a problem in the MATLAB bvp testset format using the MATLAB ODE solvers  
%   
% [sol,time] = bvpMtest(problem,method,tol,plotsel,opts,varargin) 
%   
% INPUT  
% problem = handle  of the BVP TEST SET problem  
%     
% solver = string with the calling name of the MATLAB solver, or one of the solvers in bvptwp  
% tol = relative tolerance  
% plotsel = index of component to be plotted (optional)   
% opts = solver parameters set with ODESET or BVPTWPSET (optional)   
%     
% OUTPUT  
% the structure sol given by the solver with
% sol.x  : mesh point of the approximate solution  
% sol.y  : solution computed at sol.x  
% sol.solver : solver used  
% sol.mescd : mixed significant digit of the approximate solution  
% sol.iflbvp  : error flag (0 = successfull run, 1 = error) 
% sol.condpar : information about the conditioning parameters (only for twpbvpc_m, twpbvpc_l, acdcc)
%       sol.condpar.kappa : conditioning of the bvp  in inf norm
%       sol.condpar.kappa1 : conditioning related to
%                            changes in the initial values (inf norm)
%       sol.condpar.kappa2 : conditioning related to
%                              the Green's function (inf norm)
%       sol.condpar.gamma1 : conditioning related to
%                            changes in the initial values (1 norm)
%       sol.condpar.sigma : stiffness parameter
%       sol.condpar.stabcond : 1 or 0 if the conditioning
%                              parameters stabilized.
% sol.stats  : Computational cost statistics and information about
%                       the maximum scaled error computed
%  
%   
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
            


method = solver;
if  strcmp(method,'acdc')||strcmp(method,'acdcc')
    epsmin = varargin{1};
    [prob,f,g,fdy,gdy,esolu,setoutput,settolerances] = problem(varargin{2:end});
else
    [prob,f,g,fdy,gdy,esolu,setoutput,settolerances] = problem(varargin);
end
[problm,type,ncomp,Linear,numjac,numbcjac,Vectorized,JVectorized,solinit] = prob();
tolvec = settolerances(tol);

if nargin<4
    plotsol = [];
end
[solref,printsolout,nindsol,indsol] = setoutput(plotsol);

if ~isempty(opts) || nargin < 5
    opt = opts;
else
    opt = bvptwpset('nmax',10000);
end
opt = bvpset(opt,'Vectorized',Vectorized);
if ~numjac
    opt = bvpset(opt,'FJacobian',fdy);
end
if ~numbcjac
    opt=bvpset(opt,'BCJacobian',gdy);
end


if ~isempty(strfind(method,'twpbvp')) ||  ~isempty(strfind(method,'acdc'))
    opt=bvptwpset(opt,'Solver',method,'reltol',tolvec,'Linear',Linear,'JVectorized',JVectorized);
else
    opt=bvpset(opt,'abstol',tolvec,'reltol',tol);
end
if strfind(method,'acdc')
    opt=bvptwpset(opt,'Lambdamin',epsmin);
end



try
    
    
    if  ~isempty(strfind(method,'twpbvp')) || ~isempty(strfind(method,'acdc'))     
        solver = 'bvptwp';
        tic
        sol = feval(solver,f,g,solinit,opt);
        iflbvp = sol.iflbvp;
        time=toc;
    else
        solver = method;
        tic
        sol = feval(solver,f,g,solinit,opt);
        iflbvp=0;
        time=toc;
    end
    
    if (solref && iflbvp <=0)
        if strcmp(method,'acdc')||strcmp(method,'acdcc')
            [Exact] = esolu(sol.x,sol.lambda);
        else
            [Exact] = esolu(sol.x);
        end
     
      
        merror = 0;    
        for i = 1:ncomp
            merror = max(merror,max(abs(Exact(i,:) - sol.y(i,:))./(1 + abs(Exact(i,:)))));
        end
        sol.error = merror;
        sol.mscd = -log10(merror); 
       
        
    elseif  iflbvp <=0
        opt=bvptwpset(opt,'reltol',tolvec,'Linear',Linear,'JVectorized',JVectorized);
        if strcmp(method,'acdc')||strcmp(method,'acdcc')
            opt=bvptwpset(opt,'LambdaMin',sol.lambda,'LambdaStart',sol.lambda);
        end
        
        solver = 'bvptwp';
        method = 'twpbvpc_l';
        opt=bvptwpset(opt,'Solver',method);
        
        opt.Nmax = opt.Nmax*2;
        fixpnt = (sol.x(1:end-1) + sol.x(2:end))/2;
        solinit1.fixpnt = sort(cat(2,fixpnt,sol.x));
        solinit1.y=sol.y;
        solinit1.x=sol.x;
        sol1 = feval(solver,f,g,solinit1,opt);
        sol2 = interpu(sol1.x,sol1.y,sol.x);
        
     
        merror = 0;        
        for i = 1:ncomp
            merror = max(merror,max(abs(sol2(i,:) - sol.y(i,:))./(1 + abs(sol2(i,:)))));
        end
        sol.error = merror;
        sol.mscd = -log10(merror); 
        
    end
    if printsolout
        for i = 1:nindsol
            figure(i)
            plot(sol.x,sol.y(indsol(i),:))
            title(sprintf('Problem %s, of type %s',problm,type))
        end
    end
catch me
    disp('The program failed to obtain solution')
    time=0;
    sol.error =0;
    sol.iflbvp=1; 
    fprintf('%s: %s',me.identifier, me.message);
    
end

end

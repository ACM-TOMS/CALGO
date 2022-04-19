function [X,Y,feas,stopFlag,objectiveVal] = solve_semi_definite_program(A,b,c,K,reltol,options)
%
% SOLVE_SEMI_DEFINITE_PROGRAM solve the primal SDP min c'x : Ax=b, x in K.
%
% [X,Y,feas,stopFlag,objectiveVal] = SOLVE_SEMIDEFINITE_PROGRAM (A,b,c,K,reltol)
%   Solve the primal problem.
%       min c'x     s.t. Ax = b, x in K
%   With the given relative tolerance.
%
%   Vector X is the primal solution, and vector Y is the dual solution.
%   The output feas (binary) is true when the problem is feasible, and 
%   false else. 
%   The output stopFlag is equal to 0 if all is ok, 1 if the solver is sure
%   that the problem is feasible/unfeasible, but X or Y is inacurate, and
%   equal to 2 if there is a complete failure.
%   Finally, objectiveVal is the value of the objective function.
%
% [...] = SOLVE_SEMIDEFINITE_PROGRAM (A,b,c,K,reltol,options)
%       Does the same as above but with specified parameters described in
%       fields of the structure options. See JSRSETTINGS and below for
%       available parameters and options.
%
%   options.SDPsolver.solver
%       String or function handle. If it is a string, then a
%       default interface already implemented is used. Else, the
%       function handle is used. More details below.
%
%     # Case where solver is a string :
%
%       By default, we solve the SDP with SeDuMi, but it is possible to
%       specify another solver in options.solver.
%       Available interfaces are:
%               - SeDuMi 
%               - SDPT3 
%       If solver's options must be modified, edit
%       options.SDPsolver.solverOptions to to set solver's options.
%       For example, to set the parameter 'beta' to 0.5 in SeDuMi, we have 
%       to write :
%       >> options.SDPsolver.solver = 'sedumi';
%       >> options.SDPsolver.solverOptions.beta = 0.5;
%       And the call
%       >> [X,Y,feas,stopFlag]=solve_semi_definite_program(A,b,c,K,reltol,options)
%       gives to SeDuMi the desired parameters. To see all available
%       parameters for solvers, type HELP "Solver_Used".
%
%     # Case where solver is a function handle (advanced users) :
%
%   	Else, if options.SDPsolver.solver is a function handle, then
%   	[ ... ] = options.SDPsolver.solver(A,b,c,K,reltol,options.SDPsolver.solverOptions) 
%   	is used to solve the SDP problem. The field SDPsolver.solverOptions
%   	is the option given to the interface. For example, if the custom
%   	interface requires a structure with the field maxIter, then we have
%   	to write :
%       >> options.SDPsolver.solver = @customInterface;
%       >> options.SDPsolver.solverOptions.maxIter = ...;
%       And the call
%       >> [ ... ] = solve_semi_definite_program(A,b,c,K,tol,reltol,options)
%       gives to the interface the desired parameters.
%
%       WARNING : In the case where a custom interface is used, numerical
%       errors must be treated inside!
%
%       WARNING : In the case where a custom interface is used, be careful 
%       to give the same output than this interfaces :
%       X : vectorized primal variables.
%       Y : vectorized dual variables.
%       feas = 0 if infeasible, 1 if feasible.
%       stopFlag
%           = 0 if all is ok
%           = 1 if the problem is feasible/infeasible, but primal
%               and/or dual solutions are not accurate enough with
%               respect to the given tolerance.
%           = 2 if complete failure (the solver can not to decide if
%               the problem is feasible or not).
%       objectiveVal : the value of the objective function.
%
% Sedumi can be downloaded here :
% http://perso.uclouvain.be/raphael.jungers/sites/default/files/sedumi.zip
%

stopFlag = [];
feas = false; % by default, the problem is infeasible.

solver_detected = 0;

if(nargin <6)
    options = jsrsettings;
    defaultInterface = true;
else
    if ischar(options.SDPsolver.solver)
        defaultInterface = true;
        options.SDPsolver.solver = lower(options.SDPsolver.solver);
    elseif isa(options.SDPsolver.solver, 'function_handle')
        defaultInterface = false;
        solver_detected = 1;
    else
        error('Invalid solver.');
    end
end

if(defaultInterface)
    
    
    %%%%%%%%%%%%%%%%%
    % SeDuMi Solver %
    %%%%%%%%%%%%%%%%%
    if strcmpi(options.SDPsolver.solver,'sedumi')
        
        % Check existence of SeDuMi
        s = which('sedumi.m');
        if isempty(s)
            error(['Solver sedumi could not be found. Please install SeDuMi.'...
                'See http://sedumi.ie.lehigh.edu/']);
        end
        
        % Check version of SeDuMi
        if not(check_version_sedumi)
            warning(['solve_semi_definite_program has detected that an outdated version of SeDuMi is used. ',...
                'It could cause the method to crash. We strongly advise to download the lastest ',...
                'version of SeDuMi on', ...
                'http://perso.uclouvain.be/raphael.jungers/sites/default/files/sedumi.zip.'])
        end
        
        solver_detected = 1;
        
        % Check options of SeDuMi
        if not(isfield(options.SDPsolver,'solverOptions'))
            options.SDPsolver.solverOptions = [];
        end
        
        pars = options.SDPsolver.solverOptions;
        
        if not(isfield(pars,'fid'))
            pars.fid = 0; % Silent.
        end
        
        pars.eps = reltol;
        
        % Disable rank deficient warning.
        s = warning('off','MATLAB:rankDeficientMatrix');
        % Solve SDP
        [X,Y,info] = sedumi(A,b,c,K,pars);
        %Return warning state
        warning(s);
        stopFlag = 0;
        objectiveVal = [];
        if(info.pinf == 1 && ~isfield(info,'numerr'))
            % case where sedumi detects the infeasibility in the
            % pre-process
            stopFlag = 0;
            feas = false;
        else
        
            % Postprocess
            if( info.numerr >= 2 ) % Total failure
                warning('SeDuMi : Numerical errors occurred');
                stopFlag = 2;
            else
                if(info.numerr == 1) % X and Y incacurate
                    stopFlag = 1;
                    objectiveVal = mean([c(:)'*X, b(:)'*Y]); 
                end

                % Check feasability.
                if(info.pinf == 0 && info.dinf == 0 && info.feasratio > 0) % Feasible
                    feas = true;
                    objectiveVal = mean([c(:)'*X, b(:)'*Y]); 
                else
                    feas = false; % Infeasible
                end
            end
        end
    end
    
    
    
    %%%%%%%%%%%%%%%%
    % SDPT3 Solver %
    %%%%%%%%%%%%%%%%
    
    if strcmpi(options.SDPsolver.solver,'sdpt3')
        solver_detected = 1;
        s = which('sdpt3.m');
        if isempty(s)
            error(['Solver SDPT3 could not be found. Please install SDPT3 properly.'...
                'See http://www.math.nus.edu.sg/~mattohkc/sdpt3.html']);
        end
        
        if not(isfield(options.SDPsolver.solverOptions,'printLevel'))
            options.SDPsolver.solverOptions.printlevel = 0; % silent
        end
        
        % tolerance
        options.SDPsolver.solverOptions.inftol = reltol;
        options.SDPsolver.solverOptions.rmdepconstr = 1;
        options.SDPSolver.solverOptions.warning = 0;
        
        
        [blk,At,C,bb] = read_sedumi(A,b,c,K); % Sedumi to SDPT3
        [objectiveVal,X,Y,~,info] = sdpt3(blk,At,C,bb, options.SDPsolver.solverOptions);
        objectiveVal = mean(objectiveVal);
        % resize output (cell to vec).
        primalVar = zeros(size(A,1),1);
        idx = 1;
        for i=1:length(X)
            submatrixX = X{i};
            primalVar(idx : (idx+numel(submatrixX)-1)) =  submatrixX(:) ;
            idx = idx + numel(X{i});
        end
        
        dualVar = Y(:) ;
        
        if(info.termcode > 0) % Primal or dual infeasible
            feas = 0;
        else
            feas = 1;
        end
        X = primalVar;
        Y = dualVar;
    end
    
    
end

if not(defaultInterface)
    [X,Y,feas,stopFlag,objectiveVal] = options.SDPsolver.solver(A,b,c,K,tol,options.SDPsolver.solverOptions);
end

if(solver_detected == 0)
    error('Invalid solver.')
end

    function out = check_version_sedumi()
        % If version is Ok, then out = 1. 
        olderCorrectVersion = 20130724;
        VersionSedumi = load('Version.txt');
        VersionSedumi = VersionSedumi(2);
        out = VersionSedumi >= olderCorrectVersion;
        
    end

end
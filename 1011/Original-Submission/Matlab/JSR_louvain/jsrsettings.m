function opts = jsrsettings(varargin)
% 
% JSRSETTINGS Constructs option structure for the algorithms
%
% OPTIONS = JSRSETTINGS with no input arguments returns
%   setting structure with default values
%
% OPTIONS = JSRSETTINGS('NAME1',VALUE1,'NAME2',VALUE2,...) creates a
%   solution options structure OPTIONS in which the named properties have
%   the specified values.  Any unspecified properties have default values.
%   It is sufficient to type only the leading characters that uniquely
%   identify the property.  Case is ignored for property names.
%
% OPTIONS = JSRSETTINGS(OLDOPTS,'NAME1',VALUE1,...) alters an existing options
%   structure OLDOPTS.
%
%
% JSRSETTINGS PROPERTIES :
%
% GENERAL
% 
%   saveinIt         - Save normal output in VALUE.mat at each iteration, default : (0). 
%                      Value must be a string, see example below
%   saveEnd          - Save outputs in VALUE.mat (only at termination)
%
%   verbose          - 0 the method is silent, 1 normal, >1 loud, (1)
%
%   maxTime          - Approximate max time in seconds one allows the
%                      algorithm to run. At the end of each iteration the  
%                      algorithm checks the elapsed time, if it has run for
%                      more than maxTime or maxTime - (time of last iteration)
%                      it stops iterating.
%                      Termination due to reaching of maxTime is signaled
%                      with a msg and in info.status. (Inf)
%       
%   logfile          - Value 0 writes no log-file (default value).
%                      Value 1 prints all possible messages from algorithms
%                      to a text file named "log_<nameOfAlgorithm>" in current
%                      directory, overwriting the file at each call. 
%                      Setting a string value allows to specify a name for the file
%                      in which to write the log.
%                      One can also specify a valid FID created by fopen to write at 
%                      the end of that file (in that case, functions will
%                      not close the file).
%
% BALANCEDCOMPLEXPOLYTOPE
%
%   options.bcp, see help jsr_norm_balancedComplexPolytope
%
% BALANCEDREALPOLYTOPE
% 
%   options.brp, see help jsr_norm_balancedRealPolytope
%
% BRUTEFORCE
%
%   options.bruteforce, see help jsr_prod_bruteForce
%   NOTE : options.saveinIt & options.maxTime not available with bruteForce
%
% CONITOPE
%   
%   options.conitope, see help jsr_norm_conitope
%
% ELLIPSOID
%
%   options.ellips, see help jsr_conic_ellipsoid
%
% GRIPENBERG
%   
%   options.grip, see help jsr_prod_Gripenberg
%
% ITMETH
%
%   options.itMeth, see help itComp
%
% JSR
%   
%   options.jsr, see help jsr
%
% LINEAR
% 
%   options.linear, see help jsr_conic_linear
%
% LINEARRELAXATION2D
%
%   options.linrel, see help jsr_norm_linearRelaxation2D
%
% LOWERBRUTEFORCE
% 
%   options.lowbrut, see help jsr_prod_lowerBruteForce
% 
% MAXRELAXATION
% 
%   options.maxrel, see help jsr_norm_maxRelaxation
%
% MAXRELAXATION2D
%
%   options.maxrel2D, see help jsr_norm_maxRelaxation2D
%
% PRUNINGALGORITHM
%
%   options.pruning, see help jsr_prod_pruningAlgorithm
%                            
% SEMIDEFINITE
% 
%   options.semidef, see help jsr_lift_semidefinite
%
% SOS
%
%   options.sos, see help jsr_opti_sos
% 
% pathcomplete
%
%   options.pathcomplete, see help jsr_pathcomplete
%   
%
% ------------------------------------------------------
% Example :
%          
%   >> opts = jsrsettings('saveinIt','conitope_test1','conitope.maxiter',42) 
%   >> [bounds, prodOpt, prodV, info] = jsr_norm_conitope(M,opts)
%              
%  Executes algorithm with maxiter = 42 and saves current output values 
%  in the file conitope_test1.mat                              
% ------------------------------------------------------      

% (Ripped from sdpsettings.m by Johan L�fberg)

% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    help jsrSettings
    return;
end


Names = {
    % General
    'saveinIt'
    'saveEnd'
    'verbose'
    'logfile'
    'maxTime'
    
    % balancedComplexPolytope
    'bcp.maxiter'
    'bcp.plotBounds'
    'bcp.plotpopHist'
    'bcp.plotTime'
    
    % balancedRealPolytope
    'brp.maxiter'
    'brp.plotBounds'
    'brp.plotpopHist'
    'brp.plotTime'
    
    % bruteForce
    'bruteforce.maxdepth'
    'bruteforce.normfun'
    'bruteforce.plotBounds'
    
    % conitope
    'conitope.initProd'
    'conitope.maxiter'
    'conitope.tol'
    'conitope.per'
    'conitope.plotBounds'
    'conitope.plotpopHist'
    'conitope.reinitBall'
    'conitope.plotTime'
    
    % ellipsoid
    'ellips.initBounds'
    'ellips.maxiter'
    'ellips.tol'
    
    % grip
    'grip.delta'
    'grip.normfun'
    'grip.maxEval'
    'grip.plotBounds'
    'grip.plotTime'
    'grip.plotpopHist'
    
    %itMeth
    'itMeth.lift'
    'itMeth.length'
    'itMeth.meth'
    'itMeth.tol'
    'itMeth.plotBounds'
    'itMeth.plotTime'
    
    % jsr
    'jsr.time'
    'jsr.tol'
    'jsr.triang'
    
    % linear
    'linear.initBounds'
    'linear.maxiter'
    'linear.tol'
        
    % linearRelaxation2D
    'linrel.step'
    'linrel.lambda'
    'linrel.maxiter'
    'linrel.tol'
    'linrel.plotBounds'
    'linrel.plotEllips'
    'linrel.plotTime'
    
    % lowerBruteForce
    'lowbrut.depths'
    'lowbrut.ndisp'
    'lowbrut.plotBounds'
    'lowbrut.plotTime'
    
    % maxRelaxation
    'maxrel.N'
    'maxrel.avgfun'
    'maxrel.maxiter'
    'maxrel.tol'
    'maxrel.plotBounds'
    'maxrel.plotTime'
    
    % maxRelaxation2D
    'maxrel2D.step'
    'maxrel2D.avgfun'
    'maxrel2D.tol'
    'maxrel2D.plotBounds'
    'maxrel2D.plotEllips'
    'maxrel2D.plotTime'
    
    % pruningAlgorithm
    'pruning.maxEval'
    'pruning.delta'
    'pruning.normfun'
    'pruning.domfun'
    'pruning.plotBounds'
    'pruning.plotpopHist'
    'pruning.plotTime'
    
    % semidef
    'semidef.maxDepth'
    'semidef.plotBounds'
    'semidef.plotTime'
    
    % sos
    'sos.initBounds'
    'sos.deg'
    'sos.maxiter'
    'sos.tol'
    
    % get_graph
    'graphOptions.type'
    
    % DeBruijn
    'debruijn.dimension'
    
    % standart_opti_form
    % TODO
    
    % SDPsolver
    'SDPsolver.solver'
    'SDPsolver.solverOptions'
    
    % pathcomplete
    'pathcomplete.LbBisec'
    'pathcomplete.UbBisec'
    'pathcomplete.maxiter'
    'pathcomplete.reltol'
    'pathcomplete.abstol'
    'pathcomplete.testUb'
    'pathcomplete.testLb'
    'pathcomplete.loadIt'
    'pathcomplete.loadItFile'
};


obsoletenames ={ % Use when options have become obsolete
};

[m,n] = size(Names);
names = lower(Names);

if (nargin>0) & isstruct(varargin{1})
    opts = varargin{1};
    paramstart = 2;
else
    paramstart = 1;
    
    % General 
    opts.saveinIt = 0;
    opts.saveEnd = 0;
    opts.verbose = 1;
    opts.logfile = 0;
    opts.maxTime = Inf;
    
    % balancedComplexPolytope
    opts.bcp.maxiter = 500;
    opts.bcp.plotBounds = 0;
    opts.bcp.plotpopHist = 0;
    opts.bcp.plotTime = 0;
    
    % balancedRealPolytope
    opts.brp.maxiter = 1000;
    opts.brp.plotBounds = 0;
    opts.brp.plotpopHist = 0;
    opts.brp.plotTime = 0;
    
    % bruteForce
    opts.bruteforce.maxdepth = 4;
    opts.bruteforce.normfun = @norm;
    opts.bruteforce.plotBounds = 0;
    
    % conitope
    opts.conitope.initProd = [];
    opts.conitope.maxiter = 60;
    opts.conitope.tol = 1e-8;
    opts.conitope.per = 1;
    opts.conitope.plotBounds = 0;
    opts.conitope.plotpopHist = 0;
    opts.conitope.reinitBall = 1;
    opts.conitope.plotTime = 0;
    
    % ellipsoid
    opts.ellips.initBounds = [];
    opts.ellips.maxiter = 1000;
    opts.ellips.tol =1e-8;
    
    % grip
    opts.grip.delta = 1e-2;
    opts.grip.normfun = @norm;
    opts.grip.maxEval = 1000;
    opts.grip.plotBounds =0;
    opts.grip.plotTime = 0;
    opts.grip.plotpopHist = 0;
    
    % itComp
    opts.itMeth.lift = 1;
    opts.itMeth.length = 4;
    opts.itMeth.meth = {'ellips'};
    opts.itMeth.tol = 1e-8;
    opts.itMeth.plotBounds = 0;
    opts.itMeth.plotTime = 0;
    
    % jsr
    opts.jsr.time = 60*2;
    opts.jsr.tol  = 1e-6;
    opts.jsr.triang = 1;
    
    % linear
    opts.linear.initBounds = [];
    opts.linear.maxiter = 1000;
    opts.linear.tol =1e-8;
    
    % linearRelaxation2D
    opts.linrel.step = 2*pi/10000;
    opts.linrel.lambda = 0.42;
    opts.linrel.maxiter = 1000;
    opts.linrel.tol = 1e-10;
    opts.linrel.plotBounds = 0 ;
    opts.linrel.plotEllips = 0;
    opts.linrel.plotTime = 0;
    
    % lowerBruteForce
    opts.lowbrut.depths = 1:2:5;
    opts.lowbrut.ndisp = 0;
    opts.lowbrut.plotBounds = 0;
    opts.lowbrut.plotTime = 0;
    
    % maxRelaxation
    opts.maxrel.N = 1e4;
    opts.maxrel.avgfun = 'a';
    opts.maxrel.maxiter = 500;
    opts.maxrel.tol = 1e-8;
    opts.maxrel.plotBounds = 0;
    opts.maxrel.plotTime = 0;
    
    % maxRelaxation2D
    opts.maxrel2D.step = 2*pi/(10000);
    opts.maxrel2D.avgfun = 'a';
    opts.maxrel2D.maxiter = 1000;
    opts.maxrel2D.tol = 1e-10;
    opts.maxrel2D.plotBounds = 0;
    opts.maxrel2D.plotEllips = 0;
    opts.maxrel2D.plotTime = 0;
    
    % pruning
    opts.pruning.maxEval = 1000;
    opts.pruning.delta = 1e-8;
    opts.pruning.normfun = @norm;
    opts.pruning.domfun = 'defaultdomfun';
    opts.pruning.plotBounds = 0;
    opts.pruning.plotpopHist = 0;
    opts.pruning.plotTime = 0;
    
    % semidef
    opts.semidef.maxDepth = 6;
    opts.semidef.plotBounds = 0;
    opts.semidef.plotTime = 0;
    
    % sos
    opts.sos.initBounds = [];
    opts.sos.deg = 1;
    opts.sos.maxiter = 1000;
    opts.sos.tol = 1e-10;
    
    % get_graph
    opts.graphOptions.type = 'DeBruijn';
    
    % DeBruijn
    opts.debruijn.dimension = 1 ;
    
    % standart_opti_form
    % Not available yet.
    
    
    % SDPsolver
    opts.SDPsolver.solver = 'SeDuMi';
    opts.SDPsolver.solverOptions = [];
    
    % pathcomplete
    opts.pathcomplete.maxiter = inf;
    opts.pathcomplete.LbBisec = [];
    opts.pathcomplete.UbBisec = [];
    opts.pathcomplete.reltol = 1e-6;
    opts.pathcomplete.abstol = inf;
    opts.pathcomplete.solver = 'SeDuMi';
    opts.pathcomplete.solverOptions = [];
    opts.pathcomplete.testUb = 1;
    opts.pathcomplete.testLb = 1;
    opts.pathcomplete.loadIt = 0;
    opts.pathcomplete.loadItFile = 'saveinIt';
    
end

i = paramstart;
% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
    error('Arguments must occur in name-value pairs.');
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
    arg = varargin{i};
    
    if ~expectval
        if ~ischar(arg)
            error(sprintf('Expected argument %d to be a string property name.', i));
        end
        
        lowArg = lower(arg);
        
        j_old = strmatch(lowArg,obsoletenames);
        if ~isempty(j_old)
            % For compability... No need yet
        end
        
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
            error(sprintf('Unrecognized property name ''%s''.', arg));
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                msg = sprintf('Ambiguous property name ''%s'' ', arg);
                msg = [msg '(' deblank(Names{j(1)})];
                for k = j(2:length(j))'
                    msg = [msg ', ' deblank(Names{k})];
                end
                msg = sprintf('%s).', msg);
                error(msg);
            end
        end
        expectval = 1;                      % we expect a value next
    else
        eval(['opts.' Names{j} '= arg;']);
        expectval = 0;
    end
    i = i + 1;
end

if expectval
    error(sprintf('Expected value for property ''%s''.', arg));
end

end




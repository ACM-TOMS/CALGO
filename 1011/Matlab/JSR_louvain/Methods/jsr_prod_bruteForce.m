function [bounds, P, info] = jsr_prod_bruteForce(M, varargin)

% JSR_PROD_BRUTEFORCE Approximates the jsr using brute force.
%
%    [BOUNDS, P, INFO] = JSR_PROD_BRUTEFORCE(M)
%      returns lower and upper bounds on the jsr of M using brute force by
%      considering products of length up to 4. The set M must be a cell
%      array of matrices.
%
%    [BOUNDS, P, INFO] = JSR_PROD_BRUTEFORCE(M, MAXDEPTH, NORMFUN)
%      returns lower and upper bounds on the jsr of M using brute force by
%      considering products of length up to MAXDEPTH.
%      The parameter NORMFUN, if specified, should be a function handle.       
%
%    [BOUNDS, P, INFO] = JSR_PROD_BRUTEFORCE(M, OPTS)
%      Does the same but with values of parameters defined in the structure OPTS.
%      See below and help JSRSETTINGS for available parameters.
%
%      BOUNDS contains the lower and upper bounds on the JSR
%
%      P is a cell containing the indexes of two products attaining the
%        respective bounds, (to use with buildProduct)
%
%      INFO is a structure containing various data about the iterations :
%         info.elapsedtime
%         info.allLb          - lower bounds at each depth
%         info.allUb          - upper bounds at each depth
%         info.opts           - the structure of options used
%
%  The field opts.bruteforce (generated by jsrsettings) can be used to
%  tune the method :
%
%      bruteforce.maxdepth       - depth to which all products are
%                                  computed, (4)
%      bruteforce.normfun        - function handle to a matrix norm,
%                                  (@norm)
%      bruteforce.plotBounds     - if 1 plots the evolution of the bounds 
%                                  w.r.t. the depth, (0)
%
% NOTE : options opts.saveinIt & opts.maxTime do not work here because of 
%        recursive implementation.
%
% See also JSRSETTINGS
% 
% REFERENCES
%    R.Jungers, 
%      "The Joint Spectral Radius: Theory and Applications" 
%      Vol. 385 in Lecture Notes in Control and Information
%      Sciences, Springer-Verlag. Berlin Heidelberg, June 2009


if(nargin>1)
    
    if (length(varargin)==2 && isnumeric(varargin{1}) && isa(varargin{2},'function_handle'))
        opts = jsrsettings('bruteforce.maxdepth',varargin{1},'bruteforce.normfun',varargin{2});
    elseif (length(varargin)==1 && isnumeric(varargin{1}))
        opts = jsrsettings('bruteforce.maxdepth',varargin{1});
    elseif (length(varargin)==1 && isstruct(varargin{1}))
        opts = varargin{1};
    else
        warning('Invalid optional arguments, using default parameters !');
        opts = jsrsettings;
    end

else
    opts = jsrsettings;
end


% logfile opening
close = 1;
if (ischar(opts.logfile) )    
    logFile = fopen(opts.logfile,'wt');
    if (logFile == -1)
        warning(sprintf('Could not open file %s',opts.logfile));
    end
elseif isnumeric(opts.logfile)
    if (opts.logfile==0)
        logFile = -1;
    elseif (opts.logfile ==1)
        logFile = fopen('log_bruteForce','wt');
        if (logFile == -1)
            warning('Could not open logfile')
        end
    else
        logFile = opts.logfile;
        close =0;
    end
else
    logFile = fopen('log_bruteForce','wt');
    if (logFile == -1)
        warning('Could not open logfile')
    end
end

if (logFile~=-1)
    fprintf(logFile,[datestr(now) '\n\n']);
end

msg(logFile,opts.verbose>1,'\n \n******** Starting jsr_prod_bruteForce ******** \n \n')

% Initialization
starttime = cputime;
n = length(M{1}); % dimension of the matrices
maxdepth = opts.bruteforce.maxdepth;
normfun = opts.bruteforce.normfun;

% Recursively scan through all products of length <= maxdepth
[lowpart, uppart, lowerP, upperP] = bfRecursive(...
        M, normfun, eye(n), zeros(1, maxdepth), maxdepth, 0,logFile,opts);
[lower, idx1] = max(lowpart.^(1./(1:maxdepth)));
[upper, idx2] = min(uppart.^(1./(1:maxdepth)));
P = {lowerP{idx1}, upperP{idx2}};
% Reset product index convention
P{1} = P{1}(length(P{1}):-1:1);
P{2} = P{2}(length(P{2}):-1:1);

% Post-processing
msg(logFile,opts.verbose>0,'> Bounds on the jsr: [%.15g, %.15g]', lower, upper);
bounds = [lower, upper];
elapsedtime = cputime - starttime;

msg(logFile,opts.verbose>1,'\n End of algorithm after %5.2f s',elapsedtime)

if (logFile~=-1 && close)
    fclose(logFile);
end

info.elapsedtime = elapsedtime;
info.allLb = lowpart.^(1./(1:maxdepth));
info.allUb = uppart.^(1./(1:maxdepth));
info.opts = opts;

if(ischar(opts.saveEnd))
    save([opts.saveEnd '.mat'],'bounds','P','info')
end

if(opts.bruteforce.plotBounds)
    figure
    plot(1:maxdepth,info.allLb,'-+g',1:maxdepth,info.allUb,'-*r')
    title('bruteForce : Evolution of the bounds')
    xlabel('Depth')
    legend('Lower bound','Upper bound')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lower, upper, lowerP, upperP] = ...
        bfRecursive(M, normfun, current, word, maxdepth, curdepth,logFile,options)

lower(maxdepth) = 0;
upper(maxdepth) = 0;

% Produit de longueur non nulle
if (curdepth > 0),
    optseigs.disp = 0;
    lower(curdepth) = abs(eigs(current, [], 1, 'LM', optseigs));
    upper(curdepth) = feval(normfun, current);
    lowerP{curdepth} = word(1:curdepth);
    upperP{curdepth} = word(1:curdepth);
end

% Limite de longueur des produits non atteinte
if (curdepth < maxdepth),
    curdepth = curdepth + 1;
    for i = 1:length(M),
        word(curdepth) = i; %#ok<AGROW>
        [sublower, subupper, sublowerP, subupperP] = bfRecursive(...
                M, normfun, current*M{i}, word, maxdepth, curdepth,logFile,options);
        
        [lower, lowerM] = max([lower ; sublower], [], 1);
        [upper, upperM] = max([upper ; subupper], [], 1);
        lowerM = (lowerM == 2);
        upperM = (upperM == 2);
        lowerP(lowerM) = sublowerP(lowerM);
        upperP(upperM) = subupperP(upperM);
        
    end
end
end
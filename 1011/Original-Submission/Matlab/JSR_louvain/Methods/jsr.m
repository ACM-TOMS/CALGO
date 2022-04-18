function [bounds, info] = jsr(varargin)
%
%  [BOUNDS, INFO] = JSR(A1,A2,...,Am) or [BOUNDS, INFO] = JSR(M)
%       Launches different methods to find bounds on the the joint spectral      
%       radius of the square matrices A1,A2,...,Am or contained in the cell
%       array M. Is limited in time to 2*60 seconds.
%  
%  [BOUNDS, INFO] = JSR(M,'time',T) or [BOUNDS, INFO] = jsr(A1,...,Am,'time',T)
%       Uses the time limit T in seconds instead of default 2*60.
%  
%  [BOUNDS, INFO] = JSR(M,OPTS) or [BOUNDS, INFO] = jsr(A1,...,Am,OPTS)
%       Uses the paramaters values specified in the structure OPTS
%       generated with jsrsettings. See below and jsrsettings for available
%       options and parameters.
%  
%  BOUNDS contains the best lower and upper bounds found on the jsr 
%
%  INFO is a structure containing various data about the algorithms and
%       their outputs.
%       If the matrices could be block-triangularized then the sets of blocks
%       will be returned in the first field, the second will contain the 
%       quick bounds computed for each one of them in order to quickly
%       eliminate irrelevant blocks and the other fields will contain the
%       INFO of the recursive JSR calls on each relevant block-sets along 
%       with the bounds found.                         
%       If no block-triangularization was performed the normal output of 
%       each method called is returned in a field named after the concerned  
%       method.
%  INFO.opts contains the OPTS structure used with JSR.
%
%  The field opts.jsr (generated by jsrsettings) can be used 
%  to tune the method:
%
%      jsr.time               - time limit in seconds to be divided between 
%                               the methods called, (2*60)    
%      jsr.tol                - tolerance asked to called methods, (1e-6)
%                               (see help jsr_conic_linear, jsr_conic_ellipsoid,
%                               jsr_pprod_Gripenberg and jsr_norm_conitope  
%                               for precise meaning of this tolerance)
%      jsr.triang             - if 1 JSR tries first to block-triangularise 
%                               the matrices, (1)
%
%  Structure of the algorithm:
%   0) Checks for negative or complex entries
%     - If matrices have real nonnegative entries:
%        1.1) Tries to block-triangularise with permTriangul, if did not
%        work, try with jointTriangul. 
%          1.2) If blocks could be found, tries a quick elimination and then
%          launches JSR on the remaining sets of blocks.
%          1.2') If could not triangularise, launches pruningAlgorithm, then
%          conic_linear, then conic_ellips and, finally, jsr_norm_conitope 
%          with as initial prodOpt the one attaining lower bound of pruning.
%     - If matrices have positive & negative entries or complex entries:
%        2.1) Tries to block-triangularise with permTriangul, if did not
%        work, try with jointTriangul.
%          2.2) If blocks could be found, tries a quick elimination and then
%          launches JSR on the remaining sets of blocks.
%          2.2') If could not triangularise, launches Gripenberg, then launches 
%          method ellips, finally launches conitope with as initial prodOpt  
%          the one from Gripenberg.
%
%  See also JSRSETTINGS, jsr_prod_pruningAlgorithm, jsr_conic_linear,
%  jsr_conic_ellipsoid, jsr_prod_Gripenberg, jsr_norm_conitope,
%  jsr_norm_balancedComplexPolytope
%
% Requires SeDuMi ( http://sedumi.ie.lehigh.edu/ )
%
% Please report any bug, comment or 
% suggestion to jsr.louvain@gmail.com

if isnumeric(varargin{1})
    M = [];
    iarg = 1;
    while (isnumeric(varargin{iarg}))
        M = [M {varargin{iarg}}];
        iarg = iarg+1;
        if (iarg>nargin)
            break;
        end
    end
else
    M = varargin{1};
    iarg = 2;
end

opts = jsrsettings;

while (nargin>= iarg)
    if isstruct(varargin{iarg})
        opts = varargin{iarg};
        iarg=iarg+1;
    elseif ischar(varargin{iarg})
            opts = jsrsettings(['jsr.' varargin{iarg}],varargin{iarg+1});
            iarg = iarg+2;
    end
end

% logfile opening
close = 1;
if (ischar(opts.logfile) )    
    logFile = fopen(opts.logfile,'wt');
    if (logFile == -1)
        warning(sprintf('Could not open file %s',opts.logfile));
    end
elseif isnumeric(opts.logfile)
    if (opts.logfile ==0 )
        logFile = -1;
    elseif (opts.logfile==1)
        logFile = fopen('log_JSR','wt');
        if (logFile == -1)
            warning('Could not open logfile')
        end
    else
        logFile = opts.logfile;
        close = 0;
    end
else
    logFile = fopen('log_JSR','wt');
    if (logFile == -1)
        warning('Could not open logfile')
    end
end

if (logFile~=-1)
    fprintf(logFile,[datestr(now) '\n\n']);
end

msg(logFile,opts.verbose>1,'\n \n******** Starting jsr ******** \n \n')
starttime = cputime;

m = length(M);

if (m==1)
    msg(logFile,opts.verbose>0,'\n>>JSR: Only one matrix in set. Joint spectral radius = spectral radius. \n')
    bounds = max(rho(M))*ones(1,2);
    
    return;
end

di=size(M{1});
if (di==[1 1])
    msg(logFile,opts.verbose>0,'\n>>JSR: dimension 1. Joint spectral radius = max absolute value. \n')
    bounds=0;
for imat=1:m
    if (abs(M{imat})>bounds)
        bounds=M{imat};
    end
end
    bounds=[bounds bounds];
    info=0;
    return;
end

% Test nonnegativity & complex entries
nonneg = 1;

for imat=1:m
    if (sum(sum(M{imat}<0))>0 || ~isreal(M{imat}))
        nonneg =0;
    end
end

if (nonneg)
    % Nonnegative case
    msg(logFile,opts.verbose>0,'Matrices have nonnegative entries. \n')
    
    triang = 0;
    triangPerm = 0;
    if (opts.jsr.triang)
        % Try to block-triangularize by permutation
        [triangPerm, blocks, perm] = permTriangul(M);
        
        if ~triangPerm
            % Try to block triangularize
            [triang, blocks, B] = jointTriangul(M);
        end
    end
        
    if (triang || triangPerm)
        if triang
            msg(logFile,1,'Your matrices could be jointly block-triangularized with jointTriangul, (in %d blocks)',length(blocks));
            msg(logFile,1,'\nWARNING: \n  jointTriangul is a heuristic and could have unexpected behavior, e.g., on badly conditioned matrices.');
            msg(logFile,1,'   Please see help jointTriangul and verify the joint-triangularization.');
            msg(logFile,1,'   Consider counter-checking results with option jsr.triang to 0.\n');
            
            fprintf('Press any key to proceed.\n')
            pause
        else
            msg(logFile,1,'Matrices could be jointly block-triangularized by permutation, (in %d blocks)',length(blocks));
        end
        
        % Try to quickly eliminate blocks (give around 5 seconds max)
        nblo = length(blocks);
        lbBlo = zeros(1,nblo);
        ubBlo = zeros(1,nblo);
        
        for iblo = 1:nblo
            Msub = blocks{iblo};
            if (size(Msub{1},1)==1)
                lbBlo(iblo) = max(cell2mat(Msub));
                ubBlo(iblo) = lbBlo(iblo);
            else
                prunTime  = min(opts.jsr.time/(10*nblo),5/(nblo));
                optsQuickPrun = jsrsettings('verbose',0,'logfile',0,'maxTime',prunTime);
                bQuickPrun = jsr_prod_pruningAlgorithm(Msub,optsQuickPrun);
                lbBlo(iblo) = bQuickPrun(1);
                ubBlo(iblo) = bQuickPrun(2);
            end
        end
        
        if nargout>1
            info.allBlocks = blocks;
            info.quickBoundsBlocks = [lbBlo;ubBlo];
            
            if triangPerm
                info.triangPermutation = perm;
            else
                info.triangMat = B;
            end
        end
            
        elim = zeros(1,m);
            
        for iblo=1:nblo
            elim(ubBlo<lbBlo(iblo))=1;
        end             

        msg(logFile,opts.verbose>0,'Could quickly eliminate %d block set(s) out of %d',sum(elim),nblo);
        
        blocks(elim==1) = [];
        lbBlo(elim==1) = [];
        ubBlo(elim==1) = [];
        nblo = length(blocks);
              
        for iblo = 1:nblo
            msg(logFile,opts.verbose>0,'-----------------------------')
            msg(logFile,opts.verbose>0,'Launching jsr on block set %d',iblo);
            msg(logFile,opts.verbose>0,'-----------------------------')
            optssubJSR = jsrsettings(opts,'logfile',logFile,'jsr.time',opts.jsr.time/(nblo));
            [boundsi, infoi] = jsr(blocks{iblo},optssubJSR);
            lbBlo(iblo) = max(boundsi(1),lbBlo(iblo));
            ubBlo(iblo) = min(boundsi(2),ubBlo(iblo));
            infoi.bounds = [lbBlo(iblo), ubBlo(iblo)];
            eval(sprintf('info.block%d = infoi;',iblo));
        end
        
        % Second elimination
        if (nblo>1)
            elim = zeros(1,nblo);
            for iblo=1:nblo
                elim(ubBlo<lbBlo(iblo)) = 1;
            end
            msg(logFile,opts.verbose>0,'Could eliminate %d block set(s) out of %d',sum(elim),nblo);
            blocks(elim==1) = [];
            lbBlo(elim==1) = [];
            ubBlo(elim==1) = [];
            nblo = length(blocks);
        end
        
        % Cannot choose between blocs 
        if (nblo>1)
            msg(logFile,1,'\nJSR cannot further eliminate blocks.');
            msg(logFile,1,'First output will be double:')
            msg(logFile,1,'out.bounds    is a vector with on each row the bounds for a block set')
            msg(logFile,1,'out.blocks    is a cell with the corresponding block sets\n')
            msg(logFile,1,'The joint spectral radius of M is the max of the joint spectral radii of each')
            msg(logFile,1,'block sets.')
            bounds.bounds = [lbBlo' ubBlo'];
            bounds.blocks = blocks;
        else
            bounds = [lbBlo', ubBlo'];
        end
        
    else
        
        % Limit pruning on time and tol, not on number of evaluations
        optsprun   = jsrsettings('logfile',logFile,'verbose',opts.verbose-1,...
                                 'pruning.maxEval',Inf,'pruning.delta',opts.jsr.tol,'maxTime',opts.jsr.time/4);
        
        msg(logFile,opts.verbose>0,'\n>> JSR: Launching pruningAlgorithm\n')
        [boundsprun, Pprun, infoPrun] = jsr_prod_pruningAlgorithm(M,optsprun);
        
        bounds(1) = boundsprun(1);
        bounds(2) = boundsprun(2);
        
        optslinear = jsrsettings('logfile',logFile,'verbose',opts.verbose-1,...
                                 'linear.initBounds',bounds,'maxTime',(opts.jsr.time-(cputime-starttime))/3,'linear.tol',opts.jsr.tol);
        msg(logFile,opts.verbose>0,'\n>> JSR: Launching conic_linear\n')
        [boundslin, betterlin, xlin, infolin]  = jsr_conic_linear(M,optslinear);
        
        bounds(1) = max([bounds(1) boundslin(1)]);
        bounds(2) = min([bounds(2), boundslin(2)]);
        
        optsellips = jsrsettings('logfile',logFile,'verbose',opts.verbose-1,...
                                 'ellips.initBounds',bounds,'maxTime',(opts.jsr.time-(cputime-starttime)/2),'ellips.tol',opts.jsr.tol);
        msg(logFile,opts.verbose>0,'\n>> JSR: Launching conic_ellipsoid\n')
        [boundsellips, Pellips, infoellips] = jsr_conic_ellipsoid(M,optsellips);
        
        bounds(1) = max([bounds(1) boundsellips(1)]);
        bounds(2) = min([bounds(2), boundsellips(2)]);
        
        optsConitope = jsrsettings('logfile',logFile,'verbose',opts.verbose-1,...
                                 'maxTime',(opts.jsr.time-(cputime-starttime)),'conitope.initProd',Pprun{1},'conitope.tol',opts.jsr.tol);
        msg(logFile,opts.verbose>0,'\n>> JSR: Launching conitope\n')
        [boundsConitope, prodOptconi, prodVconi, infoconi]  = jsr_norm_conitope(M ,optsConitope);

        bounds(1) = max([bounds(1) boundsConitope(1)]);
        bounds(2) = min([bounds(2), boundsConitope(2)]);
        
        
        if nargout>1
            info.pruningAlg.bounds = boundsprun;
            info.pruningAlg.P = Pprun;
            info.pruningAlg.info = infoPrun;
            
            info.linear.bounds = boundslin;
            info.linear.better = betterlin;
            info.linear.x = xlin;
            info.linear.info = infolin;
            
            info.ellips.bounds = boundsellips;
            info.ellips.P = Pellips;
            info.ellips.info = infoellips;
            
            info.conitope.bounds = boundsConitope;
            info.conitope.prodOpt = prodOptconi;
            info.conitope.prodV = prodVconi;
            info.conitope.info = infoconi;
        end
        
        
    end
    
% End of nonnegative part
else
    
    triang = 0;
    triangPerm = 0;
    if (opts.jsr.triang)
        % Try to block-triangularize by permutation
        [triangPerm, blocks, perm] = permTriangul(M);
        
        if ~triangPerm
            % Try to block triangularize
            [triang, blocks, B] = jointTriangul(M);
        end
    end
        
    if (triang || triangPerm)
        if triang
            msg(logFile,1,'Your matrices could be jointly block-triangularized with jointTriangul, (in %d blocks)',length(blocks));
            msg(logFile,1,'\nWARNING: \n  jointTriangul is a heuristic and could have unexpected behavior, e.g., on badly conditioned matrices.');
            msg(logFile,1,'   Please see help jointTriangul and verify the joint-triangularization.');
            msg(logFile,1,'   Consider counter-checking results with option jsr.triang to 0.\n');
            
            fprintf('Press any key to proceed.\n')
            pause
        else
            msg(logFile,1,'Matrices could be jointly block-triangularized by permutation, (in %d blocks)',length(blocks));
        end
        
        % Try to quickly eliminate blocks
        nblo = length(blocks);
        lbBlo = zeros(1,nblo);
        ubBlo = zeros(1,nblo);

        for iblo = 1:nblo
            Msub = blocks{iblo};
            if (size(Msub{1},1)==1)
                lbBlo(iblo) = max(abs(cell2mat(Msub)));
                ubBlo(iblo) = lbBlo(iblo);
            else
                mdbrute  = ceil(log(10)/log(m));
                optsQuickBrute = jsrsettings('verbose',0,'logfile',0,'bruteforce.maxdepth',mdbrute);
                bQuickBrute = jsr_prod_bruteForce(Msub,optsQuickBrute);
                lbBlo(iblo) = bQuickBrute(1);
                ubBlo(iblo) = bQuickBrute(2);
            end
        end

        if nargout>1
            info.allBlocks = blocks;
            info.quickBoundsBlocks = [lbBlo;ubBlo];
            
            if triangPerm
                info.triangPermutation = perm;
            else
                info.triangMat = B;
            end
        end

        elim = zeros(1,m);

        for iblo=1:nblo
            elim(ubBlo<lbBlo(iblo))=1;
        end

        msg(logFile,opts.verbose>0,' Could quickly eliminate %d block set(s) out of %d',sum(elim),nblo);

        blocks(elim==1) = [];
        lbBlo(elim==1) = [];
        ubBlo(elim==1) = [];
        nblo = length(blocks);

        for iblo = 1:nblo
            msg(logFile,opts.verbose>0,'-----------------------------')
            msg(logFile,opts.verbose>0,'Launching jsr on block set %d',iblo);
            msg(logFile,opts.verbose>0,'-----------------------------')
            optssubJSR = jsrsettings(opts,'logfile',logFile,'jsr.time',opts.jsr.time/(nblo));
            [boundsi, infoi] = jsr(blocks{iblo},optssubJSR);
            lbBlo(iblo) = max(boundsi(1),lbBlo(iblo));
            ubBlo(iblo) = min(boundsi(2),ubBlo(iblo));
            infoi.bounds = [lbBlo(iblo), ubBlo(iblo)];
            eval(sprintf('info.block%d = infoi;',iblo));
        end
        
        % Second elimination
        if (nblo>1)
            elim = zeros(1,nblo);
            for iblo=1:nblo
                elim(ubBlo<lbBlo(iblo)) = 1;
            end
            msg(logFile,opts.verbose>0,'Could eliminate %d block set(s) out of %d',sum(elim),nblo);
            blocks(elim==1) = [];
            lbBlo(elim==1) = [];
            ubBlo(elim==1) = [];
            nblo = length(blocks);
        end
        
        % Cannot choose between blocks 
        if (nblo>1)
            msg(logFile,1,'\nJSR cannot further eliminate blocks.');
            msg(logFile,1,'First output will be double:')
            msg(logFile,1,'out.bounds    is a vector with on each row the bounds for a block set')
            msg(logFile,1,'out.blocks    is a cell with the corresponding block sets\n')
            msg(logFile,1,'The joint spectral radius of M is the max of the joint spectral radii of each')
            msg(logFile,1,'block sets.')
            bounds.bounds = [lbBlo' ubBlo'];
            bounds.blocks = blocks;
        else
            bounds = [lbBlo, ubBlo];
        end

    else
        % If could not block-triangularise

        optsgrip   = jsrsettings('logfile',logFile,'verbose',opts.verbose-1,'grip.delta',opts.jsr.tol,'grip.maxEval',3000);
        
        msg(logFile,opts.verbose>0,'\n>> JSR: Launching Gripenberg\n')
        [boundsGrip, prodGrip,allProdsGrip, infoGrip] = jsr_prod_Gripenberg(M,optsgrip);

        bounds(1) = boundsGrip(1);
        bounds(2) = boundsGrip(2);
        
        optsEllips = jsrsettings('logfile',logFile,'verbose',opts.verbose-1,...
                                 'ellips.initBounds',bounds,'maxTime',(opts.jsr.time-(cputime-starttime))/2,'ellips.tol',opts.jsr.tol);

        msg(logFile,opts.verbose>0,'\n>> JSR: Launching conic_ellipsoid\n')
        [boundsellips, Pellips, infoellips] = jsr_conic_ellipsoid(M,optsEllips);
                
        bounds(1) = max([bounds(1) boundsellips(1)]);
        bounds(2) = min([bounds(2), boundsellips(2)]);
        
        optsConitope = jsrsettings('logfile',logFile,'verbose',opts.verbose-1,...
            'maxTime',(opts.jsr.time-(cputime-starttime)),'conitope.initProd',prodGrip,'conitope.tol',opts.jsr.tol);
        msg(logFile,opts.verbose>0,'\n>> JSR: Launching conitope\n')
        [boundsConitope, prodOptconi, prodVconi, infoconi]  = jsr_norm_conitope(M ,optsConitope);
        
        bounds(1) = max([bounds(1) boundsConitope(1)]);
        bounds(2) = min([bounds(2), boundsConitope(2)]);
        
        
        if nargout>1
            info.grip.bounds = boundsGrip;
            info.grip.prodOpt = prodGrip;
            info.grip.allProds = allProdsGrip;
            info.grip.info = infoGrip;
            
            info.ellips.bounds = boundsellips;
            info.ellips.P = Pellips;
            info.ellips.info = infoellips;
            
            info.conitope.bounds = boundsConitope;
            info.conitope.prodOpt = prodOptconi;
            info.conitope.prodV = prodVconi;
            info.conitope.info = infoconi;
        end
    end
 

end

if (~isstruct(bounds))
    msg(logFile,opts.verbose>0,'\n>>JSR: bounds on the jsr: [%.15g, %.15g]',bounds(1),bounds(2));
else
    msg(logFile,opts.verbose>0,'\n>>JSR: Please analyse further the relevant block sets in order to find jsr(M) = max(jsr(block{1}),...,jsr(block{q}))')
end

info.optsJSR = opts;

elapsedtime = cputime-starttime;
msg(logFile,opts.verbose>1,'\nEnd of JSR after %5.2f s',elapsedtime);

if (logFile ~= -1 && close)
    fclose(logFile);
end

if ischar(opts.saveEnd)
   if nargout>1
       save([opts.saveEnd '.mat'],'bounds','info')
   else
       save([opts.saveEnd '.mat'],'bounds')
   end
end
    
end
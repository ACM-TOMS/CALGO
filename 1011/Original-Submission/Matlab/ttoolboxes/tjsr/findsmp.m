function [ cand, nearlycand, info ] = findsmp( varargin )
% [ cand, nearlycand, info ] = findsmp( T, [algorithm], ['smaxp'|'sminp'], [options] )
% Searches for s.m.p.-candidates (s.min.p's and s.max.p's). (Search for s.min.p's is still experimental)
% Two algorithm types available: 'Gripenberg type algorithms' and the 'Genetic algorithm'
%
% Input:
% ======
%   T                       Cell array of square matrices of the same size.
%   algorithm               string, default='modifiedgripenberg', name of algorithm to use
%                           Implemented algorithms: 
%                               'gripenberg',
%                               'lowgripenberg', 'highgripenberg', 'randomgripenberg'
%                               'modifiedgripenberg' (default), 
%                               'bruteforce', 'necklacebruteforce',
%                               'genetic' (not recommended)
%   
% Options for both algorithms:
% ============================
%   'vpa'                           Uses vpa
%   'double'                        Tries to convert input to double prior computing
%   'verbose',val                   default=1, Defines the verbose level.
%   'nosimplify'                    Does not simplify products cand and nearlycand
%   'maxtime',val                   default=inf, Maximal time used for computation
%   'maxeval'                       default=inf, Maximum number of evaluations (approximate)
%   'bound',val                     default=empty, Searches until 
%                                       a) the bounds of the joint/lower spectral radius is in ( val(1), val(2) ); or
%                                       b) it is proven that the bounds cannot be fulfilled anymore
%   'nearlycandidate',val           default=0.99, Nearlycandidates must have spectral radius larger than val*bounds(1), 
%   'shortnearlycandidate',val      default=1, Removes all nearlycandidates whose ordering is longer than the val*maximal-length-of-candidates-ordering.
%                                       Note that 'genetic' algorithm does not search for nearly-candidates at the moment, thus this option has no effect for 'genetic' algorithm
%   'maxsmpdepth',val               maximal length of products which is searched for. Can be
%                                       arbitrary high (>100) for 'gripenbergmodifed' algorithm
%                                       high (<15) for 'gripenberg'
%                                       small (<12) for 'brute force'
%                                       Default value depends on algorithm used
%   
% 
% Options for Gripenberg type algorithms':
% =============================================
%   One can use either pre-defined options corresponding to one specific algorithm, and/or set all/any options by hand    
%   
%   'norm',functionhandle           Handle to a norm function. Default: @norm (i.e. 2-norm)
%   'rho',functionhandle            Handle to a spectral radius function. Default: @rho                     
%   'delta',val                     default=depends on algorithm, relative delta used in the Gripenberg Algorithm.
%   'N',val                         scalar or 1x3 vector of doubles, default=depends on algorithm, Number of kept products in each step
%                                       val(1)  ... number of products with smallest matrix norm kept
%                                       val(2)  ... number of products with largest matrix norm kept
%                                       val(3)  ... number of products randomly kept
%   'minsmpdepth',val               default=1, Minimal length of products
%   'nearlycanddelta',val           default=0.99, Maximal relative difference of spectral radius between candidates and nearly-candidates
%   'maxnumnearlycandidate',val     default=10, Maximum number of nearly candidates returned. If number of nearly candidates is larger, 'nearlycanddelta' is decreased
%   'sminp' | 'smaxp'               default='smaxp', Defines whether to search for s.min.p's or s.max.p's
%   'hardworking',val               default=1, sets 'maxsmpdepth' to 'hardworking' times <last-smpcandidate-length> each time a new candidate is found
%   'epsilon',val                   default=1e-10, epsilon used for comparing spectral radii
%
%       Pre-defined algorithms:
%           'gripenberg'/'grip'             Standard Gripenberg algorithm
%                                           Gives lower and upper bounds
%                                           Does not miss candidates
%                                           NN = [inf; inf; inf]; delta=0.95,
%
%           'lowgripenberg'/'lowgrip'       Modified Gripenberg algorithm keeping NN small products, delta=1,
%
%           'highgripenberg'/'highgrip'     Modified Gripenberg algorithm keeping NN large products, delta=1,
%
%           'modifiedgripenberg'/'modgrip'  Modified Gripenberg algorithm keeping NN/2 small and NN/2 large products, delta=1,
%                                           Gives good (usualy tight) lower, but bad upper bounds
%           'randomgripenberg'/'randgrip'   Modified Gripenberg algorithm keeping NNrandom products, delta=1,
%                                           Gives good lower, but bad upper bounds
%           'bruteforce'/'bf'               Brute force algorithm,
%                                           Very slow, does not miss candidates
%                                           Only algorithm which is proven to work for searching s.min.p's
%           'necklacebruteforce'/'nlbf'     Brute force algorithm, only computing products in different short necklace classes
%                                           Does not compute upper bounds
%                                   
% 
%
% Options for genetic algorithm:            
% ================================            
%               Written by Chia-Tche Chang, 2011. 
%               There is a bug and the returned value is sometimes wrongly normalized.                                           
%               Does not work for 's.min.p'
%       'popsize',val                   (default 300, minimum 10) population size.
%       'maxgen',val                    (default 1000) maximum number of generations.
%       'maxstall',val                  (default 15) maximum number of stalling iterations before increasing the maximum product length.
%       'maxtotstall',val               (default 100) maximum number of stalling iterations before terminating the algorithm.
%       'mutantprop',val                (default 0.3) mutation probability of a given product.
%       'muteprop',val                  (default 0.2) mutation proportion of a given product.
% 
% Output: 
% =======
%   cand                    Ordering of the products with highest spectral radius.
%   nearlycand              Orderings of products with spectral radius greater than val*bounds(1). Can be empty, depending on the algorithm.
%   info                    (struct) Additional info, depends on the algorithm used.
%                               info.time         time in seconds needed for computation
%                               info.jsrbound     (Interval) estimate on the JSR/LSR based on the spectral radius of the found candidates and norm estimates
%                               info.spectralgap  relative difference between info.jsrbound and second biggest eigenvalue found (from nearly-candidates)
%                               info.count        approximate number of computed matrices
%
% Note:
% ============
%   The Gripenberg algorithms are parallelised, the genetic algorithm is not.
%   There is a bug in the genetic algorithm and the returned value is sometimes wrongly normalized.
%
% E.g.: [ c, nc, info ] = findsmp( {[1 -1; 3 -2], [1 3; -1 -1]}, 'maxsmpdepth',15 )
%
% See also: tjsr
%
% Written by: tommsch, 2018

% Changelog: tommsch, 2019-11-15    Option 'sparse' removed
%                                   Option 'vpa' changed
%            tommsch, 2020-03-25    Totally rewritten

% XX Fuer modgrip selbst herausfinden lassen, wann vpa verwendet werden soll
% XX Legende fuer modgrip hinzufuegen
% XX searchonlyonecandidate works only for smaxp's at the moment
% XX bf can be implemented by gripenberg with delta=inf

%parse input
%%%%%%%%%%%%%%%%%%%%%
if( parsem('help',varargin) );
    varargin={{1},'help'}; end;
[~,varargin,method] =               parsem( {'modgrip','gripenberg','grip', ...
                                             'lowgripenberg','lowgrip','highgripenberg','highgrip','randomgripenberg','randgrip','modifiedgripenberg','modgrip', ...
                                             'bruteforce','bf','necklacebruteforce','nlbf', ...
                                             'genetic','gen'}, varargin );
[~,varargin,minmax] =               parsem( {'smaxp','sminp'}, varargin );
[nosimplify,varargin] =             parsem( 'nosimplify', varargin );
[vpaflag,varargin] =                parsem( 'vpa', varargin ); 
[doubleflag,varargin] =             parsem( {'double','d'}, varargin ); 
[verbose,varargin] =                parsem( {'verbose','v'}, varargin, 1, 'expect',@isnumeric ); %Verbose level
[shortnearlycand,varargin] =        parsem( {'shortnearlycandidate','shortnc'}, varargin, 1, 'expect',{'clcl',[0 inf]} ); %only keeps nearlycandidates whose length is less or equal than val times the length of the shortest candidate 
[maxnumnearlycand,varargin] =       parsem( {'maxnumnearlycandidate','maxnumnc'}, varargin, 10, 'expect',{'clcl',[0 inf]} ); %maximal number of nearlycandidates. If there are more, then delta is set to (delta+1)/2
[nearlycanddelta,varargin] =        parsem( {'nearlycandidate','nearlycanddelta','ncdelta','ncd','nc'}, varargin, .99, 'expect',{'opcl',[0 1]} ); %delta for nearlycandidates
[epsilon,varargin] =                parsem( {'epsilon','eps'}, varargin, 1e-10, 'expect',{'clcl',[0 1e-2]} ); 
[rhofun,varargin] =                 parsem( {'rho','r'}, varargin, @rho, 'expect',@(x) isa(x,'function_handle') );
[normfun,varargin] =                parsem( {'norm','n'}, varargin, [], 'expect',@(x) isa(x,'function_handle') || isempty(x) || isnumeric(x) );
[maxsmpdepth,varargin] =            parsem( 'maxsmpdepth', varargin, [], 'expect',@(x) isempty(x) || x>=1 ); %maximum length of products to be computed
if( ~isempty(normfun) && isnumeric(normfun) ); 
    normfun = @(x) norm( x, normfun ); end;
if( isempty(normfun) )
    switch method;
        case {'gripenberg','grip','lowgripenberg','lowgrip','highgripenberg','highgrip','modifiedgripenberg','modgrip','randomgripenberg','randgrip','bruteforce','bf'}; normfun = @norm;
        case {'necklacebruteforce','nlbf'}; normfun = @(M) inf; end; end;


M = varargin{1}; 
varargin(1) = [];
J = numel( M );
dim = size(M{1},1);

%preprocessing
%%%%%%%%%%%%%%%%%%%%
vprintf( 'Search candidate-smp: ', 'imp',[1 verbose] );
vprintf( '\n', 'imp',[2 verbose] );

if( vpaflag );
    Msym = cell( 1, J );
    for j = 1:numel( M )
        Msym{j} = sym( M{j} ); end; end;
    
if( doubleflag );
    for j = 1:J
        M{j} = double( M{j} ); end; end;

%select method
if( isequal(method,'genetic') || isequal(method,'gen') );
    if( ~isequal(minmax,'smaxp') );
        error( 'Genetic algorithm only works for ''smaxp''.'); end;
    [cand,nearlycand,info] = findsmp_genetic( M, maxsmpdepth, verbose, varargin{:} );
    o = [cand nearlycand];
    P = cellfun( @(x) tbuildproduct_fast(M,x), o, 'UniformOutput', false);
    R = cellfun( @rho, P );
    lb = info.jsrbound(1);
    ub = info.jsrbound(end);
    delta = [0; 0];
else;
    [NN,varargin] = parsem( 'N', varargin, [] ); %the number of kept products in each step equals ~NN
    if( isempty(NN) );
        NN = ceil( sqrt(dim*J)*10 ); end;
if( isempty(maxsmpdepth) );
        switch method
            case {'gripenberg','grip'}; maxsmpdepth = 30;
            case {'lowgripenberg','lowgrip'}; maxsmpdepth = 100;
            case {'highgripenberg','highgrip'}; maxsmpdepth = 100;
            case {'modifiedgripenberg','modgrip'}; maxsmpdepth = 100;
            case {'randomgripenberg','randgrip'};  maxsmpdepth = 100;
            case {'bruteforce','bf'}; maxsmpdepth = 10; 
            case {'necklacebruteforce','nlbf'}; maxsmpdepth = 15; end; end;    
    if( numel(NN)==1 );
        switch method
            case {'gripenberg','grip'}; NN = [inf; inf; inf]; 
            case {'lowgripenberg','lowgrip'}; NN = [NN; 0; inf]; 
            case {'highgripenberg','highgrip'}; NN = [0; NN; inf]; 
            case {'modifiedgripenberg','modgrip'}; NN = [ceil( NN/2 ); ceil( NN/2 ); inf];
            case {'randomgripenberg','randgrip'}; NN = [inf; inf; NN];
            case {'bruteforce','bf','necklacebruteforce','nlbf'}; NN = [inf; inf; inf]; end; end;
    [delta,varargin] = parsem( 'delta', varargin, 1 ); %delta from Gripenberg Algorithm        
    if( isrow(delta) );
        delta = delta.'; end;
    if( numel(delta) == 1);
        switch minmax;
            case {'sminp'}; delta = [delta; 0];
            case {'smaxp'}; delta = [0; delta]; end; end;
    switch method;
        case {'gripenberg','grip','lowgripenberg','lowgrip','highgripenberg','highgrip','modifiedgripenberg','modgrip','randomgripenberg','randgrip'}; delta_now = [inf; inf];
        case {'bruteforce','bf','necklacebruteforce','nlbf'}; delta_now = [0; 0]; end;
    choosematrix = @(Nt,Rt,lb,ub,delta_now) choosematrix_lhr( Nt, Rt, lb, ub, delta_now, delta, NN, epsilon );
        
    %computeNtRt
    computeNtRt = @(Pt,ot,lb) computeNtRt_grip( Pt, ot, lb, normfun, rhofun );
    
    %makePtot/computelbub
    switch method;
        case {'gripenberg','grip','lowgripenberg','lowgrip','highgripenberg','highgrip','modifiedgripenberg','modgrip','randomgripenberg','randgrip'}; makePtot = @makePtot_idx;
        case {'bruteforce','bf'}; makePtot = @makePtot_bf;
        case {'necklacebruteforce','nlbf'}; makePtot = @makePtot_nl; end;
    switch minmax
            case {'sminp'}; computelbub = @computelbub_sminp;
            case {'smaxp'}; computelbub = @computelbub_smaxp; end;
        
    [P,o,N,R,lb,ub,delta] = findsmp_worker( M, makePtot, computeNtRt, computelbub, choosematrix, delta_now, maxsmpdepth, verbose, varargin{:} ); end; %#ok<ASGLU>


%post processing
%%%%%%%%%%%%%%%%%%

%identify candidates and nearlycandidates
if( ~nosimplify )
    switch minmax
        case {'sminp'}; c_idx = R<=ub(end)+epsilon;
        case {'smaxp'}; c_idx = R>=lb(end)-epsilon; end;
   
    vprintf( 'Simplify %i candidates. ', nnz(c_idx), 'imp',[2 verbose] );    
    cand = simplify_ordering( o(c_idx) );
    while( true )
        switch minmax
            case {'sminp'}; nc_idx = ( R<ub(end)/nearlycanddelta & ~c_idx & cellfun( @numel, o )<=max( cellfun(@numel,cand) )*shortnearlycand );
            case {'smaxp'}; nc_idx = ( R>lb(end)*nearlycanddelta & ~c_idx & cellfun( @numel, o )<=max( cellfun(@numel,cand) )*shortnearlycand ); end;
        if( any(nc_idx) );
            vprintf( 'Simplify %i nearly-candidates. ', nnz(nc_idx), 'imp',[2 verbose] ); end;
        nearlycand = simplify_ordering( o(nc_idx) );
        if( numel(nearlycand)<maxnumnearlycand );
            break;
        else;
            nearlycanddelta = (nearlycanddelta+.5)/1.5; end; end; end;
vprintf( '\n', 'imp',[1 verbose] );

%compute spectralgap
try
    spectralgap = lb(end)/max( R(R<lb(end)/(1+epsilon)) );
catch
    spectralgap = min( R(R>ub(end)*(1+epsilon)) )/ub(end); end;

%Text output
vprintf( 'Candidates: \n%v\n', cand, 'imp',[2 verbose] );
if(~isempty(nearlycand) );
    vprintf( 'Nearly candidates: \n%v\n', nearlycand, 'imp',[2 verbose] ); end;
vprintf( '\n', 'imp',[2 verbose] )
vprintf( 'Bounds on the jsr : [%.15g, %.15g]\n', lb(end), ub(end), 'imp',[1 verbose] );
vprintf( 'Spectral gap: %.15g\n', spectralgap, 'imp',[1 verbose] );
if( any(delta) )
    vprintf( 'Relative Gripenberg-Delta: %.15g/%.15g\n', delta(1,end), delta(2,end), 'imp',[2 verbose] ); end;

%construct return variables
info.o = o;
info.jsrbound = [lb(end) ub(end)];
info.spectralgap = spectralgap;

end


function [P, o, N, R, lb, ub, delta_now] = findsmp_worker( M, makePtot, computeNtRt, computelbub, choosematrixfun , delta_now, maxsmpdepth, verbose, varargin  )
    
    %parse input
    [minsmpdepth,varargin] =            parsem( 'minsmpdepth', varargin, 1, 'expect',{'clop',[1 inf]} );
    [bound,varargin] =                  parsem( {'bound','searchonlyonecandidate','one','minJSR','sufficientbound','b'}, varargin, [], 'expect',@(x) numel(x)==2 || numel(x)==0 ); %stops algorithm after one candidate with spectral radius val has been found.
    [maxtime,varargin] =                parsem( 'maxtime', varargin, inf, 'expect',{'clcl',[0 inf]} );
    [maxeval,varargin] =                parsem( 'maxeval', varargin, inf, 'expect',{'clcl',[0 inf]} );
    [hardworking,varargin] =            parsem( 'hardworking', varargin, 1, 'expect',{'opop',[0 inf]} );
    parsem( varargin, 'test' );

    %initialize variables
    smpdepth = max( minsmpdepth-1, 0 );
    J = length( M );
    o = num2cell( mixvector(1:J,minsmpdepth), 1);
    P = cellfun( @(x) tbuildproduct_fast(M,x), o, 'UniformOutput',false );
    [N,R] = computeNtRt( P, o, 0 );
    [lb,ub] = computelbub( 0, inf, N, R, [1; 1] );
    [idx,delta_now] = choosematrixfun( N, R, lb, ub, delta_now );
    %start algorithm
    starttime = clock;
    numeval = 0;
    while( true )
        %check termination conditions
        smpdepth = smpdepth+1;
        if( ~any(idx) ); break; end;
        if( maxsmpdepth && smpdepth>=maxsmpdepth ); break; end; %maxsmpdept==0 means \infty
        if( ~isempty(bound) && (lb(end)>bound(1) && ub(end)<bound(2) || lb(end)>bound(2)  || ub(end)<bound(1)) ); break; end;
        if( etime(clock,starttime)>=maxtime ); break; end;
        if( numeval>=maxeval ); break; end;        

        %make new matrices and orderings
        [Pt,ot] = makePtot( M, P, o, idx );
        numeval = numeval + numel(Pt);

        %compute norm and rho
        vprintf( 'Number of matrices with length %i to compute: %6i,  \t', smpdepth, numel(Pt), 'imp',[2 verbose] );
        [Nt,Rt] = computeNtRt( Pt, ot, lb(end) );

        %compute bounds
        [lbt,ubt] = computelbub( lb(end), ub(end), Nt, Rt, delta_now(:,end) );
        if( lbt > lb(end) || ubt < ub(end) )
            vprintf( 'New bounds: [%.15g, %.15g]\n', lbt, ubt, 'imp',[1 verbose] ); 
            maxsmpdepth = max( [smpdepth*hardworking maxsmpdepth] ); end;

        %choose matrices
        [idx,delta_nowt] = choosematrixfun( Nt, Rt, lbt, ubt, delta_now(:,end) );        
        idx = [false(1,numel(N)) idx]; %#ok<AGROW>
        
        %add to results
        P = [P Pt]; %#ok<AGROW>
        o = [o ot]; %#ok<AGROW>
        N = [N Nt]; %#ok<AGROW>
        R = [R Rt]; %#ok<AGROW>
        lb = [lb lbt]; %#ok<AGROW>
        ub = [ub ubt]; %#ok<AGROW>
        delta_now = [delta_now delta_nowt]; %#ok<AGROW>
        vprintf( '\n', 'imp',[2 verbose] ); 
    end

end

function [ Pt, ot ] = makePtot_idx( M, P, o, idx )
    J = numel(M);
    if( issym(idx) )
        idx = isAlways( idx ); end;
    n = nnz(idx);
    %compute new matrices
    Pt = cell( 1, n*numel(M) );
    i = 0;
    for m = 1:J
        for p = find(idx);
            i=i+1;
            Pt{i} = M{m}*P{p}; end; end;
    ot = num2cell( [repmat([o{idx}],[1 J]); reshape(repmat(1:J,[n 1]),1,[])], 1 ); 
end

function [ Pt, ot ] = makePtot_bf( M, ~, o, ~ )
    ot = mixvector( 1:numel(M), max(cellfun(@numel,o))+1 );
    ot = num2cell( ot, 1 );
    Pt = cellfun( @(x) tbuildproduct_fast(M,x), ot, 'UniformOutput',false );
end

function [ Pt, ot ] = makePtot_nl( M, ~, o, ~ )
    ot = tgenNecklaces( max(cellfun(@numel,o))+1, numel(M) ); %TODO: computation of necklaces is expensive and should be accelerated
    ot = num2cell( ot, 1 );
    Pt = cellfun( @(x) tbuildproduct_fast(M,x), ot, 'UniformOutput',false );
end

function [ Nt, Rt ] = computeNtRt_grip( Pt, ot, lb, normfun, rhofun ); %#ok<INUSL>
    %compute norm and spectral radius
    exponent = 1./cellfun( @numel, ot );
    nPt = numel( Pt );
    if( nPt<200 )
        Nt = cellfun( normfun, Pt );
        Rt = cellfun( rhofun, Pt );
    else
        [Nt,Rt] = deal( zeros(1,numel(Pt)) );
        parfor i=1:numel(Nt);
            Nt(i) = normfun( Pt{i} ); %#ok<PFBNS>
            Rt(i) = rhofun( Pt{i} ); end; end;%#ok<PFBNS>
    Nt = Nt.^exponent;
    Rt = Rt.^exponent;
end 

function [ lb, ub ] = computelbub_smaxp( lb, ub, Nt, Rt, delta_now );
    %compute lower and upper bound
    lb = max( [lb, max(Rt)] );
    ub = min( [ub, max([ lb/delta_now(2) Nt])] );
end

function [ lb, ub ] = computelbub_sminp( ~, ub, ~, Rt, ~ );
    %compute lower and upper bound
    lb = 0;
    ub = min( [ub Rt] );
end

%choose matrix functions
%=======================
% Chooses products which will get children in the next round
% Input:
%   Nt/Rt               norms/spectral radii of last round
%   lb/ub               current lower/upper bound of JSR
%   delta_now           currently computed minimal relative delta, which is used to compute upper bound of JSR
%   delta               relative delta from Gripenberg algorithm
%   Nlow/Nhigh/Nrand    number of products with lowest/highest norm kept, and number of randomly matrices kept
% Output:
%   idx                 chosen products
%   delta_now           new minimal relative delta
%
function [ idx, delta_now ] = choosematrix_lhr( Nt, ~, lb, ub, delta_now, delta, Nlhr, epsilon ) %N-high-low-rand

    %preprocess
    delta_now = min( delta_now, delta );
    n = numel( Nt );
    
    %sort
    [Ntsort,ia] = sort( Nt );
    [~,ib] = sort( ia );
    
    %process upper part
    idxlast = Ntsort>lb*(1-epsilon)/delta(2); %XX Check if this is correct for sminp's
    idxlastbefore = idxlast;
    last = find( idxlast, 1 );
    
    if( isempty(last) );
        last = n; end;
    if( ~isinf(Nlhr(1)) && ~isinf(Nlhr(2)) );
        idxlast(last+Nlhr(1):n-Nlhr(2)) = false; end;
    if( Nlhr(3)<nnz(idxlast) );
        idxlast(randperm( nnz(idxlast), nnz(idxlast)-Nlhr(3) )+last-1) = false; end;    

    
    %process lower part
    idxfirst = Ntsort<ub/(1-epsilon)*delta(2); %XX Here is a bug
    first = find( idxfirst==0, 1, 'first' )-1;
    if( isempty(first) );
        first = n; end;    
    if( ~isinf(Nlhr(1)) && ~isinf(Nlhr(2)) );
        idxfirst(1+Nlhr(1):first-Nlhr(2)) = false; end;    
    if( Nlhr(3)<nnz(idxfirst) );
        idxfirst(randperm( nnz(idxfirst), nnz(idxfirst)-Nlhr(3))) = false; end;    
    
    %compute new delta
    if( ~isequal(idxlast,idxlastbefore) );
        val = find( idxlast==0, 1, 'last' )+1;
        delta_now(2) = minm( [delta_now(2) lb/minm(Ntsort(val:n)) delta(2)] ); end;    
    delta_now(1) = 0; %At the moment there is no knowledge about lower bounds
    
    idx = idxfirst | idxlast;
    
    %un-sort
    idx = idx(ib); %XX ia oder ib?

end



function [cand, nearlycand, info] = findsmp_genetic( M, maxsmpdepth, verbose, varargin )
%[cand, nearlycand, info] = tjsr_genetic(M, [options])
%
%Inputs:
%   M      -  a non-empty cell array containing square matrices of identical dimensions.
%
% Options:
%       'verbose',val    (default 1) defines the printing level:
%                        <  1 - no output is printed on the screen;
%                        >= 1 - prints the final lower bound and the associated product at the end of the algorithm;
%                        >= 2 - prints the current lower bound and the associated product at each generation of the algorithm;
%                        >= 3 - prints the current maximum product length and the number of stalling iterations at this current maximum product length at each generation.
%      'popsize',val       (default 150, minimum 10) population size.
%      'maxgen',val        (default 500) maximum number of generations.
%      'maxstall',val      (default 70) maximum number of stalling iterations before increasing the maximum product length.
%      'maxtotstall',val   (default 200) maximum number of stalling iterations before terminating the algorithm.
%                        
%      'mutantprop',val    (default 0.3) mutation probability of a given product.
%      'muteprop',val      (default 0.2) mutation proportion of a given product.
%       
% Output:
%   bound               a lower bound on the jsr of M.
%   c                   a product of matrices corresponding to this lower bound:
%                       the product A1*A2*A3 is expressed as [3 2 1] in this order.
%   nc                  Empty
%   info.elapsedtime    cpu time used to run the algorithm.
%   info.population     the complete population at the end of the algorithm.
%
% Written by Chia-Tche Chang, 2011
% Interface modified by tommsch

[searchonlyonecandidate,varargin] =         parsem( {'minJSR','searchonlyonecandidate','one','bound'}, varargin, 0, 'expect',@(x) numel(x)==1 ); %stops algorithm after one candidate has been found. To be used together with minJSR
[POPSIZE,varargin] =                        parsem( 'popsize', varargin, 500, 'expect',{'clop',[1,inf]} );
[MAXGEN,varargin] =                         parsem( {'maxgen'}, varargin, 200, 'expect',{'clop',[1,inf]} );
[MAXSTALLING,varargin] =                    parsem( {'maxstall','maxstalling'}, varargin, 70, 'expect',{'clop',[1,inf]} );
[MAXTOTALSTALLING,varargin] =               parsem( {'maxtotstall','maxtotalstalling'}, varargin, 1000, 'expect',{'clop',[1,inf]} );
[MUTANTPROP,varargin] =                     parsem( 'mutantprop', varargin, 0.3, 'expect',{'clcl',[0,1]} );
[MUTEPROP,varargin] =                       parsem( 'muteprop', varargin, 0.2, 'expect',{'clcl',[0,1]} );
[maxtime,varargin] =                        parsem( 'maxtime', varargin, inf, 'expect',{'clcl',[0,inf]} );
[maxeval,varargin] =                        parsem( 'maxeval', varargin, inf, 'expect',{'clcl',[0,inf]} );
starttime = clock; %starttime of tommsch

parsem( varargin, 'test' );

vprintf( 'The genetic algorithm has a bug, and sometimes reports spectral radii which are normalized wrongly!\nUse this algorithm with care.\n', 'cpr','err', 'imp',[0 verbose] );

% Initialization
STARTTIME = cputime; %starttime of Chiang
m = numel( M );
n = size( M{1}, 1 );
k = ceil( log(POPSIZE)/log(m+1) );
scaler = ((m+1).^(k-1:-1:0))';

% Cache
cache = findsmp_genetic_genCache( M, k, m, n );
ncache = size(cache, 2);
stalling = 0;
totalstalling = 0;
bound = searchonlyonecandidate^k;
bestpop = [];
for i = 2:ncache;
    value = max( abs(eig(mat(cache(:, i)))) );
    if( value > bound );
        bound = value;
        bestpop = zeros( 1, k );
        key = i-1;
        for j = 1:k;
            bestpop(j) = floor( key/scaler(j) );
            key = key - bestpop(j)*scaler(j); end; end; end;
bestpop = bestpop(bestpop ~= 0 );
bound = bound^(1/k); %XX Here is a bug. 
                        %T=overlapMat needs this line. 
                        %For a=1/12*[3 3 4 3 3 4 3 3 4 3 3]';M=-3;D=[-2 -1 0];S=getS('a',a,'M',M,'D',D);[T,Om] =transitionmatrix(S);U=constructU(T,1,'sym');T=restrictmatrix(T,U);
                        %the line is wrong.
X = eye( n );
testbestpop = bestpop;
for i = 1:length( testbestpop )-1;
    X = X * M{testbestpop(i)};
    testval = max(abs(eig(X)))^(1/i);
    if( testval >= bound );
        bestpop = testbestpop(1:i);
        bound = testval; end; end;
if( verbose >= 2 );
    fprintf(' \n');
    fprintf('Starting population: init lower bound on the JSR = %.15g with product: %s\n', tif( isempty(bestpop),0,bound), num2str(bestpop) ); end;

% Initial population
CURLENGTH = 2*k;
POPULATION = floor( (m+1)*rand(POPSIZE, CURLENGTH) );

% Genetic evolution
gen=0;
counter=0;
while(true)
    counter=counter+gen*MAXSTALLING;
    if( counter>maxeval); 
        break; end;
    gen=gen+1; 
    if( ~searchonlyonecandidate && gen>=MAXGEN ); 
        break; end;
    % Evaluation
    if( searchonlyonecandidate && ~isempty(bestpop) ); 
        break; end;
    if( etime(clock,starttime)>=maxtime ); 
        break; end;
    [nbeyes, idx] = sort( sum(POPULATION == 0, 2), 'descend' );
    POPULATION = POPULATION(idx, :);
    POPULATION(nbeyes == CURLENGTH, 1) = ceil( m*rand(1) );
    nbeyes = min( nbeyes, CURLENGTH-1 );
    fitness = zeros( POPSIZE, 1 );
    cached = floor( CURLENGTH/k );
    cachekey = zeros( POPSIZE, cached );
    for j = 1:cached;
        cachekey(:, j) = POPULATION(:, k*j-k+1:k*j)*scaler; end;
    cachekey = cachekey + 1;
    for i = 1:POPSIZE;
        X = eye(n);
        for j = 1:cached;
            X = X * mat(cache(:, cachekey(i, j)) ); end;
        for j = k*cached+1:CURLENGTH;
            if( POPULATION(i, j) ~= 0 );
                X = X * M{POPULATION(i, j)}; end; end;
        fitness(i) = max( abs(eig(X)) )^(1/(CURLENGTH-nbeyes(i))); end;
    [fitness, idx] = sort( fitness, 'descend' );
    POPULATION = POPULATION(idx, :);
    
    % Local optimization
    if( fitness(1) > bound );
        stalling = 0;
        totalstalling = 0;
        bound = fitness(1);
        bestpop = POPULATION(1, :);
        bestpop = bestpop(bestpop ~= 0 );
        testbestpop = bestpop;
        lenbestpop = length( testbestpop );
        CELLX = cell( lenbestpop + 1, 1 );
        X = eye(n);
        CELLX{1} = X;
        for i = 1:lenbestpop;
            X = X * M{testbestpop(i)};
            CELLX{i+1} = X;
            testval = max( abs(eig(X)) )^(1/i);
            if( testval > bound || (testval == bound && i < length(bestpop)) );
                bestpop = testbestpop(1:i);
                bound = testval; end; end;
        lenbestpop = length( testbestpop );
        X = eye( n );
        localimprove = 0;
        for i = lenbestpop:-1:1;
            testval = max( abs(eig(CELLX{i}*X)) )^(1/(lenbestpop-1));
            if( testval > bound );
                bestpop = testbestpop([1:i-1 i+1:end] );
                bound = testval;
                localimprove = localimprove + 1; end;
            for j = 1:m;
                testval = max( abs(eig(CELLX{i}*M{j}*X)) )^(1/lenbestpop);
                if( testval > bound );
                    bestpop = testbestpop;
                    bestpop(i) = j;
                    bound = testval;
                    localimprove = localimprove + 1; end; end;
            if( lenbestpop < CURLENGTH );
                for j = 1:m;
                    testval = max( abs(eig(CELLX{i+1}*M{j}*X)) )^(1/(lenbestpop+1));
                    if( testval > bound );
                        bestpop = [testbestpop(1:i)  j  testbestpop(i+1:end)];
                        bound = testval;
                        localimprove = localimprove + 1; end; end; end;
            X = M{testbestpop(i)}*X; end;
        if( localimprove );
            POPULATION = [bestpop, zeros(1, CURLENGTH-length(bestpop) ); POPULATION(1:end-1, :)]; end;
    else
        stalling = stalling + 1;
        totalstalling = totalstalling + 1; end;
    
    % Stopping criterion
    if( verbose >= 3 );
        fprintf( 'Generation #%3d: STA = %2d, LEN = %2d, lower bound = %.15g with product: %s\n', gen, stalling, CURLENGTH, tif( isempty(bestpop),0,bound), num2str(bestpop) );
    elseif( verbose >= 2 );
        fprintf( 'Generation #%3d:  current lower bound on the JSR = %.15g with product: %s\n', gen, tif( isempty(bestpop),0,bound), num2str(bestpop) ); end;
    if( ~searchonlyonecandidate && totalstalling >= MAXTOTALSTALLING );
        break; end;
    if( stalling >=  MAXSTALLING );
        stalling = 0;
        CURLENGTH = CURLENGTH + 1;
        if( CURLENGTH>maxsmpdepth )
            break; end;
        POPULATION = [POPULATION, zeros(POPSIZE, 1)]; end; %#ok<AGROW>
    
    
    % Selection and crossover
    NB_ELITE = max(2, min(3, floor(POPSIZE/50)) );
    NB_SPAWN = max(4, floor(POPSIZE/50) );
    NB_SWAP = min(POPSIZE - NB_ELITE - NB_SPAWN, floor(POPSIZE/2) );
    NB_MIX = POPSIZE - NB_ELITE - NB_SPAWN - NB_SWAP;
    ID_SWAP = ceil((POPSIZE/2)*rand(NB_SWAP, 2) );
    ID_MIX = ceil(POPSIZE*rand(NB_MIX, 2) );
    
    POP_ELITE = POPULATION(1:NB_ELITE,:);
    POP_SPAWN = ceil(m*rand(NB_SPAWN, CURLENGTH) );
    POP_SWAP = zeros(NB_SWAP, CURLENGTH);
    for i = 1:NB_SWAP;
        cut = ceil(CURLENGTH*rand(1) );
        POP_SWAP(i, :) = [POPULATION(ID_SWAP(i, 1), 1:cut) POPULATION(ID_SWAP(i, 2), cut+1:end)]; end;
    POP_MIX = zeros(NB_MIX, CURLENGTH);
    for i = 1:NB_MIX;
        mixer = floor(2*rand(CURLENGTH, 1))';
        POP_MIX(i, :) = POPULATION(ID_MIX(i, 1), :) .* mixer + POPULATION(ID_MIX(i, 2), :) .* (1-mixer); end;
    
    POPULATION = [POP_ELITE ; POP_SPAWN ; POP_SWAP ; POP_MIX];
    
    % Mutation
    MUTESTR = ceil(CURLENGTH * MUTEPROP);
    for i = 2:POPSIZE;
        if (rand(1) < MUTANTPROP);
            POPULATION(i, ceil(CURLENGTH * rand(1, MUTESTR))) = floor((m+1)*rand(1, MUTESTR) ); end; end;
end

% Termination
elapsedtime = cputime - STARTTIME;
%bestpop = tjsr_genetic_deperiod(bestpop); %we use the function simplify_ordering instead

% POST PROCESSING FOR tommsch-Interface GENERATION  %
if( isempty(bestpop) ); 
    bound=0; end;
cand = {flip(bestpop).'}; %change order and make it to column vector, since tommsch-programs need that
nearlycand = {};
info.time = elapsedtime;
info.population = POPULATION;
info.jsrbound = [bound inf];
info.count = counter;

if (verbose >= 1)
    fprintf('Algorithm terminated with lower bound on the JSR = %.15g with product: %s\n', bound, num2str(bestpop) ); end;

end

function cache = findsmp_genetic_genCache( M, k, m, n ) %  CACHE GENERATION  %
    if( k <= 1 );
        cache = zeros( n*n, m+1 );
        cache(:, 1) = vec( eye(n) );
        for i = 1:m;
            cache(:, i+1) = vec( M{i} ); end;
        return; end;
    cache = repmat( findsmp_genetic_genCache(M, k-1, m, n), 1, m+1 );
    gsize = (m+1)^(k-1);
    for i = 1:m;
        cache(:, gsize*i+1:gsize*(i+1)) = kron( eye(n), M{i}) * cache(:, gsize*i+1:gsize*(i+1) ); end;
end

function [ cand ] = simplify_ordering( cand )
    %Simplify candidates and remove duplicates
    if( ~isempty(cand) );
        cand=cellfun( @reducelength, cand, 'UniformOutput', 0 ); %reduce length
        LENGTHMAX = max( cellfun(@(x) max(size(x)),cand) )+1;   
        for i = 1:size( cand, 2 ); 
            cand{i}(LENGTHMAX,1) = 0; end;  
        cand = cell2mat( cand );                                
        cand(end,:) = [];                                     
        cand = unique( cand', 'rows' )';                        
        cand = num2cell( cand, 1 );
        cand = cellfun( @(x) removezero(x,1), cand, 'UniformOutput', 0 ); end;
end

function v = vec(M)
% Converts matrices to vectors
% Input: %   M       matrix
% Output: %   v       vector
    v = reshape(M,numel(M),1);
end

function M = mat(v,n)
% Input:    v       vector of length 2^n
%          [n]     (optional) the value of n.
%                  If n is not given, it is computed
% Output:   M       the converted matrix
if( nargin < 2);
    n = floor(sqrt(length(v)) ); end;
M = reshape(v,n,n);
end



function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   

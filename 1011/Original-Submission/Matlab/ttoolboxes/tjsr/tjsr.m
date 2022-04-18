 function [ JSR, type, alltype ] = tjsr( varargin )
% [ JSR, info, allinfo ] = tjsr( M, [options] ) 
% Computes the joint spectral radius of the set of matrices M.
%
% Input:
%   M       the set of matrices
%   
% Most important options:
%   'maxsmpdepth'               maximal length of s.m.p.-candidates to search for
%   'plot',val                  Plot-Output, val='norm', 'polytope','L', ...
%   'nearlycandidate',val       Relative difference between s.m.p.-candidates and nearly-s.m.p.s
%   'ordering',val              Ordering of a s.m.p.-candidate, given as COLUMN vector!
%   'nobalancing'               if no balancing shall be done
%   'balancingvector',val       Vector of as many entries as there are cyclic-roots (including from extra-vertices)
%   'invariantsubspace',val     Whether to search for invariant subspaces or not, val='none','joint','perm','trans'
%   'delta',val                 Accuracy of returned value. Default=1. If val<1 algorithm is faster, but gives only bounds.
%   'verbose'                   verbose level. Default=2
%
%   See the source-code of tjsr_parseinputtjsr for all available options.
%
% Data Output:
%   JSR             Interval containing the JSR or exact value of the JSR
%   info            a struct containing nearly all data which was generated during the computation
%
% E.g.: tjsr(tgallery('rand_rho',3,2,100),'plot','polytope','verbose',3)
%
% See also: tjsr_plotoutput, tjsr_parseminput
%
% Written by: tommsch, 2018
% For more information write to: <a href="tommsch@gmx.at">tommsch@gmx.at</a>

% Changelog: tommsch, 2019-05-27,   added option 'maxmemory'
%                                   less output for verbose level 2
%            tommsch, 2019-11-08,   Changed behaviour of option 'nobalancing'. Argument now is mandatory. If nobalancing=-1, balancing is done whenever a balancing vector was found.
%                                   Changed behaviour of option 'delta'. Option now does not disable balancing. Instead it sets <'nobalancing',-1>
%                                                                        If delta==-1 and validateupperbound>0, then min(type.lambda/type.opt.validateupperbound*1.02,1);

% Things to be done 
% XX Preproccesing must be done preworker, because otherwise there are problems when invariant subspaces are found                                                                                                                  
% XX Fix intermediate bounds!
% XX for points whose location could not be determined, we still can estimate the norm by norm(p)/smallest-diameter-of-the-polytope
% XX type.opt.nearlycandidate and type.opt.testspectralradius and type.opt.testeigenplane: Namen aendern zu eps...
% XX complex case machen
% XX See Jungers paper, about finiteness conjecture for binary matrices. Conditions under which one can eliminiate matrices from a set.
% XX Ich habe JSR(1) und lambda vermischt!
% XX fastnorm: Vertices die ausserhalb sind erst hinzufuegen, wenn sie auch getestet worden waeren.
% XX Ich estimate immer die gleichen vertices
% XX Test ob alle orderings Spaltenvektoren sind. Testen fuer: orderings von restarts, orderings berechnet, orderings uebergeben per parameter
% XX Alle matrizen transponieren, damit ich sparse arrays machen kann
% XX Test if tree corresponding to a candidate gets vertices added, and not those of nearly-candidates or extra-vertices. Otherwise restart and make nearly-candidates and extra-vertices smaller
% XX implement termination-criteria with eigenplanes, and its search for new smp-candidates
% XX "Duration tree (including proof)": andere Mitteilung wenn kein proof gemacht wird
% XX eps wann zwei smp gleichen spektralradius haben ueberdenken
% XX Test einfuegen fuer zero JSR
% XX merge options: 'nobalancing', 'balancingvector', 'balancing'.
% XX Add option specifiyng maximal number of threads used for norm computation. 
% XX If numberofthreads==1, run computation in mainthread
% XX Add estimate of remaining time for accuracy 0.9999

type = struct; %struct which contains all information
type.arguments_raw = varargin;
type.info.infotext = '';
type.info.errortext = '';

if( parsem('help',varargin) );
    varargin = {{[1]},'help'}; end;

% Test for most basic reasonable input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( nargin==0 );
    JSR = [];               
    type.info.errorcode = TJSR_INPUTERROR; 
    type.info.infotext = 'No input data'; 
    alltype = {struct}; 
    return; end;
if( isempty(varargin{1}) ); 
    JSR = [];               
    type.info.errorcode = TJSR_INPUTERROR; 
    type.info.infotext = 'Empty input data';
    alltype = {struct}; 
    return; end;
if( ~iscell(varargin{1}) ); 
    JSR = rho( varargin{1} ); 
    type.info.errorcode = TJSR_NOERROR; 
    return; end;

 %#ok<*LINPROG>
parg=varargin;

% Parse Options
%%%%%%%%%%%%%%%%%%
[val,parg] =                     parsem( 'clc', parg ); %clears the console before starting the algorithm
if( val );  
    clc; end;          
[type.opt.verbose,~] =           parsem( {'verbose','v'}, parg, 1, 'expect',@isnumeric );                 %defines the verbose-level: 0=no print output, 1=default, 2=lots of print output, 3=???
[type.opt.save,~] =              parsem( 'save', parg, 0, 'expect',{0 1} );                          %saves the struct type to a file after the algorithm terminated
[type.opt.maxnumrestart,parg] =  parsem( 'maxnumrestart', parg, 10, 'expect',{'opcl',[0 inf]} );                %maximum number of restarts
[type.opt.pauseonreset,parg] =   parsem( 'pauseonreset', parg, 0, 'expect',{0 1} );                  %pauses the algorithm and waits for a keypress after a reset has occured
[type.opt.proof,parg] =          parsem( 'proof', parg, 0, 'expect',{0 1} );                         %tests the constructed polytope at the end again. This test is done without any tricks to accelerate the computation. 
                                                                                %If val==1, the original matrices are used, if val==2, the processed/normalized matrices are used.
                                                                                %''proof'' works only if there are no invariant subspaces
[val,parg] =                     parsem( 'nopreprocess', parg );                      %Does not preprocess matrices (make positive, remove duplicates, add transposed, perturbate matrices) %Preprocess options
type.opt.preprocess=~val;
% Pre-processing steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
vprintf( 'once', -1 );
if( type.opt.preprocess); 
    parg{1} = preprocessmatrix( parg{1}, 'v',type.opt.verbose-1 ); end;  %preprocess matrices
alltype = cell(0); %cell which contains all returned structs  'type' from iterated runs
numrestart = 0; %counts how often the algorithm restarts

% Start
%%%%%%%%%%%%%%%%
while( true )
    
    % Start Preworker
    %%%%%%%%%%%%%%%%%%%%%
    type = tjsr_preworker( type, parg{:} );
    
    % Post-process Output
    %%%%%%%%%%%%%%%
    JSR = type.JSR;                   
    [JSR, parg, type, breakflag] = parseerrorcode( JSR, parg, type ); %parse error code and errors which can happen in this function
    if( nargout==3 ); 
        alltype{end+1} = type; end; %#ok<AGROW> %if all types are wanted, save them here.
    if( breakflag ); 
        break; end;
    if( type.opt.pauseonreset ); 
        type.info.infotext = vprintf( 'Restart. Press any key to continue\n', 'cpr','err', 'str',type.info.infotext ); 
        pause; end;
        
    numrestart = numrestart+1;
    if(numrestart>type.opt.maxnumrestart); 
        type.info.infotext = vprintf( 'Number of maximal restarts reached.\n', 'cpr','err', 'imp',[1 type.opt.verbose], 'str',type.info.infotext ); 
        break; end;
end

% Post-processing
%%%%%%%%%%%%%%%%%
type = orderfields( type ); %bring the elements in type in alphabetic order
type.opt = orderfields( type.opt ); %bring the elements in type.opt in alphabetic order
type.counter = orderfields( type.counter ); %bring the elements in type.counter in alphabetic order
type.cyclictree = orderfields( type.cyclictree ); %bring the elements in type.counter in alphabetic order
type.info = orderfields( type.info ); %bring the elements in type.counter in alphabetic order

if( type.opt.proof );
    try;    
        type = tjsr_proofpolytope( varargin{1}, type, type.opt.proof );
    catch;  
        type.info.infotext = vprintf( '\nProof failed. \n', 'cpr','err', 'imp',[0 type.opt.verbose], 'str',type.info.infotext );
        type.info.errortext = vprintf( '\nProof failed. \n', 'str',type.info.errortext,'npr'); end; end;

if( type.opt.save ); 
    save( ['tjsr_save_' datestr(now,'yyyy_mm_dd_HH_MM_SS') '.mat'], 'type' ); end; %save the output
tjsr_printerrortext( type, alltype );

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         TJSR_PREWORKER
% sets everything up for tjsr_worker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function type = tjsr_preworker(type, varargin)

% PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set some fields in type
%%%%%%%%%%%%%%%%%%%%%%%%%%
[M, type] = tjsr_parseminput( type, varargin ); %parse input and set type
type.info.infotext = vprintf( 'Preworker Start. Time: %s\n', datetime('now'), 'imp',[2 type.opt.verbose], 'str',type.info.infotext );
type.info.infotext = vprintf( 'Preprocessed matrices:\n%v\n', M, 'imp',[3 type.opt.verbose], 'str',type.info.infotext );
type = tjsr_checkoptions(type); %make tests if options are all known and sensible

type.info.matrixtype = identifymatrix( M );
if( type.info.matrixtype.nan );
    type.JSR = NaN;
    type.info.infotext = vprintf('Matrices contain NaNs.\n', 'imp',[1 type.opt.verbose], 'str',type.info.infotext ); 
    type.info.errorcode = TJSR_INPUTERROR;
    return; end;
if( ~type.info.matrixtype.finite );
    type.JSR=Inf;
    type.info.infotext=vprintf('Matrices contain Infs.\n', 'imp',[1 type.opt.verbose], 'str',type.info.infotext ); 
    type.info.errorcode=TJSR_NOERROR;
    return; end;
type.info.infotext=vprintf('Number of matrices: %i, Dimension: %i\n', type.counter.nummatrix, type.info.dim, 'imp',[2 type.opt.verbose], 'str',type.info.infotext ); 

if( type.counter.nummatrix==1 );
    type.JSR = trho(M{1});
    type.cyclictree.ordering = {[1]}; type.cyclictree.smpflag=[0]; type.cyclictree.orho=1;
    type.info.infotext = vprintf('JSR = %15.12g \n', type.JSR, 'imp',[1 type.opt.verbose], 'str',type.info.infotext ); 
    type.info.errorcode = TJSR_NOERROR;
    type.candidate = 1;
    return; end;

if( isempty(type.JSR) && isempty(type.cyclictree.ordering) ); %make rough estimate for JSR    
    type.JSR = estimatejsr( M ); %XX replace with findsmp
    type.info.infotext = vprintf('Rough estimate for JSR: %v\n',type.JSR, 'imp',[1 type.opt.verbose], 'str',type.info.infotext );
elseif( isempty(type.JSR) )
    type.JSR = [0 inf]; end; 

if( type.info.dim==1 ); %1-d case
    [type.JSR, val] = max( abs(cell2mat(M)) ); 
    type.cyclictree.ordering = {[val]}; 
    type.cyclictree.smpflag = [0]; 
    type.cyclictree.orho = 1;
    type.info.infotext = vprintf( 'JSR = %15.12g \n', type.JSR, 'imp',[1 type.opt.verbose], 'str',type.info.infotext ); 
    type.info.errorcode = TJSR_NOERROR; 
    return; end;

% Search for invariant subspaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[MM, basis] = invariantsubspace( M, type.opt.invariantsubspace, 'v',type.opt.verbose-2 );  %search for invariantsubspaces
type.counter.numblock = size( MM, 2 );
if( size(MM,2)>1 );
    type.basis = basis;
    type.info.infotext = vprintf( 'Invariant Subspaces found. Algorithm may behave wrong if matrices are badly scaled.\n', 'cpr',[.7 .4 .1]', 'imp',[2 type.opt.verbose], 'str',type.info.infotext );
    type.info.infotext = vprintf( 'Number of Blocks: %i\n', type.counter.numblock, 'imp',[2 type.opt.verbose], 'str',type.info.infotext );
    type.info.infotext = vprintf( 'Basis: \n%v\n', basis, 'imp',[3 type.opt.verbose], 'sze',[type.opt.showquantity, type.info.dim+type.counter.numblock], 'str',type.info.infotext );
    type.info.infotext = vprintf( 'Blocks: \n%v\n', MM, 'imp',[3 type.opt.verbose], 'sze',[type.opt.showquantity, type.info.dim+type.counter.numblock], 'str',type.info.infotext );
    type.block = cell( 1, type.counter.numblock );
    for m = 1:type.counter.numblock; 
        type.block{m}.JSR = estimatejsr( MM{m} ); end; %estimate JSR %XX replace with something more robust
    %set arguments for function call
    blockarg = varargin;     
    for m = 1:type.counter.numblock
        type.info.infotext = vprintf('=============================================\n=============================================\n=============================================\n', 'cpr',[0.0 0.3 0.3], 'imp',[3 type.opt.verbose], 'str',type.info.infotext );
        type.info.infotext = vprintf('=============================================\nCompute Block: %i',m,'cpr',[0.0 0.3 0.3], 'imp',[1 type.opt.verbose], 'str',type.info.infotext );
        blockarg{1} = MM{m}; %copy arguments, and replace matrices with the ones from the block
        val = max( cellfun(@(x) min(x.JSR),type.block) );
        if( max(type.block{m}.JSR)<val ); 
            type.info.infotext = vprintf( ' Block %i ommited due to estimate for JSR.\n', m, 'imp',[1 type.opt.verbose], 'str',type.info.infotext );
            continue; end;
        type.info.infotext = vprintf( '\n', 'imp',[1 type.opt.verbose], 'str',type.info.infotext );
        [~,blockarg] = parsem( 'validateupperbound', blockarg, val, 'set' ); %Algorithm terminates if upper estimate of JSR is less than validateupperbound.
        type.block{m} = type; 
        type.block{m} = rmfield( type.block{m}, 'block' );
        type.block{m} = tjsr_preworker( type.block{m}, blockarg{:} ); %XX I have the feeling that tjsr() should be called instead, since this is a totally new run
        type.info.infotext = vprintf( 'JSR of block %i: %v\n', m, type.block{m}.JSR, 'imp',[1 type.opt.verbose], 'str',type.info.infotext );
        if( type.block{m}.info.errorcode>0 );
            vprintf( '-------------------------------\nError in Block %i\n', m, 'cpr','err', 'imp',[1 type.opt.verbose], 'str',type.info.infotext );
            type.info.errorcode = TJSR_ERROR_INVARIANTSUBSPACE;
            break; end; end; 
    type = tjsr_settype_block( type ); %combine informations from blocks into type.
    if( type.info.errorcode~=TJSR_ERROR_INVARIANTSUBSPACE );
        type.info.errorcode = TJSR_NOERROR_INVARIANTSUBSPACE; end;
    return; end;

% MAKE CANDIDATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[type,flag] = tjsr_makecandidate( M, type );
if( flag ); 
    type.info.errorcode = TJSR_TOOMUCHCANDIDATE;
    return; end; %parse errorcode from makecandidate;
if( type.lambda==0 || ~isfinite(type.lambda)); 
    type.info.errorcode = TJSR_ZEROJSR; 
    return; end;
if(size(type.cyclictree.smpflag,2)==0); 
    type.info.errorcode=TJSR_NOCANDIDATEFOUND; 
    return; end;
for j = 1:type.counter.nummatrix; 
    M{j} = M{j}/type.lambda; end; %normalize matrices
type.M_normalized = M;

% BALANCING
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( type.counter.numordering==1);
    type.cyclictree.balancingvector=ones(size(type.cyclictree.ordering,2),1);
elseif( type.opt.nobalancing==1);
    type.info.infotext=vprintf('No balancing due to user input/default options.\n', 'imp',[2 type.opt.verbose], 'str',type.info.infotext );
    type.cyclictree.balancingvector=ones(size(type.cyclictree.ordering,2),1);
elseif(~isempty(type.cyclictree.balancingvector)); 
    %do nothing
elseif(anym(isnan([type.cyclictree.v0s{:}])));
    type.info.infotext=vprintf('No balancing due to orthogonal eigenplane/eigenvector.\n', 'imp',[2 type.opt.verbose], 'str',type.info.infotext );
    type.cyclictree.balancingvector=ones(size(type.cyclictree.ordering,2),1);
else
    type.info.infotext=vprintf('Balance %i Trees. ',type.counter.numordering, 'imp',[1 type.opt.verbose], 'str',type.info.infotext );
    type.info.infotext=vprintf('Balancing depth: %i - ',type.opt.balancingdepth, 'imp',[2 type.opt.verbose], 'str',type.info.infotext );
    type.info.infotext=vprintf('\n', 'imp',[3 type.opt.verbose], 'str',type.info.infotext );
    type.balancing=cell(1,type.counter.numordering);
    
    val=tliftproduct(M, type.opt.balancingdepth);
    for i=1:type.counter.numordering
        if( type.cyclictree.smpflag~=2)
            val2=arrayfun(@(x) tbuildproduct(M,type.cyclictree.ordering{i}(1:x,1))*type.cyclictree.v0{i}/type.cyclictree.orho(i)^x(end),0:length(type.cyclictree.ordering{i}(:,1)),'UniformOutput',0);
            val2=[val2{:}];
        else
            val2=type.cyclictree.v0{i}; end
        val2=cellfun(@(x) x*val2,val,'UniformOutput',0);
        type.balancing{i}.cyclictree.V{1}=[val2{:}];
        type.balancing{i}.cyclictree.v0{1}=type.cyclictree.v0{i};
        type.balancing{i}.cyclictree.v0s{1}=type.cyclictree.v0s{i}; end;
    [type.cyclictree.balancingvector, val]=tjsr_computebalancingvector(type);
    
    if( type.opt.nobalancing==-1 && (isempty(type.cyclictree.balancingvector) || val < 1)); 
        type.info.infotext=vprintf('No balancing vector found.\n', 'imp',[2 type.opt.verbose], 'str',type.info.infotext );
        type.cyclictree.balancingvector=ones(size(type.cyclictree.ordering,2),1);
    elseif(isempty(type.cyclictree.balancingvector) || val < 1)
        type.info.infotext=vprintf('No balancing vector found.\n','cpr','err', 'imp',[1 type.opt.verbose], 'str',type.info.infotext );
        type.info.errorcode=TJSR_NOBALANCINGVECTORFOUND; 
        return;
    else
        type.info.infotext=vprintf('Balancing vector found.\n', 'imp',[1 type.opt.verbose], 'str',type.info.infotext ); end;
    type.info.infotext=vprintf('\n\n==================================\n\n', 'imp',[3 type.opt.verbose], 'str',type.info.infotext,'cpr',[0.1 0.5 .2]); 
    type.info.infotext=vprintf('Balancing vector: %v\n',type.cyclictree.balancingvector, 'imp',[2 type.opt.verbose], 'str',type.info.infotext ); 
    if( type.opt.waitafterbalancing); 
        type.info.infotext=vprintf('Balancing finished. Press any key to continue.\n', 'imp',[1 type.opt.verbose], 'str',type.info.infotext );
        pause; end; end;
if(length(type.cyclictree.balancingvector)~=size(type.cyclictree.v0,2)); 
    type.info.infotext =vprintf('Length of balancing vector is not the same as number of orderings.\n','cpr','err', 'imp',[1 type.opt.verbose], 'str',type.info.infotext );
    type.info.errortext=vprintf('Length of balancing vector is not the same as number of orderings.\n','str',type.info.errortext,'npr'); end;
for i=1:min(type.counter.numordering,length(type.cyclictree.balancingvector));
    type.cyclictree.v0{i}=type.cyclictree.v0{i}*type.cyclictree.balancingvector(i); %balance v0 and v0s.
    type.cyclictree.v0s{i}=type.cyclictree.v0s{i}/type.cyclictree.balancingvector(i);
    type.info.infotext=vprintf('v0(%i):   \t%v\n',i,type.cyclictree.v0{i}, 'imp',[3 type.opt.verbose], 'str',type.info.infotext,'sze',[type.opt.showquantity,type.info.dim]); end;
for i=1:length(type.counter.numordering);
    type.info.infotext=vprintf('v0s(%i):  \t%v\n',i,type.cyclictree.v0s{i}, 'imp',[3 type.opt.verbose], 'str',type.info.infotext,'sze',[type.opt.showquantity,type.info.dim]); end;

% MAKE CYCLIC ROOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
type=tjsr_makecyclicroot(M,type); %makes the root of the trees. Sets V, Vs, o, normval, Rho, child

% MULTIPLY MATRICES WITH DELTA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(isequal(type.opt.delta,-1) && type.opt.validateupperbound>0);
    type.opt.delta=min(type.lambda/type.opt.validateupperbound*1.02,1); end;
for j=1:type.counter.nummatrix; 
    M{j}=M{j}*type.opt.delta; end %multiply matrices by delta

if( type.info.algorithm==TJSR_COMPLEXFUNCT); 
    type.info.errorcode = TJSR_COMPLEX; 
    return; end; %XX
type.counter.starttreetime = clock; %start_treetime
% START WORKER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
type = tjsr_worker(M, type);

% POST-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( type.opt.save); 
    type = tjsr_saveoutput(type,3); end;
type.counter.totaltime = etime(clock,type.counter.starttime);
type.counter.treetime = etime(clock,type.counter.starttreetime);
type.counter.treetime = etime(clock,type.counter.starttreetime); 
type.info.infotext = vprintf('\n', 'imp',[1 type.opt.verbose], 'str',type.info.infotext ); 
type.info.infotext = vprintf('Total time: %f s\n',type.counter.totaltime, 'imp',[2 type.opt.verbose], 'str',type.info.infotext ); 
type.info.infotext = vprintf('Time used for building the tree: %f s\n',type.counter.treetime, 'imp',[2 type.opt.verbose], 'str',type.info.infotext ); 
type.info.infotext = vprintf('Number of steps conducted/including steps in linprog: %i/%.3gM\n',type.counter.numstepsmall,type.counter.numstepbig/1000000, 'imp',[2 type.opt.verbose], 'str',type.info.infotext ); 
type.counter.numberofvertex = nnz([type.cyclictree.norm{:}]>1-type.opt.epspolytope);
type.info.infotext = vprintf('Number of vertices of polytope: %i\n',type.counter.numberofvertex, 'imp',[1 type.opt.verbose], 'str',type.info.infotext ); 
type.info.infotext = vprintf('Products which give lower bounds of JSR: \n%v\n',type.cyclictree.ordering{1}, 'imp',[1 type.opt.verbose], 'str',type.info.infotext ); 

end

function [type] = tjsr_worker(M, type)
%Worker uses only value type.lambda for computations.
%It must use type.JSR only for output stuff!

type.info.infotext=vprintf('Start constructing the polytope.\n\n', 'imp',[3 type.opt.verbose], 'str',type.info.infotext ); 

if( type.counter.numordering>1); 
    type = tjsr_recomputeroot(type); end; %recompute norms of all starting vertices, to find vertices which are inside %XX

if( type.lambda>type.JSR(1))
    type.JSR(1)=type.lambda; end;

type.info.infotext=vprintf('JSR = [ %15.12g, %15.12g ], norm= %13.12g, ', type.JSR(1), type.JSR(2), type.cyclictree.normlvl(end), 'imp',[1 type.opt.verbose], 'str',type.info.infotext );
while( true )
    % reset variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    leveltime=clock; %temporary variable for the clock-value
    
    % plot output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tjsr_plotoutput(type); %plot output  
    
    % check if we terminate the algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    val=whos('type'); val=val.bytes;
    if( val>type.opt.maxmemory); 
        type.info.errorcode = TJSR_OUTOFMEMORY; 
        break; end;
    if( type.opt.validatelowerbound+abs(type.opt.epspolytope)<type.lambda );      
        type.info.errorcode = TJSR_VALIDATELOWERBOUND;    
        break; end;
    if( type.opt.validateupperbound-abs(type.opt.epspolytope)>type.JSR(2) );      
        type.info.errorcode = TJSR_VALIDATEUPPERBOUND;    
        break; end;
    if( etime(clock,type.counter.starttime)>=type.opt.maxtime );                  
        type.info.errorcode = TJSR_MAXTIMEREACHED;        
        break; end;
    if( etime(clock,type.counter.starttreetime)>=type.opt.maxtreetime );          
        type.info.errorcode = TJSR_MAXTREETIMEREACHED;    
        break; end;
    if( type.counter.iteration>=type.opt.maxiteration );                          
        type.info.errorcode = TJSR_MAXITERATION;          
        break; end;
    if( type.counter.numstepsmall>=type.opt.maxstepnumber );                      
        type.info.errorcode = TJSR_MAXSTEPNUMBERREACHED;  
        break; end;
    if( nnz([type.cyclictree.norm{:}]>1-type.opt.epspolytope )>=type.opt.maxvertexnumber); 
        type.info.errorcode = TJSR_MAXVERTEXNUMBERREACHED;
        break; end;
    %Test for normal termination is done in the middle of this section
    
    % generate all valid new vertices and select the ones whose norms we compute
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    type = tjsr_generatenewvertex( M, type );   
    type = tjsr_estimatenorm( type );
   
    % select vertices of which we compute the exact norm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if( ~isnan(type.info.errorcode) ); 
        return; end;
    [ type, selectidx,  numselected, VV, VVidx, nVVbig ] = tjsr_selectvertex( type ); %choose vertices whose norm shall be computed. It is allowed to select parents and children.
    
    % test if we terminate or whether we have to change epspolytope
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(numselected(1)==0 && numselected(2)==0); %test if there are no vertices to test left
        if( type.opt.epspolytope>type.opt.epslinprog ); 
            type.JSR = type.lambda;
            type.info.errorcode = TJSR_NOERROR; 
            break;
        else
            type.opt.epspolytope = type.opt.epspolytope/2+2*type.opt.epslinprog; %if epspolytope is larger than epslinprog, we make epspolytope smaller and do not terminate
            type.info.infotext = vprintf('Set ''epspolytope'' to %i\n',type.opt.epspolytope, 'imp',[2 type.opt.verbose], 'str',type.info.infotext ); end; end;
    
    % compute norms
    %%%%%%%%%%%%%%%%%%
    type.info.infotext = vprintf( '#test: %i/%i, #V:%i/%i, ', numselected(1), sum(numselected), size(VV,2), nVVbig, 'imp', [1 type.opt.verbose], 'str', type.info.infotext );  
    vertextotest = [type.cyclictree.V{:}]; 
    vertextotest = vertextotest(:,selectidx);
    if( ~type.opt.alwaysout ); 
        [ normval, iterations ] = computepolytopenorm( vertextotest, VV, type.info.algorithm, type.opt.numcore, type.opt.epslinprog, type.opt.verbose, type.opt.solver );
    else
        normval = inf*ones(1,size(vertextotest,2));
        iterations = 0; end;
    
    % Save norms and other stuff 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    type.cyclictree.norm =      savetocellarray(normval,selectidx, type.cyclictree.norm);
    type.counter.numstepsmall = type.counter.numstepsmall+size(vertextotest,2);
    type.counter.numstepbig =   type.counter.numstepbig+sum(iterations);
    type.counter.iteration =    type.counter.iteration+1;
    type = tjsr_computemaxnormval( type, VVidx );
    
    type.cyclictree.livingvertex(end+1) = sum(numselected);
    type.cyclictree.timelvl(end+1) = etime(clock,leveltime);
    type.info.infotext = vprintf('\n', 'imp',[1 type.opt.verbose], 'str',type.info.infotext );  
    val=whos('type'); type.info.infotext=vprintf('Size of struct ''type'': %.2g MiB \n',val.bytes/1024/1024, 'imp',[3 type.opt.verbose], 'str',type.info.infotext );
    type.info.infotext = vprintf('numstepsmall/stepbig so far: %i/%.3gM \n',type.counter.numstepsmall,type.counter.numstepbig/1000000, 'imp',[3 type.opt.verbose], 'str',type.info.infotext );
    type.info.infotext = vprintf('Total number of vertices: %i \n',sum([type.cyclictree.L{:}]), 'imp',[3 type.opt.verbose], 'str',type.info.infotext );     
    type.info.infotext = vprintf('Time: %.1fs/%.1fs ',type.cyclictree.timelvl(end),etime(clock,type.counter.starttreetime), 'imp',[3 type.opt.verbose], 'str',type.info.infotext ); 
    type.info.infotext = vprintf('Vertices per tree: \n', 'imp',[3 type.opt.verbose], 'str',type.info.infotext ); 
    for iii = 1:type.counter.numordering; 
        if( type.cyclictree.L{iii}(end)~=0 || type.opt.verbose>=4); 
            type.info.infotext=vprintf('%r\n',removezero(type.cyclictree.L{iii},'right'), 'imp',[3 type.opt.verbose], 'str',type.info.infotext ); end; end;
    if( type.opt.naturalselection>0); 
        type.info.infotext=vprintf('JSR = [ %15.12g, %15.12g ], norm= %13.12g, ', type.JSR(1), type.JSR(2), type.cyclictree.normlvl(end), 'imp',[1 type.opt.verbose], 'str',type.info.infotext ); end;
    
    type.info.infotext=vprintf('\n\n', 'imp',[3 type.opt.verbose], 'str',type.info.infotext );  
    if( type.opt.save ); 
        type = tjsr_saveoutput( type, type.opt.save ); end; end;

end



function [JSR, parg, type, breakflag] = parseerrorcode(JSR, parg, type)
    verbose=parsem({'verbose','v'},parg,1);
    breakflag=0; %indicates that we make no restart of the algorithm anymore
    [~,parg]=parsem('JSR',parg,type.JSR,'set');   %since we want to reuse the old computations
    switch type.info.errorcode %parse errorcode      
        %STRANGE ERRORS
        case TJSR_INPUTERROR;
            type.info.infotext=vprintf('Input error.\n','cpr', 'err', 'imp',[0,verbose], 'str',type.info.infotext );
            breakflag=1;
        case TJSR_UNKOWNERROR;
            type.info.infotext=vprintf('Unkown error occured.\n', 'cpr','err', 'imp',[0,verbose], 'str',type.info.infotext );
            breakflag=1;
        case TJSR_STRANGEERROR
            type.info.infotext=vprintf('Abnormal program termination or other strange things happend.\n', 'cpr','err', 'imp',[0,verbose], 'str',type.info.infotext );
            breakflag=1;
            
        %WRONG TERMINATION; ERRORNUMBERS ARE POSITIVE        
        case TJSR_NOCANDIDATEFOUND;
            type.info.infotext=vprintf('No candidates found. If nothing works, add <searchonlyonecandidate>.\n','cpr','err', 'imp',[1,verbose], 'str',type.info.infotext );
            [~,parg]=parsem('maxsmpdepth',parg,ceil(type.opt.maxsmpdepth*1.2+10),'set'); %increase maxsmp
            [~,parg]=parsem('findsmp_N',  parg,type.opt.findsmp_N*2,'set');
            [~,parg]=parsem('minJSR',parg,[]);
            
        case TJSR_CANDIDATEISNOSMP
            type.info.infotext=vprintf('Candidate is no smp.\n','cpr','err', 'imp',[1,verbose], 'str',type.info.infotext );
            [~,parg]=parsem('maxsmpdepth',parg,ceil(type.opt.maxsmpdepth*1.2+10),'set');
            [~,parg]=parsem('findsmp_N',parg,10+type.opt.findsmp_N*2,'set');
            [~,parg]=parsem('minJSR',parg,[type.lambda+type.opt.epsequal JSR(2)],'set');
            [~,parg]=parsem('ordering',parg,[]);
            
        case TJSR_BETTERORDERINGFOUND 
            type.info.infotext=vprintf(['Product with higher spectral radius found. Length: %i, \nOrdering: %v,\n  rho=rho(candidate)*(1+%i)=%f\n'],...
                length(type.info.errorinformation{2}), type.info.errorinformation{2}, type.info.errorinformation{1}/type.lambda-1, type.info.errorinformation{1},'cpr','err', 'imp',[1,verbose], 'str',type.info.infotext );
            type.info.infotext=vprintf('Restart Algorithm with that ordering.\n','cpr','err', 'imp',[1,verbose], 'str',type.info.infotext );
            [~,parg]=parsem('smpflag',parg,[]);
            [~,parg]=parsem('v0',parg,[]);
            [~,parg]=parsem('v0s',parg,[]);
            [~,parg]=parsem('balancingvector',parg,[]);
            [~,parg]=parsem('maxsmpdepth',parg,10+length(type.info.errorinformation{2}),'set');
            [~,parg]=parsem('findsmp_N',parg,10+type.opt.findsmp_N*2,'set');
            [~,parg]=parsem('minJSR',parg,[type.info.errorinformation{1} JSR(2)],'set');
            [~,parg]=parsem('ordering',parg,{type.info.errorinformation{2}},'set');
            JSR=parsem('JSR',parg,[0 inf]); JSR(1)=type.info.errorinformation{1}; 
            [~,parg]=parsem('JSR',parg,JSR,'set');

        case TJSR_NOBALANCINGVECTORFOUND 
            if( type.opt.maxsmpdepth>50);
                type.info.infotext=vprintf('No balancing vector found. Set <''nobalancing'',1> and <''delta'',0.9999>.\n','cpr','err', 'imp',[1,verbose], 'str',type.info.infotext );
                type.info.errortext=vprintf('No balancing vector found. Set <''nobalancing'',1> and <''delta'',0.9999>.\n','str',type.info.errortext,'npr');
                [~,parg]=parsem('nobalancing',parg,1,'set');
                [~,parg]=parsem('delta',parg,.9999,'set');
                [~,parg]=parsem('nearlycandidate',parg,1,'set');
            else;
                type.info.infotext=vprintf('No balancing vector found.  Increase ''maxsmpdepth'', increase ''findsmp_N'', increase ''nearlycandidate''.\n','cpr','err', 'imp',[1,verbose], 'str',type.info.infotext );
                type.info.errortext=vprintf('No balancing vector found. Increase ''maxsmpdepth'', increase ''findsmp_N'', increase ''nearlycandidate''.\n','str',type.info.errortext,'npr');
                [~,parg]=parsem('maxsmpdepth',parg,ceil(type.opt.maxsmpdepth*1.1+2),'set');
                [~,parg]=parsem('findsmp_N',parg,ceil(3+type.opt.findsmp_N*1.3),'set');
                [~,parg]=parsem('nearlycandidate',parg,(type.opt.nearlycandidate+3)/4,'set');
            end
            
        case TJSR_COMPLEX 
            type.info.infotext=vprintf('Complex eigenvectors/extravertices. Algorithm not yet applicable.\n', 'imp',[1,verbose], 'str',type.info.infotext );
            breakflag=1;
            
        case TJSR_MAXTIMEREACHED 
            type.info.infotext=vprintf('''maxtime'' reached.\n','cpr','err', 'imp',[1,verbose], 'str',type.info.infotext );
            breakflag=1;
            
        case TJSR_MAXTREETIMEREACHED
            type.info.infotext=vprintf('''maxtreetime'' reached.\n','cpr','err', 'imp',[1,verbose], 'str',type.info.infotext );
            breakflag=1;
            
        case TJSR_MAXSTEPNUMBERREACHED
            type.info.infotext=vprintf('''maxstepnumber'' reached.\n','cpr','err', 'imp',[1,verbose], 'str',type.info.infotext );
            breakflag=1;
            
        case TJSR_MAXVERTEXNUMBERREACHED 
            type.info.infotext=vprintf('''maxvertexnumber'' reached.\n','cpr','err', 'imp',[1,verbose], 'str',type.info.infotext );
            breakflag=1;
            
        case TJSR_MAXTREEDEPTHREACHED 
            type.info.infotext=vprintf('''maxtreedepth'' reached.\n','cpr','err', 'imp',[1,verbose], 'str',type.info.infotext );
            breakflag=1;
            
        case TJSR_TESTEIGENPLANE 
            type.info.infotext=vprintf('Candidate is no smp. <Vs,newv> = %.15g > 1\n','cpr','err',type.info.errorinformation{1}, 'imp',[1,verbose], 'str',type.info.infotext );
            [~,parg]=parsem('maxsmpdepth',parg,ceil(type.opt.maxsmpdepth*1.2+10),'set');
            [~,parg]=parsem('findsmp_N',parg,type.opt.findsmp_N*2,'set');
            [~,parg]=parsem('minJSR',parg,[type.info.errorinformation{1} JSR(2)],'set');
            [~,parg]=parsem('ordering',parg,[]);
            
        case TJSR_MAXITERATION
            type.info.infotext=vprintf('''maxiterations'' reached.\n','cpr','err', 'imp',[1,verbose], 'str',type.info.infotext );
            breakflag=1;
            
        case TJSR_TOOMUCHCANDIDATE
            %take length such that at most 5 orderings are taken
            val=cellfun(@(x) nnz(x), type.cyclictree.ordering);
            val=sum(val<(1:max(val))',2);
            val=find(val>5,1);
            type.info.infotext=vprintf('A lot of candidates found. Maybe all orderings are smp''s. Set <''maxsmpdepth'',%i>.\n',val,'cpr','err', 'imp',[1,verbose], 'str',type.info.infotext );
            [~,parg]=parsem('maxsmpdepth',parg,val,'set');
            [~,parg]=parsem('epspolytope',parg,0,'set');
            [~,parg]=parsem('testspectralradius',parg,0,'set');
            %[~,parg]=parsem('proof',parg,1,'set');
            
        case TJSR_NOWORKER 
            %This option is actually not available anymore
            type.info.infotext=vprintf('Exit since <''numworker'',0> is set.\n', 'imp',[1,verbose], 'str',type.info.infotext );
            breakflag=1;
            
        case TJSR_ERROR_INVARIANTSUBSPACE 
            type.info.infotext=vprintf('  An error occured during computation of the JSR for an invariant subspace.\n  The algorithm may not recover from this problem.\n  If so, then probably because an option is used which does not work for matrices with joint invariant subspaces.\n  Investigate the (block-)matrices, use different options, or compute the JSR directly for the blocks.\n','cpr',[.6 .4 .1], 'imp',[1,verbose], 'str',type.info.infotext );
            for i=1:length(type.block); 
                if(isfield(type.block{i},'info') &&isfield(type.block{i}.info,'errorcode') && type.block{i}.info.errorcode>0); 
                    break; end; end; %read errorcode from block where the error occured %XX untested            
            [JSR, parg, type.block{i}, breakflag] = parseerrorcode(JSR, parg, type.block{i});
            
        case TJSR_ZEROJSR
            type.info.infotext=vprintf('JSR of candidate is maybe zero or not finite. This algorithm cannot handle this case.\n','cpr','err', 'imp',[1,verbose], 'str',type.info.infotext );
            breakflag=1;
            
        case TJSR_OUTOFMEMORY
            type.info.infotext=vprintf('Out of memory.','cpr','err', 'imp',[1,verbose], 'str',type.info.infotext );
            breakflag=1;            
        
        %CORRECT TERMINATION; ERRORNUMBERS ARE NEGATIVE
        case TJSR_VALIDATEUPPERBOUND 
            type.info.infotext=vprintf('JSR is less than validateupperbound=%i. Algorithm can terminate\n',type.opt.validateupperbound,'cpr',[0,0.5,0], 'imp',[1,verbose], 'str',type.info.infotext );
            breakflag=1; 
            
        case TJSR_VALIDATELOWERBOUND 
            type.info.infotext=vprintf('JSR is greater than validatelowerbound=%i. Algorithm can terminate\n',type.opt.validatelowerbound,'cpr',[0,0.5,0], 'imp',[1,verbose], 'str',type.info.infotext );
            breakflag=1;
            
        case {TJSR_NOERROR, TJSR_NOERROR_INVARIANTSUBSPACE, TJSR_EXACTVALUEFOUNDDURINGBALANCING }
            type.info.infotext=vprintf('Algorithm terminated correctly ','cpr',[0,0.5,0], 'imp',[1,verbose], 'str',type.info.infotext );
            switch type.info.errorcode
                case TJSR_EXACTVALUEFOUNDDURINGBALANCING;   type.info.infotext=vprintf('during balancing.\n','cpr',[0,0.5,0], 'imp',[1,verbose], 'str',type.info.infotext );
                case TJSR_NOERROR_INVARIANTSUBSPACE;        type.info.infotext=vprintf('in one invariant subspace..\n','cpr',[0,0.5,0], 'imp',[1,verbose], 'str',type.info.infotext );
                case TJSR_NOERROR;                          type.info.infotext=vprintf('\n','cpr',[0,0.5,0], 'imp',[1,verbose], 'str',type.info.infotext );
            end
            if( type.opt.epspolytope>0 && type.opt.delta>=1); 
                type.info.infotext=vprintf('Exact value found.\n','cpr',[0,0.5,0], 'imp',[1,verbose], 'str',type.info.infotext );
                JSR=JSR(1);
                type.JSR=JSR;
                type.info.infotext=vprintf('JSR = %15.12g\n',JSR(1), 'imp',[1,verbose], 'str',type.info.infotext );
            else; 
                type.info.errorcode=TJSR_NOERROR_APPROX;
                type.info.infotext=vprintf('Interval found.\n','cpr',[0,0.5,0], 'imp',[1,verbose], 'str',type.info.infotext );
                if(length(JSR)==1); JSR=[JSR JSR/type.opt.delta]; end;
                type.JSR=JSR;
                type.info.infotext=vprintf('JSR = [%15.12g, %15.12g]\n',JSR(1), JSR(2), 'imp',[1,verbose], 'str',type.info.infotext );
            end            
            breakflag=1;
            
        otherwise % THE VERY BAD TYPE Of TERMINATION
            type.info.infotext=vprintf('Unkown error occured. I have to quit.\n','cpr','err', 'imp',[1,verbose], 'str',type.info.infotext );
            breakflag=1;
            
    end;
end

%%%%%%%%%%%%%%%%%%%%%%%
% CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%
% ERR NUMBERS
%%%%%%%%%%%%%%%%%%%%%%%%%
function val = TJSR_NOWORKER;                       val=-80; end %worker was not started due to user input
function val = TJSR_VALIDATEUPPERBOUND;             val=-60; end %jsr is less than validateupperbound
function val = TJSR_VALIDATELOWERBOUND;             val=-50; end %jsr is larger than validatelowerbound
function val = TJSR_EXACTVALUEFOUNDDURINGBALANCING; val=-40; end 
function val = TJSR_NOERROR_INVARIANTSUBSPACE;      val=-20; end 
function val = TJSR_NOERROR;                        val=-10; end
function val = TJSR_NOERROR_APPROX;                 val=-5;  end
%
function val = TJSR_INPUTERROR;                     val=0;   end
function val = TJSR_UNKOWNERROR;                    val=inf; end
function val = TJSR_STRANGEERROR;                   val=nan; end
%
function val = TJSR_NOCANDIDATEFOUND;               val=10;  end
function val = TJSR_CANDIDATEISNOSMP;               val=20;  end %
function val = TJSR_NOBALANCINGVECTORFOUND;         val=30;  end
function val = TJSR_BETTERORDERINGFOUND;            val=60;  end
function val = TJSR_MAXTIMEREACHED;                 val=70;  end
function val = TJSR_MAXTREETIMEREACHED;             val=75;  end
function val = TJSR_MAXSTEPNUMBERREACHED;           val=80;  end
function val = TJSR_MAXVERTEXNUMBERREACHED;         val=90;  end
function val = TJSR_MAXTREEDEPTHREACHED;            val=100; end 
function val = TJSR_TESTEIGENPLANE;                 val=110; end %the stopping criterion has striken, thus the candidate is no smp
function val = TJSR_MAXITERATION;                   val=120; end
function val = TJSR_TOOMUCHCANDIDATE;               val=130; end %findsmp found too much candidates %this function is also defined in tjsr_makecandidate!!!
function val = TJSR_ERROR_INVARIANTSUBSPACE;        val=170; end %some error happend in an invariant subspace
function val = TJSR_ZEROJSR;                        val=180; end %JSR is zero. This algorithm cannot handle this case
function val = TJSR_OUTOFMEMORY;                    val=190; end
%
function val = TJSR_COMPLEX;                        val=1000;end


function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 
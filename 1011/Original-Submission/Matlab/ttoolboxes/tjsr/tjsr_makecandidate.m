function [type,toomuchcandidateflag] = tjsr_makecandidate(M, type);
% [type,toomuchcandidateflag] = tjsr_makecandidate(M, type);
% This function belongs to tjsr!
% Input:
% =====
%   The function reads the values in 
%   .opt.*  (in particular: .
%                           
%   .cyclictree.smpflag
%   .cyclictree.ordering
%   .type.cyclictree.v0
%   .cyclictree.v0s
%   
% Output:
% =======
%   computes the values type.
%       .ordering, .v0, .v0s, 
%       .numcandidate, .numnearlycandidate, .numextravertex, .numordering, 
%       .smpflag 
%       .info.algorithm
%   toomuchcandidateflag            set if there are too much candidates. In this case, this function may return early
%
% Behaviour:
% ==========
%   if v0 is given, then it is assumed that type.cyclictree.ordering is given already for multiple leading eigenvectors.
%   if v0 is empty, then otherwise
%   if v0 is given, then v0s must be given too
%   if smpflag is not given, then it is assumed all given orderings are smps
%   if multiplicity is not given, and not computed, it is set to one.
%
% Written by tommsch, 2018

% XX Test ob alle orderings Spaltenvektoren sind. Testen fuer: orderings von restarts, orderings berechnet, orderings uebergeben per parameter
% XX Teste spectralradius der Kandidaten. Dann kann ich mir smpflag sparen zu uebergeben.
% XX Wenn mehrere orderingens fuer einen cand gegeben sind, teste ob alle diese orderings den gleichen eigenvektor haben

    % Initialization
    %%%%%%%%%%%%%%%%%%%%%%

    toomuchcandidateflag = false;
    type.info.infotext = vprintf( 'Make smp-candidates. ', 'imp',[2 type.opt.verbose], 'str',type.info.infotext );
    type.info.infotext = vprintf( '\n', 'imp',[3 type.opt.verbose], 'str',type.info.infotext );
    dim = size( M{1}, 2 ); %the dimension
    if( iscell(type.cyclictree.smpflag) ); %test if smpflag is a vector (and not a cell array)
        error('''smpflag'' must be given as a vector.' ); end; 
    
    oo_mult = {};
    noo_mult = {};
    oo_single = {};
    noo_single = {};
    
    % Pre-processing
    %%%%%%%%%%%%%%%%%%
    if( ~isempty(type.cyclictree.ordering) && isempty(type.cyclictree.smpflag) ); %if smpflag is not given and orderings are not empty, then all orderings are candidates
        type.cyclictree.smpflag = zeros( 1, size(type.cyclictree.ordering,2) ); end; 
    
    % compute: oo_mult, v0_mult, v0s_mult, multiplicity
    %          noo_mult, nv0_mult, nv0s_mult, nmultiplicity
    %          
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   We first store everything in temporary variables: oo_single/_mult, noo_single/_mult, v0mult, v0smult, extravertex_, ... .
    
    
    if( isempty(type.cyclictree.ordering) ) %if no ordering is given, compute orderings (i.e. candidates and nearlycandidates)
        [oo_single, noo_single, type.info.findsmp] = findsmp( M, 'maxsmpdepth',type.opt.maxsmpdepth, 'v',type.opt.verbose-1, 'N',type.opt.findsmp_N, 'minJSR',type.opt.minJSR, 'nearlycandidate',type.opt.nearlycandidate ); 
        if( size(oo_single,2)>type.opt.maxnumcandidate*type.info.dim ); %if too much candidates were found
            toomuchcandidateflag = true;
            return; end; 
        
    elseif( isempty(type.cyclictree.v0) ) %no eigenvectors are given
        oo_single = type.cyclictree.ordering(type.cyclictree.smpflag==0); %orderings are given, but no eigenvectors. We still have to compute multiple eigenvectors v0 and v0s for each orderings
        noo_single = type.cyclictree.ordering(type.cyclictree.smpflag==1); 
        if( min(cellfun(@(x) size(x,2),oo_single))>max(cellfun(@(x) size(x,1),oo_single))^type.counter.nummatrix );
            type.info.infotext = vprintf( 'The array of orderings is probably transposed. Each ordering must be a column!\n', 'cpr','err', 'imp',[0 type.opt.verbose], 'str',type.info.infotext ); end;
        
    elseif( ~isempty(type.cyclictree.v0) && isempty(type.cyclictree.v0s) ) %orderings are given, v0 is given, v0s is not given. We set v0s to zero.  We do not have to compute eigenvectors any more.
        oo_mult = type.cyclictree.ordering(type.cyclictree.smpflag==0); 
        noo_mult = type.cyclictree.ordering(type.cyclictree.smpflag==1);    
        v0_mult = type.cyclictree.v0(type.cyclictree.smpflag==0); 
        nv0_mult = type.cyclictree.v0(type.cyclictree.smpflag==1); 
        v0s_mult = v0_mult; %make v0s_mult the same size as v0_mult
        [v0s_mult{:}] = deal( zeros(dim,1) );
        nv0s_mult = noo_mult; 
        [nv0s_mult{:}] = deal(zeros(dim,1));
        if( isempty(type.cyclictree.multiplicity) ); 
            multiplicity = ones( 1, sum(type.cyclictree.smpflag==0) ); 
            nmultiplicity = ones( 1, sum(type.cyclictree.smpflag==1) ); end;
        if( min(cellfun(@(x) size(x,2),oo_mult))>max(cellfun(@(x) size(x,1),oo_mult))^type.counter.nummatrix );
            type.info.infotext = vprintf( 'The array of orderings is probably transposed. Each ordering must be a column!\n', 'cpr','err', 'imp',[0 type.opt.verbose], 'str',type.info.infotext ); end;
            
    else %orderings are given and leading eigenvectors. We do not have to compute eigenvectors any more.
        oo_mult = type.cyclictree.ordering(type.cyclictree.smpflag==0); 
        noo_mult = type.cyclictree.ordering(type.cyclictree.smpflag==1);    
        v0_mult = type.cyclictree.v0(type.cyclictree.smpflag==0); 
        nv0_mult = type.cyclictree.v0(type.cyclictree.smpflag==1); 
        v0s_mult = type.cyclictree.v0s(type.cyclictree.smpflag==0);     
        nv0s_mult = type.cyclictree.v0s(type.cyclictree.smpflag==1);    
        if( isempty(type.cyclictree.multiplicity)); 
            multiplicity = ones( 1, sum(type.cyclictree.smpflag==0) ); 
            nmultiplicity = ones( 1, sum(type.cyclictree.smpflag==1) ); end;
        if( min(cellfun(@(x) size(x,2),oo_mult))>max(cellfun(@(x) size(x,1),oo_mult))^type.counter.nummatrix );
            type.info.infotext = vprintf( 'The array of orderings is probably transposed. Each ordering must be a column!\n', 'cpr','err', 'imp',[0 type.opt.verbose], 'str',type.info.infotext ); end; end;
   
    
    %compute multiple eigenvectors if needed
    if( isempty(oo_mult) ); %if we still need to compute (multiple) eigenvectors
        type.info.infotext = vprintf( 'Compute eigenvectors. \n', 'imp',[3 type.opt.verbose], 'str',type.info.infotext );
        [oo_mult, v0_mult, v0s_mult, multiplicity] = leadingeigenvector( M, oo_single, ...
                               'epsequal',type.opt.epsequal, ...
                               'nomultipleeigenvector', type.opt.nomultipleeigenvector, ...
                               'v',type.opt.verbose-1, ...
                               'complexeigenvector',type.opt.complexeigenvector, ...
                               'cycle',1 );
        [noo_mult, nv0_mult, nv0s_mult, nmultiplicity] = leadingeigenvector( M, noo_single, ...
                               'epsequal',type.opt.epsequal, ...
                               'nomultipleeigenvector', type.opt.nomultipleeigenvector, ...
                               'v',type.opt.verbose-1, ...
                               'complexeigenvector',tif(isequal(type.info.algorithm,TJSR_COMPLEXFUNCT),type.opt.complexeigenvector,3), ...
                               'cycle',1 ); end;
    
    
    if( type.opt.nearlycandidate==0 );  %Remove nearly candidates if option is given
        noo_mult = {}; 
        nv0_mult = {}; 
        nv0s_mult = {}; 
        nmultiplicity = []; end;
    
    % Check if eigenvectors are positive for cone-case
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v0property = identifymatrix( [v0_mult nv0_mult type.cyclictree.extravertex] );
    if( type.info.matrixtype.nonneg && ~v0property.nonneg ); %Even if the matrices are nonnegative, there could be multiple leading eigenvectors, of which one has negative entries.
        if( isequal(type.info.algorithm,0) || isempty(type.info.algorithm) )
            if( isempty(type.info.algorithm) ); 
                type.info.infotext = vprintf( 'Warning: Nonnegative matrices with non-nonnegative leading eigenvector. Non-nonnegative eigenvectors are removed and|or rounded towards zero.\nAlgorithm (P) is used (cone-norm). To remove this warning add <''algorithm'',''P''/''R''/''C''> as option.\n', 'cpr','err', 'imp',[1 type.opt.verbose], 'str',type.info.infotext ); end;
            for i = 1:size( v0_mult, 2 ); 
                v0_mult{i}(v0_mult{i}>-10*eps & v0_mult{i}<0)=0; 
                v0s_mult{i}(v0s_mult{i}>-10*eps & v0s_mult{i}<0)=0; end;
            idx = cellfun( @(x) any(x<0),v0_mult );
            v0_mult(idx) = []; 
            v0s_mult(idx) = []; 
            oo_mult(idx) = []; 
            multiplicity(idx) = [];
            for i = 1:size( nv0_mult, 2 ); 
                nv0_mult{i}(nv0_mult{i}<-10*eps & nv0_mult{i}<0) = 0; 
                nv0s_mult{i}(nv0s_mult{i}<-10*eps & nv0s_mult{i}<0) = 0; end;
            idx = cellfun( @(x) any(x<0), nv0_mult ); 
            nv0_mult(idx) = []; 
            nv0s_mult(idx) = []; 
            noo_mult(idx) = []; 
            nmultiplicity = []; end; end;
    
    % Select eigenvectors
    %%%%%%%%%%%%%%%%%%%%%
    type.info.infotext = vprintf( 'Set to choose from (these are not the final candidates!):\n', 'imp',[3 type.opt.verbose], 'str',type.info.infotext );
    type.info.infotext = vprintf( 'Maximum multiplicity of leading eigenvectors: %i\n', max(multiplicity), 'imp',[3 type.opt.verbose], 'str',type.info.infotext );
    type.info.infotext = vprintf( 'Candidates %v\n', oo_mult, 'imp',[3,type.opt.verbose], 'str',type.info.infotext);
    type.info.infotext = vprintf( 'Classify leading eigenvectors. \n', 'imp',[3 type.opt.verbose], 'str',type.info.infotext );
    %XX maybe we have to make sure that all orderings in oo_mult and noo_mult are column vectors
    [v0, v0s, oo, smpflag, mult, type] = tjsr_selectordering( [v0_mult nv0_mult], [v0s_mult nv0s_mult], [oo_mult noo_mult], [zeros( 1, length(v0_mult) ) ones( 1, length(nv0_mult) )], [multiplicity nmultiplicity], type ); %select smp-candidates
    orho = cellfun( @(x) max(abs(eig(tbuildproduct(M, removezero(x(:,1),'all')))))^(1/length(removezero(x(:,1),'all'))), oo ); 
    orho = orho./max( orho ); 
    orho(smpflag==0) = 1;
    
    
    type = tjsr_setproblemtype( type, v0property ); %sets the value of type.info.algorithm
    
    %compute spectral radius of candidate
    % XX I think, this should be made outside of this function
    if( isempty(oo_mult) );
        type.lambda = [];
        return;
    else
        type.lambda = trho( tbuildproduct(M,oo_mult{1}(:,1)))^(1/length(oo_mult{1}(:,1)) ); end;
    if( type.lambda==0 ); 
        %this case is handled outside of this function
        return; end;
    
    % Add Extravertices
    %%%%%%%%%%%%%%%%%%%%%%%
    if( ~isempty(type.cyclictree.extravertex) && ~iscell(type.cyclictree.extravertex) ); %if extravertices are given with the option 'extravertex'
        error('''extravertex'' must be given as a cell array of column vectors.'); end;
    extravertex_ = [type.cyclictree.extravertex type.cyclictree.v0(type.cyclictree.smpflag==2)];
    autoextravertex_ = [];
    if( type.opt.autoextravertex );
        autoextravertex_ = extravertex( [extravertex_{:}], M, oo, v0, ...
                                        't',type.opt.autoextravertex, 's',1/type.lambda, 'pt',type.info.algorithm, 'c', 'v',type.opt.verbose ); end;
    extravertex_ = [extravertex_ autoextravertex_];
    nextravertex_ = numel( extravertex_ );
    type.cyclictree=rmfield(type.cyclictree,'extravertex');
    
    % Post-processing
    %%%%%%%%%%%%%%%%%
    type.cyclictree.orho =           [orho NaN*ones( 1, nextravertex_ )];
    type.cyclictree.smpflag =        [smpflag 2*ones( 1, nextravertex_ )];
    type.cyclictree.ordering =       [oo repcell( [], 1, nextravertex_ )];
    type.cyclictree.v0 =             [v0 extravertex_];
    type.cyclictree.v0s =            [v0s num2cell( zeros(type.info.dim, nextravertex_), 1 )];
    type.cyclictree.multiplicity =   [mult ones( 1, nextravertex_ )];
    type.cyclictree.oclass =         type.cyclictree.ordering;
    type.cyclictree.maxlengthordering = cellfun ( @(x) size(x,1), type.cyclictree.ordering );
    
    type.counter.numcandidate = nnz( type.cyclictree.smpflag==0 );
    type.counter.numnearlycandidate = nnz( type.cyclictree.smpflag==1 );
    type.counter.numextravertex = nnz( type.cyclictree.smpflag==2 );
    type.counter.numordering = length( type.cyclictree.smpflag );    
    
    type.info.infotext = vprintf( 'Selected candidates (each column is one candidate)\n%v\n', type.cyclictree.ordering(type.cyclictree.smpflag==0), 'imp',[2 type.opt.verbose], 'str',type.info.infotext );
    if( ~isempty(noo_mult)); 
        type.info.infotext = vprintf( 'Selected nearlycandidates: \n%v\n', type.cyclictree.ordering(type.cyclictree.smpflag==1), 'imp',[2 type.opt.verbose], 'str',type.info.infotext ); end;
    type.info.infotext = vprintf( 'Number of candidates (with multiples)/ nearlycandidates (with multiples) / extravertices: %i / %i / %i\n', ...
                                  nnz(type.cyclictree.smpflag==0), nnz(type.cyclictree.smpflag==1), nnz(type.cyclictree.smpflag==2),...
                                  'imp',[2 type.opt.verbose], 'str',type.info.infotext );
    type.info.infotext = vprintf( 'v0: \n%v\n',v0, 'imp',[3 type.opt.verbose], 'str',type.info.infotext );
    type.info.infotext = vprintf( 'v0s: \n%v\n',v0s, 'imp',[3 type.opt.verbose], 'str',type.info.infotext );    
     
    if( (any(cellfun(@(x) any(isnan(x)),type.cyclictree.v0s)) || any(type.cyclictree.multiplicity>1)) && type.opt.testeigenplane>-inf ); 
        type.info.infotext = vprintf( 'Eigenplane-termination-criterium cannot be used.\n', 'imp',[1 type.opt.verbose], 'str',type.info.infotext );
        type.opt.testeigenplane = -inf; end;
    
    type.info.infotext = vprintf( '\n', 'imp',[2 type.opt.verbose], 'str',type.info.infotext );
    
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 
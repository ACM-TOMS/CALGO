function varargout = blf( varargin )
% [ c, PM, xyzv, oo ] = blf( [oo], S, [options] )
% Plots the basic limit function of multiple subdivision schemes.
%
% Input: 
%       [oo]                        Defines the ordering in which the subdivision operators are applied.
%                                   If S consists of more than one subdivision operator, this argument must be given.
%       S                           This argument is passed to getS and getS must return a cell array of subdivision operators.
%
% Options:
%       'diff',val                  1 x dim - vector or integer, computes partial derivatives (finite differences)
%       'plot',cellarray            cell array of arguments passed to plotm
%                                   if cellarray==0, nothing is plotted.
%       'verbose',val               verbose level
%       'start'                     starting sequence
%       'iteration',val             stops after val iterations. Default: depends on the subdivision operators S
%       'maxiteration',val          (option subject to be changed) Maximum number of iterations. Default=50
%       'numpoint',val              (option subject to be changed) Target-number of points in the sequence to be computed. Default=30000
%       'maxnumpoint',val           (option subject to be changed) Maximum number of points in the sequence to be computed. Default=300000
%       'removezero',val            default: 1, Removes zeros from the output sequence. Values greater equal 2 may lead to strange behaviour and should not be used.
%       'wavelet'                   (experimental) Uses the mask of the first subdivision scheme (w.r.t. oo) with alternating signs as starting sequence.
%       'oo',val                    
%
% Output:
%       c                           sequence of the basic limit function.
%       PM                          (experimental) dilation matrix for returned mesh. Returned matrix may be transposed due to some bug
%       xyzv                        grid of function values v to the corresponding point xyz in \RR^s
%       oo                          ordering which was used
%
% E.g.: [ c, PM, xyzv, oo ] = blf('2_butterfly','iteration',3,'diff',1)
%
% See also: tile, getS
%
% Written by: tommsch, 2016
% For more information write to: <a href="tommsch@gmx.at">tommsch@gmx.at</a>


% XX There is something wrong with the ordering of multiple schemes.
% XX Probably PM must be transposed somewhere
% XX Output of [c,PM,x]=blf(multiplyS(S2,S1),'iteration',1);
%    and       [c12,PM12,x12]=blf({[],[1 2]},[S2;S1],'iteration',2);
%    is not the same
% XX Fuer scalar-valued subdivision und vectorvalued data verallgemeinern
% XX Fuer periodische verallgemeinern
% XX Add option ''periodic''

%#ok<*ALIGN>

% Parse input
[verbose,varargin] = parsem( {'verbose','v'}, varargin, 1 );
[plotval,varargin] = parsem('plot',varargin,{});
[diff,varargin] = parsem('diff',varargin,0);
[c,varargin]=parsem( {'start','c'},varargin,[]);
[removezero,varargin] = parsem('removezero', varargin, 1 );
[iteration,varargin] = parsem( {'iteration','it'}, varargin, [] ); 
[maxiteration,varargin] = parsem( {'maxiteration','maxit'}, varargin, 60 );
[numpoint,varargin]= parsem( {'numpoint','npt'}, varargin, 30000 );
[maxnumpoint,varargin] = parsem( {'maxnumpoint','maxnpt'}, varargin, [], 'expecte',{'opop',[0 inf]} );
% [periodic,varargin] = parsem('periodic',varargin);
[wavelet,varargin] = parsem( 'wavelet', varargin );
[oo,varargin] = parsem( {'ordering','oo','o'}, varargin, [] );
[vectordim,varargin] = parsem( {'vdim','vectordim'}, varargin, 1 );


if( size(varargin,2)>=2 && ~isS(varargin{1}) );
    S = getS( varargin{2} ); 
    oo = constructordering( varargin{1} ); 
    varargin(1:2) = [];
elseif( ~isempty(oo) )
    S = getS( varargin{1} );
    oo = constructordering( oo );
    varargin(1) = [];
else
    S = getS( varargin{1} );
    oo = constructordering( size(S,1), 0, 'random',100 );    
    varargin(1) = [];
end


a = S(:,1); 
M = S(:,2);
dim = size( M{1}, 1 );     %dimension of the hypersurface
if( isempty(c) ); 
    if( wavelet );
        idx = ordering2vector( oo, 1 );
        c = a{idx};
        val = sizem( c.c );
        val = arrayfun( @(x) 1:x, val, 'UniformOutput',false );
        val = (-1).^(convm( val{:}, 'outer' ));
        c.c = (c.c).*val;
    else;
        c = sequence( 1, zeros(dim,1) ); end;
elseif( ~isa(c,'sequence') )
    c = sequence( c ); end;



%compute how often we shall iterate
if( isempty(iteration) ); %compute how many points the outcome will have    
    iteration = 0;
    val = numel( c );
    oovector = ordering2vector( oo, maxiteration );
    while( true )
        if( isempty(oovector) ); 
            break; end;
        if( numel(oovector)<=iteration ); 
            iteration = iteration-1; 
            break; end;
        iteration = iteration+1;
        val = val*abs( det(S{oovector(iteration),2}) )*dim + numel( S{oovector(iteration),1} );
        if( val>numpoint || iteration>maxiteration ); 
            iteration = iteration-1; 
            break; 
        else; 
            iteration = iteration+1; end; end; end;

parsem( varargin, 'test' );

if( c.ndims>dim );
    warning( 'blf:dimension', 'Starting sequence has larger dimension than subdivision scheme.' ); end;

oovector = ordering2vector( oo, iteration ); %change So to simple vector
if( size(oovector,2)<iteration );
    warning( 'blf:oolength', 'Length of ordering shorter than iteration.' );
    iteration = size( oovector, 2 ); end;

vprintf( 'ordering: %r\n', oovector, 'imp',[1 verbose] );

%do the work
PM = eye(dim);
for n = 1:iteration;
    if( issym(c.c) ); 
        vprintf( '|', 'imp',[1 verbose] ); end;
    if( n>iteration );
        break; end;
    if( ~isempty(maxnumpoint) && numel(c)>maxnumpoint ); 
        warning( 'blf:toomanypoint', 'Too many points. blf terminates early.' ); 
        break; end;
    if( isempty(maxnumpoint) && numel(c)>300000 )
        break; end;
    k = oovector( n );    
    if( k==0 ); 
        continue; end;

    c = conv( a{k}, upsample(c,M{k}) );
    if( removezero>=1 ); 
        if( numel(c)/nnz(c)>4 ); 
            c = c.shrink; 
            end; end; %remove zeros if there are a lot of zeros in c
    %vdisp(c); %DEBUG
    PM = PM*M{k};
    if( verbose>=4 );  
        plot( c.c ); 
        drawnow; 
        pause; end; end;

if( issym(c.c) ); 
    vprintf( '\n', 'imp',[1 verbose] ); end;
if( removezero>=1 ); 
    c = c.shrink; end;

%check if one should set <'removezero',0>
%if(nnz(c)<(dim+1)*numel(c) && removezero); vprintf('You maybe should set <''removezero'',0>.\n','imp',[1 verbose]); end;

%Postprocessing
if( any(diff) );
    c = diffsequence( c, diff ); end; 

if( ~iscell(c) ); 
    val = c; 
    c = {}; 
    c{1} = val; end; %make to cell array

%compute xyzv
if( nargout>=3 || ~isequal(plotval,0) );
    xyzv = cell( size(c) );
    for i = 1:numel( c );
        val = compresscoordinates( c{i}.c, PM, 'idx',c{i}.idx );
        xyzv{i} = val;
        if( removezero>=2 ); 
            xyzv{i}(:,xyzv{i}(end,:)==0) = []; end; end; end;

%Plot
if( ~isequal(plotval,0) );     
    switch numel( c ); %if c has more output values
        case {1,2,3}; 
            GRID = {1 numel( c )};
        case {4,5,6,7,8}; 
            GRID = {2 ceil( numel(c)/2 )};
        case {9}; 
            GRID = {3 3};
        otherwise; 
            GRID = {1 numel( c )}; end
    for i = 1:numel( c )
       if( numel(c)>1 );
            subplot( GRID{:}, i ); end;
       plotm( xyzv{i}, plotval{:} ); end; end;

if( numel(c)==1 ); 
    if( nargout>=3 || ~isequal(plotval,0) ); 
        xyzv = xyzv{1}; end;
    c = c{1}; end;

if( nargout>=1 ); 
    varargout{1} = c; end;
if( nargout>=2 ); 
    varargout{2} = PM; end;
if( nargout>=3 ); 
    varargout{3} = xyzv; end;
if( nargout>=4 ); 
    varargout{4} = oo; end;
    
end



function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.
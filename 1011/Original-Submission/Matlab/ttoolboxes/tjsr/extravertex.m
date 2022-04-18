function [ Vextra, VVVextra ] = extravertex( varargin );
% [ Ve, VVe ] = extravertex( V, [options] );
% [ Ve, VVe ] = extravertex( M, oo, v0, [options] );
% [ Ve, VVe ] = extravertex( V, M, oo, v0, [options] );
% Computes vertices such that a polytope augmented with those vertices has non-empty interior
%
% Input:
% ======
%   V                       matrix of column vectors, the polytope
%   M                       1xJ cell array of square matrices
%   oo                      column vector, OR
%                           1xN cell array of column vectors, orderings of matrix products
%   v0                      column vector, OR
%                           1xN cell array of column vectors, 
%                               Computes the polytope VV from M, oo, and v0, precisely
%                               given ordering oo{1} and vector v0{1}, the algorithm computes the vectors 
%                                   M{oo{1}(1)}*v0, M{oo{1}(1)}*M{oo{1}(2)}*v0, M{oo{1}(1)}*M{oo{1}(2)}*...*M{oo{1}(end)}*v0
%
% Options:
% ==========
%   'verbose',val           default=1, verbose level
%   'threshold',val         default=.01, influences the number of added vectors. The smaller the number, the less vectors are added.
%                           If threshold=0, VVe still is a basis for R^n. If threshold is big, VVe will have more vectors than its dimension
%   'scale',val             double, default=1, factor with which the matrices M are scaled before the vectors of the polytope VV are computed
%   'polytopetype',val      integer, default=1, The polytope which shall be considered.
%                               0 (TJSR_CONEFUNCT)     cone (used in cone functional)
%                               1 (TJSR_MINKFUNCT)     symmetric polytope (used for Minkowski Norm)
%                           Note: A method for complex polytope is not yet implemented
%   'celloutput'            Function returns the vectors as a row cell array of column vectors
%
% Output:
% =======
%   Ve      matrix of extra column vectors
%   VVe     [V Vextra]
%
% Written by: tommsch, 2020 

[verbose,varargin] = parsem( {'verbose','v'}, varargin, 1 ); 
[threshold,varargin] = parsem( {'threshold','t'}, varargin, 0.01 ); 
[scale,varargin] = parsem( {'scale','s'}, varargin, 1 ); 
[polytopetype,varargin] = parsem( {'polytopetype','pt'}, varargin, 1 ); 
[celloutput,varargin] = parsem( {'celloutput','c'}, varargin ); 

VV1=[];
VV2=[];
if( size(varargin,2)==1 || size(varargin,2)==4 )
    VV1 = varargin{1};
    varargin(1) = []; end;
if( size(varargin,2)==3 )
    M = varargin{1};
    oo = varargin{2};
    v0 = varargin{3};
    varargin(1:3)=[];
    if( ~iscell(oo) || ~iscell(v0) );
        oo = {oo}; 
        v0 = {v0}; end;
    %compute whole cycle, but I don't want to change the whole program.
    for i = 1:length( M )
        M{i} = M{i}.*scale; end;
    VV2 = cell( 1, length(oo) );

    for i = 1:length( oo )
        oo_temp = arrayfun( @(x) oo{i}(1:x),  1:length(oo{i}), 'UniformOutput',0 );
        val = cellfun( @(x) tbuildproduct(M, x)*v0{i}, oo_temp, 'UniformOutput',0 );
        VV2{i} = [val{:}]; end;
    VV2 = [VV2{:}]; end;
VV=[VV1 VV2];

    
parsem( varargin, 'test' );

dim = size(VV,1);
if( polytopetype==TJSR_CONEFUNCT )
    if( isempty(VV) )
        Vextra = eye( dim );
    else
        idx = max( VV, [], 2 )<threshold;
        Vextra = eye( dim );
        Vextra = Vextra(:,idx); end;

elseif( polytopetype==TJSR_MINKFUNCT );
    if( isempty(VV) )
        Vextra = eye( dim );
    else
        [U,SVD,~] = svd( VV );
        if( ~isvector(SVD) );
            SVD = diag( SVD ); end; %make singular values to a vector
        if( length(SVD)<dim ); 
            SVD(dim) = 0; end; %fill up with zeros
        rat = SVD/max( SVD ); 
        idx = rat<threshold;
        Vextra = normalizematrix(U(:,idx),'positive',1)*max( SVD ); end;

elseif( polytopetype==TJSR_COMPLEXFUNCT );
    vprintf( '''extravertex'' for complex case not yet implemented.\n', 'cpr','err', 'imp',[1 verbose] );
    dim = size( v0{1}, 1 );
    Vextra = zeros( dim, 0 );
else
    error( 'Unkown polytopetype given.' ); end;

VVVextra = [VV Vextra];

if( celloutput )
    Vextra = num2cell( Vextra, 1 ); 
    VVVextra = num2cell( VVVextra, 1 ); end;

    
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.
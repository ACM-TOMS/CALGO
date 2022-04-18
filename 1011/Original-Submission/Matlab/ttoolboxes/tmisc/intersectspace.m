function [ intersection, dim ] = intersectspace(varargin)
% [ intersection, dim ] = intersectspace( A, B, C, ... )
% Finds a basis of the intersection of subspaces
%
% Input:    
%   A, B, C, ...    N matrices/vectors - matrices/vectors defining basis of the subspace 
%                   Bases must be column vectors!
% 
% Output:
%   intersection    basis of the intersection space
%   dim             dimension of the intersection space
%
% E.g.: intersectspace([1 -1 0 0; 1 0 -2 1; 0 1 -2 -1]',[1 1 0 0; 1 0 -1 0]')
%
% Written by: Ondrej Sluciak <ondrej.sluciak@nt.tuwien.ac.at>
% Changed by: tommsch, 2017, to handle symbolic stuff and to use column vectors as bases

%#ok<*ALIGN>

N = nargin;

for i = 1:nargin;
    varargin{i}=varargin{i}.'; end;

if( N<1 );
    error( 'intersectspace:InvalidInput', 'Invalid input.' ); end;

if( N == 1 );
    intersection = rref( varargin{1} );
    intersection = removezero( intersection, 1 ).';
    dim = rank( intersection );
    return; end;

rnk1 = rank( varargin{1} );

for i = 2:N
    rnk1 = rnk1 + rank(varargin{i});
    if( size(varargin{i},2) ~= size(varargin{i-1},2) )
        error( 'intersectspace:InvalidInput', 'Sizes of the spaces must be equal.' ); end;
end

rnk2      = rank( [repmat(varargin{1},1,N-1); blkdiag(varargin{2:end})] ); % formula by Yongge Tian: www.math-cs.ucmo.edu/~mjms/2002.2/ytian2.ps
dim = rnk1 - rnk2;                                                 % dimension of the intersection

if( dim < 1 );
    intersection = zeros( 1, size(varargin{i},2) ).';
    return; end;

if( issym(varargin{1}) || issym(varargin{2}) );
    tmp = null( [varargin{1}',varargin{2}'] );
else
    tmp = null( [varargin{1}',varargin{2}'], 'r' ); end;
tmp = rref( tmp(end-size(varargin{2},1)+1:end,:)'*varargin{2} );

for i = 3:N
    tmp = null( [tmp',varargin{i}'], 'r' );
    tmp = rref( tmp(end-size(varargin{i},1)+1:end,:)'*varargin{i} ); end;

intersection = tmp(1:dim,:).';


end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.
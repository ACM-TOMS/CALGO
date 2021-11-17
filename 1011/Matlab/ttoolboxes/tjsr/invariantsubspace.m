function [ Mret, B ] = invariantsubspace(varargin)
% [ Mret, B ] = invariantsubspace( M, ['type'], [options] )
% Searches for invariant subspaces of matrices M.
%
% Input:
%   M       Cell array of matrices
%   'type'  A string of the following
%               'none'      Nothing happens,  Mret={M}, B=eye(dim)
%               'auto'      (default)
%               'perm'      Tries to find a permutation such that M is in block-diagonal form, using tpermTriangul, originally from the JSR-toolbox by Jungers
%               'basis'     Tries to find a basis such that M is in block-diagonal form, using tjointTriangul, originally from the JSR-toolbox by Jungers
%               'trans'     If M are transition matrices, tries to find subspaces of differences U, as described in <'Smoothness of anisotropic wavelets, frames and subdivision schemes', M. Charina and V.Yu. Protasov> 
%                           First tries numerically, then symbolically, then with vpa.
%   If type is not given, then the algorithm tries in the following order: 'perm', 'basis', 'trans' (only numerically)
%
% Options:
%   'verbose',val           Verbose level, Default=0
%   'double'                Casts matrices to double before computing subspaces
%   'vpa                    Casts matrices to vpa before computing subspaces
%   'sym'                   Casts matrices to sym before computing subspaces
%   'epsilon',val           Epsilon used in algorithm 'basis'
%   'V0',val                Subspace used in algorithm 'basis'
%
% Output:
%   Mret    cell array of the blocks in the block diagonal. 
%   B       The change of basis. I.e. B'*M{i}*B has the block diagonal form
%
% E.g.: [M,B]=invariantsubspace({[1 1 1; 1 0 1; 0 1 0], [1 -1 0; 0 2 0; 1 1 2]},'verbose',2)
%
% See also: tjointTriangul, tpermTriangul, constructU, restrictmatrix
%
% Written by: tommsch, 2018, and Jungers (JSR-toolbox)

%Rename: 'trans' to 'diff'

%XX This function is written very ugly but is working well at the moment
%XX Returns Error: invariantsubspace({[1],[2]},'basis')

[V0,varargin]         = parsem({'V0'}, varargin, [] );
[epsilon,varargin]    = parsem({'epsilon','eps','e'}, varargin, 1e-12 );
[verbose,varargin]    = parsem({'verbose','v'}, varargin, 1 );
[doubleflag,varargin] = parsem({'double','d'}, varargin );
[vpaflag,varargin]    = parsem({'vpa'}, varargin );
[symflag,varargin]    = parsem({'sym'}, varargin );
M = varargin{1};
varargin(1) = [];
dim = size( M{1}, 2 );
J = numel( M );
if( isempty(varargin) ); 
    type = 'auto'; 
else; 
    type = varargin{1}; 
    varargin(1) = []; end;


parsem( varargin, 'test' );

if( doubleflag )
    for j=1:J;
        try;
            M{j}=double(M{j});
        catch;
            vprintf('Cannot convert matrix %i to double.\n', j, 'imp',[0 verbose], 'cpr','err' ); end; end; end;

if( vpaflag )
    for j = 1:J;
        try;
            M{j} = vpa( M{j} );
        catch;
            vprintf( 'Cannot convert matrix %i to vpa.\n', j, 'imp',[0 verbose], 'cpr','err' ); end; end; end;

if( symflag )
    for j = 1:J;
        M{j} = sym( M{j} ); end; end;


info = identifymatrix( M );
if( isequal(type,'auto') && info.nonneg ); 
    vprintf( 'Since the matrices are nonnegative, we only search permutations which triangulises. \n', 'imp',[2,verbose] );
    type='perm'; end;

EYE = eye(dim);
vprintf( 'Search for common invariant subspaces: ', 'imp',[2 verbose] )
triang = false;
switch type
    case 'none'
        B=EYE;
        Mret={M};    
        
    case {'perm','permutation'}
        searchperm();
        
    case 'basis'
        searchbasis();
        
    case {'trans','transition','transitionmatrix'}
        searchtrans(1); %flag says: symbolically
        
    otherwise;
        searchperm();
        if(~triang); 
            searchbasis(); end
        if(~triang); 
            s = warning( 'off', 'restrictmatrix:notinvariant' );
            searchtrans( 0 ); 
            warning( s );
            end; end; %flag says: not symbolically
            

%vprintf( '\n', 'imp',[1 verbose] )
vprintf( '\nSizes of matrices: %r\n', cellfun(@(x) size(x{1},2),Mret), 'imp',[2 verbose] );


%functions to be called

    function searchperm;
        vprintf( 'Search for permutation: ', 'imp',[2 verbose] )
        [triang,Mret,B]=tpermTriangul(M);
        if(~triang);
            vprintf('no permutation found. ', 'imp',[2 verbose] )
            B = EYE;
            Mret = {M};
            triang = false;
        else
            vprintf( 'permutation found. ', 'imp',[2 verbose] )
            B  = EYE(B,:);
            triang = true; end; 
    end

    function searchbasis;
        vprintf( 'Search for basis: ', 'imp',[2 verbose] )
        [~,Mret,B]=tjointTriangul( M, V0, epsilon, 'verbose',verbose );
        if( ~isequal(B,eye(dim)) );
            vprintf( 'basis found. ', 'imp',[2 verbose] )
            triang = true;
        else
            vprintf( 'no basis found. ', 'imp',[2 verbose] )
            triang = false; end; 
    end

    function searchtrans( symflag );
        try
            vprintf( 'Search for difference subspace ', 'imp',[2 verbose] );
            if( symflag ); 
                vprintf('(symbolically): ', 'imp',[2 verbose] ); 
            else; 
                vprintf( '(numerically): ', 'imp',[2 verbose] ); end;
            U = constructU( M, 1, 'verbose',verbose-1 );
            val = identifymatrix( M ); 
            if( val.real ); 
                U = real( U ); end;
            [MU,~,MR,NN,BASIS] = restrictmatrix( M, U, 'verbose',verbose-2 );
            if( symflag )
                if( max(max([NN{:}]))>10e-12 ); 
                    U = constructU( M, 1, 'sym', 'verbose',verbose-2 );
                    [MU,~,MR,NN,BASIS] = restrictmatrix( M, U, 'verbose',verbose-2 ); end;
                if(max(max([NN{:}]))>10e-12); 
                    U = constructU( M, 1, 'vpa',200 );
                    [MU,~,MR,NN,BASIS] = restrictmatrix( M, U, 'verbose',verbose-2 ); end; end;
            if( max(max([NN{:}]))<10e-12 )
                vprintf( 'difference subspace found. \n', 'imp',[2 verbose] );
                Mret{1} = MU; 
                Mret{2} = MR; 
                B = BASIS;
                triang = true;
            else
                vprintf( 'no difference subspace found. \n', 'imp',[2 verbose] );
                B = EYE;
                Mret = {M};
                triang = false; end;
        catch
            vprintf( 'could not search for difference subspace. \n', 'imp',[2 verbose] );
            B = EYE;
            Mret = {M};
            triang = false; end;
    end
        
end
    
function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   

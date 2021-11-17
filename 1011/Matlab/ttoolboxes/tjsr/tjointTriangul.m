function [ triang, block, B ] = tjointTriangul( varargin )
% Heuristic for joint block-triangularization. There is no guarantee of finding every jointly invariant subspace.
% [ triang, blocks, B ] = tjointTriangul( M, [epsilon], [options] )
%     If the set M can be jointly block-triangularized then this function 
%     tries to find and return the sets of diagonal blocks. For this it tries     
%     to generate an invariant subspace by applying the matrices in M to each
%     eigenvector of every matrix of M.   
% [triang, blocks, B] = tjointTriangul( M, V0, [epsilon] )
%     Does the same but trying to 
%
% Input:
%   M               Matrices to be searched for common invariant subspaces
%   V0              If given, expands the subspace spanned by the linearly independent columns of V0 
%                   If not given or empty, expands the subspaces spanned by the eigenvectors of the matrices in M.
%   epsilon         Epsilon vor floating point arithmetic to compute ranks and to decided whether a value is 0.
%
% Options: 
%   'verbose',val   Verbose level
%
% Output:
%   triang          boolean indicating if the matrices could be jointly triangularized. 
%   B               unitary matrix giving the change of basis
%   block           the blocks of the common invariant subspaces in basis B
%                  
%                  blocks{1}{i}    *          *       *    . . .    * 
%                     0    
%                              blocks{2}{i}   *       *             .
%                     0            0                                .
%                     .            .       blocks{3}{i}             .
%      B'*M{i}*B =    .            .          0         .           
%                     .            .                        .  
%                                                              .    *
%                     0            0          0     . . .       blocks{q}{i}  
%
%
% REFERENCES
%    R.Jungers, 
%      "The Joint Spectral Radius: Theory and Applications" 
%      Vol. 385 section 1.2.2.5 in Lecture Notes in Control and Information
%      Sciences, Springer-Verlag. Berlin Heidelberg, June 2009

% Changed by: tommsch, 2019_10_18, Changed help and added input variable "epsilon"

% THIS ALGORITHM COULD BE OPTIMIZED
% By putting a flag on the subspaces that have no chance of being re-divided for instance. (Tree view)

[verbose,varargin] = parsem( {'verbose','v'}, varargin, 1 );

M = varargin{1};

n = size(M{1},1);
m = length(M);

if( n==1 )
    triang = 0;
    block = M;
    B = 1;
    return; end;

if( isscalar(varargin{size(varargin,2)}) )
    epsilon = varargin{size(varargin,2)};
    varargin(size( varargin, 2 )) = [];
else
    epsilon = 1e-12; end;

if( size(varargin,2)==2 && ~isempty(varargin{2}) )
    V0 = varargin{2};    
    [triang, block, B] = subJointTriangul( M, V0, epsilon, [], [], verbose );
else    
    triang = 0;
    block = {M};
    B = eye(n);
    imat = 0;
    
    while( ~triang && imat<m )
        imat = imat+1;
        if (issparse(M{imat}) && n>50)
            [V, ~] = eigs(M{imat},n);
        else
            [V, ~] = eig(full(M{imat})); end;

        ieigvec = 0;
        if( issym(M{1}) );
            vprintf('\nNumber of Loops: %i', n, 'imp',[1 verbose] ); end;
        while( ~triang && ieigvec<n )
            if( issym(M{1}) );
                vprintf( '\nLoop %i: ', ieigvec+1, 'imp',[1 verbose] ); end;
            ieigvec = ieigvec+1;
            [triang, blo, Bas] = subJointTriangul( M, V(:,ieigvec), epsilon, [], [], verbose );
            %fprintf( '\nComputing for eigvec %d of matrix %d', ieigvec, imat ) 

            if (triang)
                block = blo;
                B = Bas; end; end; end; end;
end

function [ triang, blocks, B ] = subJointTriangul( M, V0, epsilon, step, indnew, verbose )
n = size( M{1}, 1 );
m = length( M );

if( n==1 )
    triang = 0;
    blocks = {M};
    B = 1;
    return; end;

if( isempty(step) && isempty(indnew) )
    step = 1;
    indnew = 1; end;

if(issym(M{1}))
    B = sym( zeros(n,n) );
else
    B = zeros( n, n ); end;

nv = rank( V0 ); % Number of linearly independent vectors currently found

[nvnew, Vvecnew] = appMat( M, V0, epsilon, nv, indnew, verbose );

if( issym(M{1}) ); 
    condval = 1e6; 
else; 
    condval = 1e3; end;

if( nvnew==0 && cond(V0)<condval )
    [Qr,~] = qr(V0);
    Q1 = Qr(:,1:nv);
    Q2 = Qr(:,nv+1:n);

    Mtri1 = cell(1,m);
    Mtri2 = cell(1,m);

    for imat = 1:m
        Mtri1{imat} = Q1'*M{imat}*Q1;
        Mtri2{imat} = Q2'*M{imat}*Q2; end;

    [sub_triang1, sub_blocs1, sub_B1] = tjointTriangul( Mtri1, 'verbose',verbose );
    [sub_triang2, sub_blocs2, sub_B2] = tjointTriangul( Mtri2, 'verbose',verbose );

    if( sub_triang1 )
        blocks = sub_blocs1;
        B(:,1:nv) = Q1*sub_B1;
    else
        blocks = {Mtri1};
        B(:,1:nv) = Q1; end;

    if( sub_triang2 )
        blocks = [blocks, sub_blocs2];
        B(:,nv+1:n) = Q2*sub_B2;
    else 
        blocks = [blocks, {Mtri2}];
        B(:,nv+1:n) = Q2; end;

    Q.Q1 = Q1;
    Q.Q2 = Q2; %#ok<STRNU>
    triang=1;

elseif( nv+nvnew == n || isAlways(cond(V0)>=1e3,'Unknown','false') )
    if( step==1 )
        if( issym(V0) );
            V0 = sym(zeros(n,n)); 
        else
            V0 = zeros( n, n ); end;
        V0(:,1) = Vvecnew(:,end);
        [triang, blocks, B] = subJointTriangul( M, V0, epsilon, 2, 1, verbose );
    else % Step == 2, we have generated the whole space from both ends
        triang = 0;
        blocks = M;
        B = eye( n ); end;
else
    V0 = Vvecnew;
    [triang, blocks, B] = subJointTriangul( M, V0, epsilon, step, nv+1, verbose ); end;
end


function [ nvnew, Vvecnew ] = appMat( M, V0, epsilon, nv, indNew, verbose )
m = length( M );
Vvecnew = V0(:,1:nv);
nvnew = 0;
if(issym(M{1}));
    vprintf( '(%i)',m*(nv-indNew+1), 'imp',[1 verbose] ); end
for imat = 1:m
    for ivec = indNew:nv
        vectemp = M{imat}*V0(:,ivec);
        normvec = norm( vectemp );
        if( issym(vectemp) );
                vprintf( '.', 'imp',[1 verbose] ); end
        if( isAlways(normvec > epsilon,'Unknown','true') ) %default epsilon=1e-12;
            Vvectemp = [Vvecnew(:,1:(nv+nvnew)) vectemp/normvec];
            if( issym(Vvectemp) );
                %rankTemp = rank(vpa(Vvectemp));
                rankTemp = rank(Vvectemp);
            else
                rankTemp = rank(Vvectemp,epsilon); end;
            if( rankTemp == nv+nvnew+1 )
                Vvecnew = Vvectemp;
                nvnew = nvnew+1; end; end; end; end;
end


function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 

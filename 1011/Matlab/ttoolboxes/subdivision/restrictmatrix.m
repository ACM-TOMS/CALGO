function [TA, TT, TR, NULL, BASIS, flag] = restrictmatrix(varargin)
% [ TA, TT, TR, NULL, BASIS ] = restrictmatrix( T, A, [options] )
% Restricts matrices T to a subspace A.
% Also tests, if the matrices are invariant on A.
%
% Input: 
%   T                   square matrix or cell array of square matrices
%   A                   rectangular matrix giving the columns of the subspace A, to which T will be restricted
%
% Options:
%   'verbose',val       default: 1, verbose level
%   'epsilon',val       default: 1e-11, controls the offset to test whether A is invariant or not
%   'smallsize',bool    default: true, Removes unnecessary columns of A (starting from the last column) prior to computation. Returns rank(A) x rank(A) matrix. 
%   'basis',val         default: [], column vectors which complement A to a full basis. If not given, this is constructed using unit vectors
% 
% Output:
%   TA              dimA x dimA matrices. Restriction of T to A
%   TT              dimT x dimT matrices. T in the basis of A complemented to a basis of RR^dimT
%   TR              dim(T-A) x dim(T-A) matrices. Lower right corner of matrix TT
%   NULL            Lower left corner of T. If T is A invariant, then NULL consists only of zeros.
%   BASIS           the complemented basis of A. Thus BASIS\T*BASIS == T
%   flag            Is set to 1 if subspace A is invariant under all Matrices. 0 if not. NaN if unkown.
%   
% E.g.: vdisp(restrictmatrix([1 1 0; 0 1 1; 1 0 1],[1 -1 0; 0 1 -1]'))
%
% Written by: tommsch, 2018

% Changelog: tommsch, 2019_10_25,   Added output variable 'flag'
%            tommsch, 2019_01_19,   Added option 'basis'


if( ~iscell(varargin{1}) ); %make T a cell array
    cellflag = 0; 
    TT{1} = varargin{1}; 
else; 
    cellflag = 1; 
    TT = varargin{1}; end; 
A = varargin{2}; %subspace
J = size(TT,2); %number of matrices T
dim = size(TT{1},1);
epsilon = parsem( {'epsilon','eps','e'}, varargin, 1e-11 );
verbose = parsem( {'verbose','v'}, varargin, 1 );
smallsize = parsem( {'smallsize','ss','s'}, varargin, 1 );
basis = parsem( {'basis','b'}, varargin, [] );

RANK = rank( A );
CORANK = dim-RANK;

if( smallsize ) %remove linear dependent vectors from A
    for j = size(A,2):-1:1;
        if( rank(A(:,[1:j-1 j+1:end]) )==RANK );  
            A(:,j) = 0; end; end;
    A = removezero( A, 2 ); end; %remove the zero columns

SZE = size( A, 2 ); %number of columns of A
if( isempty(basis) )
    BASIS = A;
    %make A full rank
    BASIS = [BASIS eye(dim)]; %add a base of \RR^dim 
    for j = size(BASIS,2) : -1:SZE+1; %size(A,2)-dim+1:size(A,2);
        if( rank(BASIS(:,[1:j-1 j+1:end]))==dim );  %remove linear dependent vectors.
            BASIS(:,j) = 0; end; end;
    ADD = removezero( BASIS(:,SZE+1:end), 2 ); %remove the zero columns from the unit base
    BASIS = [BASIS(:,1:SZE) ADD];
    clear ADD;
else;
    BASIS = [A basis];
    if( rank(BASIS)<dim );
        error( 'restrictmatrix:nobasis', 'Given ''basis'' does not constitute a complement to ''A''.\n' ); end; end;

if( rank(BASIS(:,1:SZE))~=RANK || rank(BASIS)~=dim );
    error( 'restrictmatrix:failure', 'restrictmatrix: Wrong computation.\n' ); end;

[TA,TR,NULL] = deal( cell(size(TT)) );
MAXERR = 0;
for j = 1:J %go through all matrices in T
    TTval = BASIS\TT{j}*BASIS;
    NULL{j} = TTval(end-CORANK+1:end,1:SZE); %the residuum in the lower left corner. Equals zero if the matrix is invariant
    TR{j} = TTval(end-CORANK+1:end,SZE+1:end); 
    TT{j} = TTval; %#ok<AGROW>
    TA{j} = TTval(1:SZE,1:SZE); end;

flag = 1;
try;
    errorprinted = 0;
    for j = 1:J
        ERR = norm( double(simplify(NULL{j})) );
        if( ERR>epsilon*dim );
            if( errorprinted==0 ); 
                warning( 'restrictmatrix:notinvariant', 'Subspace not invariant under matrices.' ); end;
            vprintf( '%i (err = %g >> 0),  ',j,sum(abs(NULL{j}(:))), 'cpr','err', 'imp',[1 verbose] );
            errorprinted = 1;
            flag = 0; end; end ;
    if( errorprinted==1 ); 
        vprintf( '\n', 'imp',[1 verbose] ); end;
catch
    warning( 'restrictmatrix:notest', 'Could not test invariance under matrices. It is likely that the matrices are not invariant. Be careful and check output NULL (4th output argument).');
    flag = NaN; end;
vprintf( 'Maximal error during computation = %g \n', MAXERR, 'imp',[2 verbose] );

if( ~cellflag )
    TA = TA{1};
    TT = TT{1};
    TR = TR{1};
    NULL = NULL{1}; end;
    
end

function dummy; end %#ok<DEFNU>  %Generates an error, if the 'end' of a function is missing.
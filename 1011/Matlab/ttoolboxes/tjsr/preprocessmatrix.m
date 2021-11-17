function [M] = preprocessmatrix(varargin)
% [ M ] = preprocessmatrix( M, [options] )
% Simplifies sets of matrices while preserving its joint spectral radius.
% For default values: the output set has the same JSR as the input set.
%
% Input:
%   M           (cell-array) matrices to be preprocessed
%
% Options: 
%   'sym',bool               default: false, Transforms matrices to double
%   'double',bool            default: false, Transforms matrices to double
%   'inverse',bool           default: false, Takes the Moore-penrose pseudo-inverse of all matrices. 
%   'addinverse',bool        default: false, Adds the transposed matrices. 
%   'transpose',bool         default: false, Transposes all matrices. 
%   'addtranspose',bool      default: false, Adds the transposed matrices. 
%   'makepositive',bool      default: true,  Makes first non-zero entry of array positiv.
%   'timestep',val           default: false, Changes matrices to I*(1-val) + M*val. 
%   'perturbate',val         default: 0,     Perturbates matrices randomly by rand(val). 
%   'removezero',bool        default: true, Removes all (but one) zero matrices. 
%   'removeduplicate,bool    default: true, Removes duplicates. Has no tolerance in comparisons. 
%   'basechange,val          default: 0,    Performs a base change.
%                               val = 0          No basechange.
%                               val = 'random'   a random base change is made.
%                               val = 1          The singular value decomposition of M{1} is used. Thus M{1} is in Jordan-normal-form. If M{1} is scaled baldy, M{2} is used, etc... .
%                               val = matrix     The matrix is used. B must be invertible. M{i}=B^(-1)*M{i}*B
%   'exponential',val        default: false, The matrix exponential of each matrix is taken.
%   'nodouble',val           default: false, Does not convert matrices to double. 
%   'verbose',val            default: 1, Verbose level.
%
% If no options are given, the matrices are processed in the order as written above. The second to last step is again 'makepositive'
% If at least one option (except 'verbose') is given, only that option is used.
%
% Output:
%   M           cell-array of preprocessed matrices
%               If no matrices are removed or added, the return array has the same size as the input array.
%               Otherwise it is a vector.
%
% E.g.: vdisp( preprocessmatrix({[-1 2; 2 3],[-1 2; 2 3]}) )
%
% Written by tommsch, 2018

%#ok<*ALIGN>

% XX Add options of: ON THE FINITENESS PROPERTY FOR RATIONAL MATRICES, RAPHAËLL JUNGERS AND VINCENT D. BLONDEL
% XX 'removeduplicate,val    Removes duplicates with tolerance val. Default val=0.

[verbose,varargin] = parsem( {'verbose','v'}, varargin, 1 );
if( size(varargin,2)==1 ) %if no options are given, we use default values
    symflag = 0;
    doubleflag = 0;
    inverse = 0;
    addinverse = 0;
    transpose = 0;
    addtranspose = 0;
    makepositiveflag = 1;
    timestep = 1;
    perturbate = 0;
    removezero = 1;
    removeduplicate = 1;
    basechange = 0;
    exponential = 0;
    nodouble = 0;
else  %if at least one option is given, only that option is used
    [symflag,varargin] = parsem( 'sym', varargin, 0 );
    [doubleflag,varargin] = parsem( 'inverse', varargin, 0 );
    [inverse,varargin] = parsem( 'inverse', varargin, 0 );
    [addinverse,varargin] = parsem( 'addinverse', varargin, 0 );
    [transpose,varargin] = parsem( 'transpose', varargin, 0 );
    [addtranspose,varargin] = parsem( 'addtranspose', varargin, 0 );
    [makepositiveflag,varargin] = parsem( 'makepositive', varargin, 0 );
    [timestep,varargin] = parsem( 'timestep', varargin, 1 );
    [perturbate,varargin] = parsem( 'perturbate',varargin, 0 );
    [removezero,varargin] = parsem( 'removezero', varargin, 0 );
    [removeduplicate,varargin] = parsem( 'removeduplicate',varargin, 0 );
    [basechange,varargin] = parsem( 'basechange', varargin, 0 );
    [exponential,varargin] = parsem( 'exponential', varargin, 0 );
    [nodouble,varargin] = parsem( 'nodouble', varargin, 0 );
end
M = varargin{1};
varargin(1) = [];
parsem( varargin, 'test' );

if( symflag );
    vprintf( 'Cast to sym. ',  'imp',[1 verbose] );
    J = numel( M );
    for j = 1:J;
        M{j} = sym( M{j} ); end; end;

if( doubleflag );
    vprintf( 'Cast to double. ',  'imp',[1 verbose] );
    J = numel( M );
    for j = 1:J;
        M{j} = double( M{j} ); end; end;

if( inverse )
    vprintf( 'Invert matrices. ', 'imp',[1 verbose] );
    J = numel( M );
    for j = 1:J;
        M{j} = pinv( M{j} ); end; end;

if( addinverse )
    J = numel( M );
    dim = size( M, 1 );
    for j = 1:J; 
        if( rank(M{j})<dim ); 
            vprintf( 'Matrix %i not invertible\n',j,'cpr','err', 'imp',[0 verbose]); end;
        M{j+J}=inv(M{j}); end; end;

if( transpose );
    J = numel( M );
    vprintf( 'Transpose matrices. ' , 'imp',[1 verbose] );
    for j = 1:J; 
        M{j}=M{j}.'; end; end;

if( addtranspose ); 
    J = numel( M );
    M = M(:);
    vprintf( 'Add transposed matrices. ', 'imp',[1 verbose] );
    M{2*J}=[]; %preallocate space
    for j = 1:J; 
        M{j+J}=M{j}'; end; end;

if( makepositiveflag )
    J = numel( M );
    vprintf( 'Make matrices positive. ', 'imp',[1 verbose] );
    for i = 1:J; 
        M{i}=makepositive(M{i});  end; end;

if( timestep~=1 )
    J = numel( M );
    EYE = eye( size(M{1},1) );
    for j = 1:J; 
        M{j} = (timestep-1)*EYE+timestep*M{j}; end; end;

if( perturbate~=0 )
    J = numel( M );
    vprintf( 'Perturbate matrices. ', 'imp',[1 verbose] );
    for i=1:J; 
        dim=size(M{i},2);
        M{i}=M{i}+randn(dim)*perturbate; end; 
end
if( removezero )
    J = numel( M );
    vprintf( 'Remove zero matrices. ', 'imp',[1 verbose] );
    for j = 2:J;
        if( ~any(M{j}) ); 
            M{j} = []; end; end;
    idx = cellfun( @isempty, M );
    if( anym(idx) ); 
        M = M( ~idx ); end; end;

if( ~nodouble )
    J = numel( M );
    for j = 1:J
        M{j} = double( M{j} ); end; end;

if( removeduplicate )
    vprintf( 'Remove duplicates. ', 'imp',[1 verbose] );
    M = uniquecell( M, 'stable' ); end;

if( basechange )
    J = numel( M );
    dim = size( M{1}, 1 );
    if( isequal(basechange,'random') ); %generate base-change matrix
        while( true )
            basechange = randn(dim);
            if( cond(basechange)<10*dim );
                break; end; end; end;
    if( isequal(basechange,1) )
        i = 0;
        while( true )
            i = i+1;
            [V,D] = eig( M{i} );
            [V,~] = cdf2rdf( V, D );
            if( cond(V)<100*size(M{1},2) ); 
                vprintf( 'Make basechange. ', 'imp',[1 verbose] );
                for j = 1:J
                    M{j} = V\M{j}*V; end;
                break; end
            if( i>=J ); 
                break; end; end;
    elseif( issquare(basechange) )
        for j = 1:J
            M{j} = basechange\M{j}*basechange; end;
    else
        vprintf( 'Wrong value for ''basechange''.', 'imp',[1 verbose] ); end; end;

if( exponential )
    J = numel( M );
    for j = 1:J; 
        M{j} = expm( M{j} ); end; end;

if( makepositiveflag )
    J = numel( M );
    vprintf( 'Make matrices positive again. ', 'imp',[1 verbose] );
    for i = 1:J; 
        M{i} = makepositive( M{i} ); end; end;

vprintf( '\n', 'imp',[1 verbose] );

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   
function [type] = identifymatrix( varargin ) 
% [type] = identifymatrix( M , [allflag], 'what' ) 
% Returns standard properties of matrices.
%
% Input: 
%   M           matrix, or cell array of matrices. Works only for non-sparse matrices
%               For symbolic matrices the output may not correct
%   allflag     (optional), integer, if set, more properties are tested
%   what        (optional), string, if set, only this property is returned
%               Note that, although only one value is returned in this case, all other values are still computed
%   
% Output: 
%   type    a struct array which holds a distinct set of properties of M.
%
%               type.nM             number of matrices
%               type.tensordim      the maximal dimension of all matrices (as given by ndimsm)
%               type.samesize       checks if all matrices have the same size
%               type.square         checks if all matrices are square-matrices
%               type.sym            checks if any matrix is symbolic
%               type.size           sizes of the matrices
%
%   If the matrices are square matrices, the following tests are also done
%               type.int            checks if all entries are integers (does check the value, not the type!)
%               type.bool           checks if all entries are boolean
%               type.pm1            checks if all entries are eiterh -1,0 or 1
%               type.nonneg         checks if all entries are nonnegative
%               type.pos            checks if all entries are positive
%               type.finite         checks if all entries are finite (excluding inf and nan)
%               type.real           checks if all entries are real
%               type.symmetric      checks if all matrices are symmetric (M=M.')
%               type.hermitian      checks if all matrices are hermitian (M=M')
%               type.columnstochastic      if all matrices are column-stochastic
%               type.rowstochastic  checks if all matrices are column-stochastic
%               type.doublestochastic      if all matrices are column-stochastic
%               type.hessl    checks if all matrices are lower Hessenberg
%               type.hessu    checks if all matrices are upper Hessenberg
%               type.tril    checks if all matrices are lower triangular
%               type.triu    checks if all matrices are upper triangular
%
%   Properties which are only set if allflag>=1
%               type.inv            checks if all matrices are invertible
%               type.normal         checks if all matrices are normal (M*M'=M'*M) (See Jungers, The JSR, Prop 2.2)
%               type.possemidef     checks if all matrices are positive semi-definite
%               type.posdef         checks if all matrices are positive definite
%               type.unitary        checks if all matrices are unitary
%               type.zero           checks if all matrices have spectral radius zero
%               type.sparsity       rel. sparsity computed over all matrices
%
%   Properties which are only set if allflag>=2
%               type.commute        checks if all matrices commute with each other 
%
% E.g.: identifymatrix( {[1 2 3; 0 2 1; 0 0 1],[2]} )
%
% Written by: tommsch, 2019

%#ok<*TRYNC>

    %If M is not a cell array make it to one
    if( ~iscell(varargin{1}) );
        M = varargin(1); 
    else; 
        M = varargin{1}; end;
    varargin(1) = [];
    
    if( ~isempty(varargin) && ischar(varargin{size(varargin,2)}) ); 
        what = varargin{size( varargin, 2 )};
        varargin(size(varargin,2)) = [];
    else
        what=[]; end;
    
    if( size(varargin,2)>=1 ); 
        allflag = varargin{1};
        varargin(1) = []; %#ok<NASGU>
    else
        allflag = 0; end;
    


    type.nM = length(M);
    type.tensordim = max(cellfun(@(x) ndimsm(x), M));
    if( ~nnz(diff(cell2mat(cellfun(@(x) sizem(x,[],type.tensordim)', M, 'UniformOutput',false)),1,2)) );  
        type.samesize = 1;      
    else;
        type.samesize = 0; end; 
    
    type.square = all( cellfun(@(x) issquare(x), M) );
    type.real   = all( cellfun(@(x) allm(isreal(x)), M) );
    type.finite = all( cellfun(@(x) allm(isfinite(x)), M) );
    type.nonneg = type.real && all( cellfun(@(x) allm(isAlways(x>=0, 'Unknown','false')), M) );
    type.pos    = type.real && all( cellfun(@(x) allm(isAlways(x>0, 'Unknown','false')), M) );    
    type.int    = type.finite && all( cellfun(@(x) allm(iswholenumber(x)), M) );
    type.bool   = type.finite && type.real && all( cellfun( @(x) allm( isAlways(x==0,'Unknown','false') | isAlways(x==1,'Unknown','false') ), M) );
    type.pm1    = type.finite && type.real && all( cellfun( @(x) allm( isAlways(x==0,'Unknown','false') | isAlways(x==1,'Unknown','false') | isAlways(x==-1,'Unknown','false') ), M) );
    type.nan    = any( cellfun(@(x) anym(isnan(x)), M) );
    type.sym    = any( cellfun(@(x) issym(x), M) );    
    type.size   = cellfun( @size,M, 'UniformOutput',false );
    
    if( allflag>=1 )
        type.sparsity = sum( cellfun(@nnz,M) )/sum( cellfun(@numel,M) ); end;
 
    if( type.square );
        if( type.sym ); 
            epsilon=0;
        else; 
            epsilon=10e-12*max( cellfun(@numel,type.size) ); end;
        
        epsilondouble = eps*max( cellfun(@max,type.size) );
        type.columnstochastic = all( cellfun(@(x) isAlways(norm(sum(x,1)-ones(1,size(x,2)))<epsilon, 'Unknown','false'),M) ) && type.nonneg;
        type.rowstochastic    = all( cellfun(@(x) isAlways(norm(sum(x,2)-ones(size(x,1),1))<epsilon, 'Unknown','false'),M) ) && type.nonneg;
        type.doublestochastic = type.columnstochastic && type.rowstochastic;        
        type.triu = all( cellfun(@(x) isequal(triu(x,0),x), M) );
        type.tril = all( cellfun(@(x) isequal(tril(x,0),x), M) );        
        type.hessu = all( cellfun(@(x) isequal(triu(x,-1),x), M) );
        type.hessl = all( cellfun(@(x) isequal(tril(x,1),x), M) );
        if( allflag>=1 && type.sym || allflag>=0 &&~type.sym )
            type.symmetric = type.finite && all( cellfun(@(x) isAlways(norm(x-x.',1)<epsilon, 'Unknown','false'), M) );
            type.hermitian = type.finite && all( cellfun(@(x) isAlways(norm(x-x' ,1)<epsilon, 'Unknown','false'), M) ); end;
        
        if(allflag>=1)
            type.normal = type.finite && all( cellfun(@(x) isAlways(norm(x*x'-x'*x,1)<epsilon, 'Unknown','false'), M) );
            try; 
                type.possemidef = type.finite && all( cellfun(@(x) allm(eig(double(x))>=0), M) ); 
            catch; 
                type.possemidef = 0; end;
            try; 
                type.posdef = type.finite && all( cellfun(@(x) allm(eig(double(x))>=epsilondouble), M)); 
            catch; 
                type.posdef=0; end;
            type.unitary    = type.finite && all( cellfun(@(x) isAlways(norm(x*x'-eye(size(x)),1)<epsilon,'Unknown','false'),M) );
            type.inv        = type.finite && all( cellfun(@(x) rcond(double(x))>epsilondouble, M) );
            type.zero       = type.finite && all( cellfun(@(x) rho(double(x))<epsilondouble, M) ); end;
        
        if( allflag>=2 )
            comm = @(A,B) A*B-B*A;
            type.commute = all( cellfun(@(x) norm(comm(M{x(1)},M{x(2)}))<epsilon,num2cell(mixvector(1:type.nM,2),1)) ); end;
    end;
    
    if( ~isempty(what) )
        try
            type = type.(what);
        catch
            warning( 'identifymatrix:wrongname', 'Given value not available' );
            type = orderfields( type ); end;
    else
        type = orderfields( type ); end;

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 
function [Tnorm] = normalizematrix(varargin)
% Tnorm = normalizematrix(T, how, [options])
% Normalizes matrices in various ways.
% Input:    T           (matrix, cell array or nested cell array of matrices) the matrices to be normalized
%           how         string which tells how the matrices shall be normalized. Default: 'colsum'
%               'colsum'            All columns sum up to 1. Columns which sum up to zero, are unchanged.
%               'rowsum'            All rows sum up to 1. --"--
%               'dirsum',dir        All lines in direction dir sum up to 1. --"--
%
%               'colmax'            Biggest value in each column in absolute value is set to 1. zero columns are unchanged.
%               'rowmax'            Biggest value in each row in absolute value is set to 1. --"--
%               'dirmax',dir        Biggest value in each line in direction dir in absolute value is set to 1. --"--
%
%               'colnorm',p         p-norm of each column is set to 1. zero columns are unchanged. p is a number or a function handle
%               'rownorm',p         p-norm of each row is set to 1. --"-- . --"--
%               'dirnorm',[dir p]   p-norm of each line in direction dir in absolute value is set to 1. --"-- . --"-- .
%
%               'rho'               Matrix is scaled such that spectral radius==1.
%               'norm',p            Matrix is scaled such that p-norm==1. p is a number or a function handle
%               'binary'            All nonzero entries are set to one
%               'dotprod',v         (experimental) Matrix is normalized such that dot(M,v)==1, M and v must be vectors. If dot(M,v)==0, then M stays unchanged.
%               'positive',dir      'makepositive' is applied to each line in direction dir
%
% Output:varargin{1} 
%   Tnorm               the normalized matrix
% 
% E.g.: normalizematrix(randn(4),'colsum')
%       vdisp(normalizematrix({randn(2), randn(2)},'rho'))
%       normalizematrix(randn(3,2), 'dirnorm', [1, inf])
%
% See also: makepositive


 %#ok<*ALIGN>

T=varargin{1};
if(isempty(T)); 
    Tnorm=T; 
    return; end;

how=varargin{2};

switch how
    case 'colsum';  
        Tnorm=nestedcellfun(@(x) normalizematrix_dirsum(x,1),                            T,'UniformOutput',false); 
    case 'rowsum';  
        Tnorm=nestedcellfun(@(x) normalizematrix_dirsum(x,2),                            T,'UniformOutput',false); 
    case 'dirsum';  
        Tnorm=nestedcellfun(@(x) normalizematrix_dirsum(x,parsem('dirsum',varargin,[])), T,'UniformOutput',false); 
        
    case 'colmax';  
        Tnorm=nestedcellfun(@(x) normalizematrix_dirnorm(x,[1,inf]),                             T,'UniformOutput',false); 
    case 'rowmax';  
        Tnorm=nestedcellfun(@(x) normalizematrix_dirnorm(x,[2,inf]),                             T,'UniformOutput',false); 
    case 'dirmax';  
        Tnorm=nestedcellfun(@(x) normalizematrix_dirnorm(x,[parsem('dirmax',varargin,[]), inf]), T,'UniformOutput',false); 
        
    case 'colnorm'; 
        Tnorm=nestedcellfun(@(x) normalizematrix_dirnorm(x,[1 parsem('colnorm',varargin,[])]),   T,'UniformOutput',false); 
    case 'rownorm'; 
        Tnorm=nestedcellfun(@(x) normalizematrix_dirnorm(x,[2 parsem('rownorm',varargin,[])]),   T,'UniformOutput',false); 
    case 'dirnorm'; 
        Tnorm=nestedcellfun(@(x) normalizematrix_dirnorm(x,parsem('dirnorm',varargin,[])),       T,'UniformOutput',false); 
    
    case 'norm';    
        Tnorm=nestedcellfun(@(x) normalizematrix_norm(x,parsem('norm',varargin,[])),          T,'UniformOutput',false); 
    case 'rho';     
        Tnorm=nestedcellfun(@(x) normalizematrix_rho(x),                                      T,'UniformOutput',false);     
    case {'binary','bin'};  
        Tnorm=nestedcellfun(@(x) normalizematrix_binary(x),                                   T,'UniformOutput',false);     
    case {'dotprod','dot'}; 
        Tnorm=nestedcellfun(@(x,y) normalizematrix_dotprod(x,y),                              T,parsem('dotprod',varargin,[]),'UniformOutput',false);     
    case {'positive','pos'};
        Tnorm=nestedcellfun(@(x) normalizematrix_positive(x,parsem({'positive','pos'},varargin,[])),  T,'UniformOutput',false); 
    otherwise;
        error('Unkown option');
end;

end

function [T] = normalizematrix_dirsum(T,dir)
    dim=ndimsm(T);
    dirsum=sum(T,dir);
    dirsum(dirsum==0)=1;
    REPEAT=ones(1,dim);
    REPEAT(dir)=size(T,dir);
    dirsum=repmat(dirsum, REPEAT);
    T=T./dirsum;
end

function [T] = normalizematrix_dirnorm(T,arg)

    dir=arg(1);
    normval=arg(2);

    if(isa(normval, 'function_handle'));
        normhandle=normval;
    else
        normhandle=@(x) norm(x,normval);
    end;


    dim=ndims(T); 
    if(iscolumn(T)); 
        dim=1; 
    end;
    val=cellfun(normhandle,num2cell(T,dir));
    val(val==0)=1;
    REPEAT=ones(1,dim);
    REPEAT(dir)=size(T,dir);
    val=repmat(val, [REPEAT 1]); 
    T=T./val;
end

function [T] = normalizematrix_norm(T,normval)
    if(isa(normval, 'function_handle'));
        T=T/normval(T);
    else
        T=T/norm(T,normval);
    end;
end

function [T] = normalizematrix_rho(T)
    T=T/trho(T);
end

function [T] = normalizematrix_binary(T)
    T=double(logical(T));
end

function [T] = normalizematrix_dotprod(T,v)
    val=dot(T,v);
    T=T/val;
end

function [T] = normalizematrix_positive(T,dir)
    val=num2cell(T,dir);
    val=cellfun(@makepositive,val,'UniformOutput',false);
    if(issym(val{1}));
        T=cell2sym(val);
    else;
        T=cell2mat(val);
    end;
end

function dummy; end %#ok<DEFNU>  %Generates an error, if the 'end' of a function is missing.
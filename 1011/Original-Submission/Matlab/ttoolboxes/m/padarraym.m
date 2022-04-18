function X = padarraym(X, r, flag)
% [ X ] = padarraym( X, r, [options] )
% Consistent behaviour for padarray,  with regards to multi-dimensional applications.
% Works for arrays and cells and symbolics.
% 
% Input:
%   X               the cell-array to be padded
%   r               how much to pad
%                       if scalar: padds in all dimension
%                       if vector: only padds in given dimension
%
% Options:
%   'post'          Only padd at the beginning
%   'pre'           Only padd at the end
%   (Default)       Padd everywhere
%
% Note:
%   If X is empty, the output is empty
%
% Output:
%   X           the padded cell-array
%   
% E.g.: padarraym({1 [1 2]; sym(2) 'a'},1,'pre')
%
% Written by: Notlikethat (stackoverflow.com/users/3156750/notlikethat), 2014
% Modified by: tommsch, 2018

if(nargin==3 && isequal(flag,'post')); 
    post=true; 
else; 
    post=false; 
end;
if(nargin==3 && isequal(flag,'pre')); 
    pre=true; 
else; 
    pre=false; 
end;

dim=max(length(r),ndimsm(X));
if(isscalar(r)); 
    r=r*ones(1,dim); end;

if(all(r==0)); 
    return; end;
if(any(r<0))
    error('''r'' must be nonnegative.'); end;

if(isempty(X))
    return;
elseif(numel(X)==1 && numel(r)==1); %workaround for wrong behaviour in this case
    flag=1;
else
    flag=0;
end

if(post)
    sz = num2cell(sizem(X,[],dim) + r);
    if(iscell(X));
        X{sz{:}} = [];     % hooray for comma-separated lists!
    else;
        X(sz{:}) = 0;     % hooray for comma-separated lists!
    end;
elseif(pre)
    sz = num2cell(sizem(X,[],dim) + r);
    if(iscell(X));
        X{sz{:}} = [];     % hooray for comma-separated lists!
    else;
        X(sz{:}) = 0;     % hooray for comma-separated lists!
    end;
    X = circshift(X, r);
else %both
    sz = num2cell(sizem(X,[],dim) + 2*r);
    if(iscell(X));
        X{sz{:}} = [];     % hooray for comma-separated lists!
    else;
        X(sz{:}) = 0;     % hooray for comma-separated lists!
    end;
    X = circshift(X, r);
end

if(flag)
    X=X.';
end

end

function dummy; end %#ok<DEFNU>  %Generates an error, if the 'end' of a function is missing.
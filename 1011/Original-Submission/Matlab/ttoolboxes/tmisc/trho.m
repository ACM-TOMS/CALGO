function rhoo = trho(M)
% [ rhoo ] = trho( M )
% Computes the spectral radius of matrices.
% If M is a cell of matrices, ans is a row vector with the spectral radius of each M{i}
%
% Input:
%   M       matrix or cell array of matrices
%
% Output:
%   rhoo    spectral radius of matrix/matrices
%
% Note:
%   In this file is also a version which uses vpa to compute rho. But it is slow, thus the code is commented out.
%
% E.g.: trho([1 2; -1 2])
%
% Taken from: Jungers, JSR-toolbox
% Modified by: tommsch, 2017


if(iscell(M))
    m = length(M);
    rhoo = zeros(1,m);

    for i = 1:m
        if(isempty(M{i}));
            rhoo(i)=-inf;
        elseif issparse(M{i})
            opts.disp = 0;
            rhoo(i) = eigs(M{i},1,'LM',opts);
        else
            rhoo(i) = max(abs(eig(M{i})));
        end
    end
else
    if(isempty(M))
        rhoo=[];
    elseif (issparse(M))
        opts.disp = 0;
        rhoo = eigs(M,1,'LM',opts);
    else
        rhoo = max(abs(eig(M)));
    end
end
end


function rhoo = trho_vpa(M,vpabound,newdigits)
% rhoo = trho(M,[vpabound],[newdigits])
% 
% Computes the spectral radius of matrix M
% if M is a cell of matrices, rhoo is a row vector with the spectral radius of each M{i}
%
% Options:
%    vpabound       If vpabound is set, then the spectral radius is computed with vpa, if abs(rho-vpabound)<epsilon
%                   epsilon  = epsilon=min([1e-12, 10^-size(M,2), max(1,1/abs(rhoo(1)-rhoo(2)))*10*eps*size(M,2)]); 
%                   If this happens, then the return value is symbolic
%                   
%    newdigits      digits used in vpa computation (if necessary)
%                   if not given, then newdigits = max(17, size(M,2)).
%
% changed by: tommsch, 2018


if(nargin<3); newdigits=0; end;
if(nargin<2); vpabound=-inf; end;

if(iscell(M))
    m = numel(M);
    rhoo = zeros(1,m);
    for i=1:m;
        rhoo(i)=trho(M{i},vpabound,newdigits);
    end
else
    %compute spectral radius numerically
    if(numel(M)==1);
        rhoo=abs(M);
        return;
    elseif (issparse(M))
        opts.disp = 0;
        rhoo = sort(eigs(M,2,'LM',opts),'descend','ComparisonMethod','abs');
    elseif(~isnumeric(M))
        rhoo=trhovpa(M, newdigits);
        vpabound=-inf;
    else
        rhoo=sort(eig(M),'descend','ComparisonMethod','abs');
    end
    
    if(vpabound<0); rhoo=rhoo(1); return; end;
    
    %test if we need vpa
    epsilon=min([1e-12, 10^-size(M,2), max(1,1/abs(rhoo(1)-rhoo(2)))*10*eps*size(M,2)]); 
    
    if(abs(rhoo(1)-vpabound)<epsilon);
        rhoo=trhovpa(M, newdigits);
    else
        rhoo=abs(rhoo(1));
    end
    rhoo=sym(rhoo);
end

end

function rhoo = trhovpa(M, newdigits)
    digitsbackup=digits;
    if(newdigits==0); newdigits=max(17,size(M,2)); end;
    digits(newdigits);
    rhoo=max(abs(eig(vpa(M))));
    digits(digitsbackup);
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 
function [varargout] = blkdiag(varargin)
%BLKDIAG (overloaded)

resul = [];
%Get the dimensions
for cont=1:nargin
    now = varargin{cont};
    if(isa(now,'rolmipvar'))
        now = now.data(1).value;
    end
    
    [rows(cont), cols(cont)] = size(now);
end

%Concatenate the matrices
for cont = 1:nargin
    aux = [];
    if (cont > 1)
        aux = zeros(rows(cont),sum(cols(1:cont-1)));
    end
    aux = [aux varargin{cont}];
    if (cont < length(cols))
        aux = [aux zeros(rows(cont),sum(cols(cont+1:end)))];
    end
    resul = [resul; aux];
end
varargout{1} = resul;
    
    
return
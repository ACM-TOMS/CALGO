function [varargout] = evalpar(varargin)
% Returns the value of the polynomial variable when assigning the given 
% values to the parameters.
%
% [var] = evalpar(poly,paramval) assigns the values on the cell array 
% PARAMVAL to the parameters of the input polynomial POLY, returning the
% result. For instance, if the variable is defined over 3 simplexes with
% vertices [2 4 3], then one right syntax would be
% [var] = evalpar(poly,{[0.5 0.5],[0.25 0.25 0.25 0.25],[0.4 0.4 0.2]})
% If PARAMVAL is a vector, instead of a cell array, then the polynomial is
% considered to be constructed from a matrix polynomially dependent on the
% parameters, being such parameters informed on PARAMVAL. 


if (nargin == 2)
    poly = varargin{1};
    paramval = varargin{2};
else
    error('Input error. Type ''help evalpar'' for more details');
end

if (length(paramval) ~= length(poly.data(1).exponent))
    error('The number of simplexes of the polynomial is different from the number of input parameters set');
end

if (iscell(paramval)) %Simplex parameters
    for cont=1:length(poly.data(1).exponent)
        if (length(poly.data(1).exponent{cont}) ~= length(paramval{cont}))
            if (sum(poly.data(1).exponent{cont}) == 0)
                for contexp = 1:length(poly.data)
                    poly.data(contexp).exponent{cont} = zeros(length(paramval{cont}));
                end
            else
                error('The numbers of vertices are inconsistent');
            end
        end
    end
else %Polynomially dependent parameters
    if (isempty(poly.bounds))
        error('The second argument must be a cell array if the input polynomial was not defined from a polynomially dependent set of bounded parameters');
    end
    for cont=1:length(paramval)
        if (length(poly.data(1).exponent{cont}) ~= 2)
            error('The second argument must be a cell array if the input polynomial was not defined from a polynomially dependent set of bounded parameters');
        end
    end
end

var = 0;
for cont=1:length(poly.data)
    param = 1;
    for contsimplex=1:length(poly.data(cont).exponent)
        if iscell(paramval)
            for contexp=1:length(poly.data(cont).exponent{contsimplex})
                param = param*paramval{contsimplex}(contexp)^(poly.data(cont).exponent{contsimplex}(contexp));
            end
        else
            alfa = (paramval(contsimplex) - poly.bounds(contsimplex,2))/(poly.bounds(contsimplex,1) - poly.bounds(contsimplex,2));
            param = param*alfa^(poly.data(cont).exponent{contsimplex}(1))*(1-alfa)^(poly.data(cont).exponent{contsimplex}(2));
        end
    end
    var = var + param*poly.data(cont).value;
end
varargout{1} = var;
   
return

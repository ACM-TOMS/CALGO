function tr = trace(varargin)
%TRACE  Sum of diagonal elements.
%   trace(poly) is the sum of the diagonal elements of poly. If the polynomial
%   is parameter-dependent, then its trace is also a parameter-dependent
%   polynomial whose monomials are the trace of the respective matrix
%   monomial.
%
%   trace(poly,label) returns the polynomial, whose monomials are the
%   computed traces, with the informed label.

%   Author: Cristiano M. Agulhari
%   2016, Feb, 4

if (nargin > 2)
    error('Input error. Type ''help trace'' for more details');
    return
end

poly = varargin{1};
if (nargin == 2)
    label = varargin{2};
else
    label = [];
end
    

if(isa(poly,'rolmipvar'))
    if (isempty(label))
        label = strcat('Tr_',poly.label);
    end
    if (length(poly.data) == 1) %Not a polynomial
        tr = rolmipvar(trace(poly.data(1).value),label,0,0);
    else
        for cont=1:length(poly.data)
            M{cont} = [];
            for contsimplex=1:length(poly.vertices)
                M{cont}{contsimplex} = poly.data(cont).exponent{contsimplex};
            end
            M{cont}{length(poly.vertices)+1} = trace(poly.data(cont).value);
        end
        
        deg = [];
        for contsimplex=1:length(poly.vertices)
            deg = [deg sum(poly.data(1).exponent{contsimplex})];
        end
        tr = rolmipvar(M,label,poly.vertices,deg);
    end
    
else
    tr = full(sum(diag(poly)));
end


return
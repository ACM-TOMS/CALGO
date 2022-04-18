function [varargout] = affine2poly(varargin)
% Procedure to compose the polynomial structs to use the parser using, as
% input, an affine description of the parameter-dependent variable 
% m(p) = m0 + p1*m1 + ... + pn*mn. The parameters p may be time-varying and
% are not necessarily in the unit simplex.
%
% [poly] = affine2poly(M,label,bounds) returns the 
% polynomial structure for a variable set M, with a given label and given 
% the bounds for the values of the parameters and of their derivatives, if 
% the parameters are time-varying. The variables M1, M2, ... , MN can be 
% concatenated in M (M = [M1 M2 ... MN]) or disposed in cells 
% (M{1} = M1, ... , M{N} = MN). Considering that parameter p1 is bounded by
% a <= p1 <= b and p2 is bounded by c <= p2 <= d, then
% bounds = [a b; c d].
%
% [poly,M] = affine2poly(rows,cols,label,bounds)
% internally defines a set of rows X cols symmetric matrix variables M with 
% a given label and given bounds related to the parameters
%
% [poly,M] = affine2poly(rows,cols,label,parametr,bounds)
% defines a variable: symmetric if parametr = 'symmetric'; full if parametr
% = 'full'; symmetric toeplitz if parametr = 'toeplitz'; symmetric hankel 
% if parametr = 'hankel'; and skew-symmetric if parametr = 'skew'.

if (ischar(varargin{2}))
    M = varargin{1};
    label = varargin{2};
    bounds = varargin{3};
elseif (ischar(varargin{3}))
    M = [];
    bounds = [];
    rows = varargin{1};
    cols = varargin{2};
    label = varargin{3};
    if ((nargin > 4) && (ischar(varargin{4})))
        parametr = varargin{4};
        bounds = varargin{5};
    else
        parametr = 'symmetric';
        bounds = varargin{4};
    end
else
    error('Syntax error. To see the proper input syntax, type ''help affine2poly''');
    return;
end

numparam = size(bounds,1);
variable = false;
if isempty(M)
    %Define the variable, based on the number of parameters
    variable = true;
    for cont=1:numparam+1
        M{cont} = sdpvar(rows,cols,parametr);
    end
end

if (numparam == 0)
    poly = poly_struct(M,label,0,0);
    varargout{1} = poly;
    if (variable)
        varargout{2} = M;
    end
    return
else
    vertices = 2*ones(1,numparam);
end

if (~iscell(M))
    %Transform to cell array
    [order,tam] = size(M);
    cols = tam/(numparam+1);
    if (round(cols) ~= cols) 
        error('There is an incorrect number of bounds, or the matrix dimensions are incorrect');
        return
    end
    for cont=1:numparam+1
        Maux{cont} = M(:,(cont-1)*cols+1:cont*cols);
    end
    clear M;
    M = Maux;
    clear Maux;
else
    if (length(M) ~= numparam+1)
        error('There is an incorrect number of bounds');
        return
    end
end


% Defining the monomials of the multi-simplex
for cont=1:numparam
    exponents{cont} = [1 0; 0 1];
end

numcoefs = 2^numparam;
contexp = ones(1,numparam+1);
T = [];
for cont=1:numcoefs
    Taux = M{1};
    for contsimplex=1:numparam
        T{cont}{contsimplex} = exponents{contsimplex}(contexp(contsimplex),:);
        Taux = Taux + ((contexp(contsimplex)-1)*bounds(contsimplex,2) + (2-contexp(contsimplex))*bounds(contsimplex,1))*M{contsimplex+1};
    end
    T{cont}{numparam+1} = Taux;
    
    contexp(1) = contexp(1) + 1;
    aux = 1;
    while ((aux <= numparam) && (contexp(aux) > size(exponents{aux},1)))
        contexp(aux) = 1;
        contexp(aux+1) = contexp(aux+1) + 1;
        aux = aux + 1;
    end
end

poly = poly_struct(T,label,vertices,ones(1,numparam));
if (sum(sum(isnan(double(M{1})))) > 0)
    for cont=1:length(poly.opcode)
        ind = find(poly.opcode{cont} == '#');
        poly.opcode{cont}(ind+1) = 'F'; % code for affine variables
    end
else
    for cont=1:length(poly.opcode)
        ind = find(poly.opcode{cont} == '#');
        poly.opcode{cont}(ind+1) = 'Z'; % code for affine matrices
    end
end
poly.nummonorig = numparam+1;
rolmip('setvar',poly);

varargout{1} = poly;
if (variable)
    vargout{2} = T;
end
return

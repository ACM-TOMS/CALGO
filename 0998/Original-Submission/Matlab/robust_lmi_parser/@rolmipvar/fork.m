function [varargout] = fork(varargin)
% Procedure to generate a polynomial variable with the same coefficients
% as the input polynomial, but on a different (multi)simplex domain.
%
% [polyout] = fork(polyin,newlabel) changes the entire (multi)simplex domain. For 
% example, if the input polynomial is defined over the simplexes [1 2 3],
% then the output polynomial will be defined over the simplexes [4 5 6].
% The label of the output polynomial is given by the string newlabel.
%
% [polyout] = fork(polyin,newlabel,targetin) receives a vector target, informing which
% simplexes are to be transformed. For example, if the input polynomial is
% defined over the simplexes [1 2 3], then the output polynomial using
% targetin = [2] will be defined over the simplexes [1 4 3].
%
% [polyout] = fork(polyin,newlabel,targetin,targetout) transforms the simplexes
% informed by targetin to the simplexes given by targetout. For example,
% if the input polynomial is defined over the simplexes [1 2 3], then the
% output polynomial using targetin = [1 3] and targetout = [3 5] will be
% defined over the simplexes [2 3 5].

if (nargin > 4)
    error('Input error. Type ''help fork'' for more details');
    return
end

polyin = varargin{1};
newlabel = varargin{2};
numsimplexes = length(polyin.vertices);
if (nargin >= 3)
    targetin = varargin{3};
    if (nargin == 4)
        targetout = varargin{4};
        if (length(targetin) ~= length(targetout))
            error('The inputs targetin and targetout must have the same number of elements');
            return
        end
    else
        targetout = numsimplexes+1:numsimplexes+length(targetin);
    end
else
    targetin = 1:numsimplexes;
    targetout = numsimplexes+1:numsimplexes+length(targetin);
end

for cont = 1:length(targetout)
    %Verify the number of vertices
    if ((targetout(cont) <= length(polyin.vertices)) && (polyin.vertices(targetin(cont)) ~= polyin.vertices(targetout(cont))) && (polyin.vertices(targetout(cont)) > 0))
        error('The target simplex must have the same number of vertices of the original simplex');
        return
    end
end

for cont = 1:length(polyin.data)
    M{cont} = [];
    for contsimplex = 1:length(polyin.vertices)
        ind = find(targetin == contsimplex);
        if (isempty(ind))
            M{cont}{contsimplex} = polyin.data(cont).exponent{contsimplex};
        else
            M{cont}{targetout(ind)} = polyin.data(cont).exponent{targetin(ind)};
            if (isempty(M{cont}{contsimplex}))
                M{cont}{contsimplex} = 0;
            end
        end
    end
    M{cont}{length(M{cont})+1} = polyin.data(cont).value;
end

for contsimplex = 1:length(polyin.vertices)
    ind = find(targetin == contsimplex);
    if (isempty(ind))
        vertices(contsimplex) = polyin.vertices(contsimplex);
        degree(contsimplex) = sum(polyin.data(1).exponent{contsimplex});
    else
        vertices(targetout(ind)) = polyin.vertices(targetin(ind));
        degree(targetout(ind)) = sum(polyin.data(1).exponent{targetin(ind)});
    end
end

if (isfield(polyin,'nummonorig'))
    nummonorig = polyin.nummonorig;
else
    nummonorig = -1;
end


%Regularize variable M
for cont1 = 1:length(M)
    for cont2 = 1:length(M{cont1})
        if (isempty(M{cont1}{cont2}))
            M{cont1}{cont2} = 0;
        end
    end
end

poly = rolmipvar(M,newlabel,vertices,degree);

if (nummonorig ~= -1)
    poly.nummonorig = nummonorig;
end
varargout{1} = poly;
return
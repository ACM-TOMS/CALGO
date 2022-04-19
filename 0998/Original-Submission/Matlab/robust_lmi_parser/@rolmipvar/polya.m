function polyout = polya(varargin)
% Implementation of the Polya's Theorem, used to reduce the
% conservativeness for the verification of polynomial positivity. 
%
% [polyout] = polya(polyin,d) applies the Polya relaxation of degree d over
% the input polynomial polyin. The input d is a vector, where the i-th
% position corresponds to the degree of the relaxation to be applied on the
% i-th simplex.
%
% [polyout] = polya(polyin) applies the Polya relaxation of degree 1 over
% all the simplexes of polyin.

if ((nargin > 2) || (nargin < 1))
    error('Input error. Type ''help polya'' for more details');
    return
else
    polyin = varargin{1};
    numsimplexes = length(polyin.data(1).exponent);
    if (nargin == 2)
        d = varargin{2};
    else
        d = ones(1,numsimplexes);
    end
    
    if (length(d) > numsimplexes)
        error('The size of d cannot be higher than the number of simplexes of polyin');
    end
    
    if (polyin.vertices == 0)
        warning('The input is not a polynomial');
    end
end


auxvet = [1 zeros(1,numsimplexes-1)];
polyout = polyin;
for cont = 1:length(d)
    %Construct unitary polynomial
    clear M;
    vertone = [1 zeros(1,polyin.vertices(cont)-1)];
    for contmon = 1:polyin.vertices(cont)
        for contsimp=1:numsimplexes
            M{contmon}{contsimp} = [0];
        end
        M{contmon}{cont} = vertone;
        M{contmon}{numsimplexes+1} = 1;
        vertone = circshift(vertone,[0 1]);
    end

    aux = rolmipvar(M,'1',auxvet*polyin.vertices(cont),auxvet);
    for contpolya = 1:d(cont)
        polyout = polyout*aux;
    end
    auxvet = circshift(auxvet,[0 1]);
end


return
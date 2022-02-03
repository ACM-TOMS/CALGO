function newPoly = powPoly(poly, pow)
%%POWPOLY   Compute power of polynomial.
%   Usage: 
%      newPoly = powPoly(poly, pow);

% newPoly = poly;
% logpow = floor(log2(pow));
% rempow = pow - 2^logpow;
% for ii = 1:logpow
%     newPoly = multiplyPoly(newPoly, newPoly);
% end
% for ii = 1:rempow
%     newPoly = multiplyPoly(newPoly, poly);
% end
    if ~checkPoly(poly)
        error('Input format is not correct');
    end
    if ceil(pow) ~= floor(pow) || pow < 0
        error('pow must be a non-negative integer.');
    end

    nn = size(poly.supports, 2);
    newPoly.supports = sparse(1, nn);
    newPoly.coef = 1;
    if pow == 0
        return;
    else
        bin = dec2bin(pow);
        if strcmp(bin(end), '1')
            newPoly = poly;
        end
        tmpPoly = poly;
        for ii = (length(bin)-1):-1:1
            tmpPoly = multiplyPoly(tmpPoly, tmpPoly, 0);
            if strcmp(bin(ii), '1')
                newPoly = multiplyPoly(newPoly, tmpPoly, 0);
            end
        end
    end
end


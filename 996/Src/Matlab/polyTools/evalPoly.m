function val = evalPoly(poly, x)
%%EVALPOLY   Evaluate the value of a polynomial at x
%   Usage:
%      val = evalPoly(poly, x);
%   Input:
%      poly: polynomial in sparsePOP format
%      x   : a vector

if ~checkPoly(poly)
    error('Input format is not correct');
end

constIdx = ~any(poly.supports,2);
approxmono = double(~constIdx);
for kk = 1:length(x)
    approxmono = approxmono .* x(kk).^poly.supports(:,kk);
end
val = poly.coef'*(approxmono+constIdx);

end%function
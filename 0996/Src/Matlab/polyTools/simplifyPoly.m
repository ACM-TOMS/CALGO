function simpPoly = simplifyPoly(poly, I01, delSW, sorted)
%%SIMPLIFYPOLY   Simplify polynomial
%
%   simpPoly = simplifyPoly(poly);
%   simpPoly = simplifyPoly(poly,I01);
%   simpPoly = simplifyPoly(poly,I01,delSW);
%
%   Input:
%      objPoly:  polynomial in sparsePOP format.
%      I01    :  Indices of 0-1 variables.
%      delSW  :  if true, this function deletes unused variables.

if nargin > 4; error('Number of input argument must be at most 3.'); end;
    
if ~(isempty(poly) || checkPoly(poly))
    error('Input format is not correct');
end
if isempty(poly); simpPoly = poly; return; end;

if nargin <= 3; sorted = true; end;
if (nargin <= 2 || isempty(delSW)); delSW = false; end;
if (nargin <= 1 || isempty(I01)); I01 = false(1, size(poly.supports, 2)); end;
if length(I01) ~= size(poly.supports, 2); error('dimensions must agree.'); end;

tmpSupports = poly.supports;

tmpSupports = sparse(tmpSupports);
tmpSupports(:, I01) = spones(tmpSupports(:, I01));

if sorted
    [simpPoly.supports, ~, ic] = grevlexUnique(tmpSupports);
else
    [simpPoly.supports, ~, ic] = fastUnique(tmpSupports);
end
simpPoly.coef = accumarray(ic,poly.coef); %full(sparse(ic,1,poly.coef));

remainIdx = simpPoly.coef ~= 0;
simpPoly.coef = simpPoly.coef(remainIdx);
simpPoly.supports = simpPoly.supports(remainIdx, :);

if delSW
    simpPoly.supports = simpPoly.supports(:, any(simpPoly.supports, 1));
end

end
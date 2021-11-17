function objPoly = quad2Poly(A, bb, c)
%%QUAD2POLY  Convert quadratic form into simplified sparsePOP format
%
%  f(x) = x'*A*x  +  bb'*x  +  c
%
%  Usage:
%    objPoly = quad2Poly(A);
%    objPoly = quad2Poly(A, bb);
%    objPoly = quad2Poly(A, [], c);
%    objPoly = quad2Poly(A, bb, c);


if nargin < 2; bb = []; end;
if nargin < 3; c = []; end;

dimVar = max([size(A), length(bb(:))]);
if (~isempty(A) && (dimVar ~= size(A, 1) || dimVar ~= size(A, 2))) ...
        || (~isempty(bb) && length(bb(:)) ~= dimVar) ...
        || (~isempty(c) && ~isscalar(c))
    error('Input is not correct');
end

if isempty(A)
    supportsA = [];
    coefA = [];
else
    U = triu(A, 1)' + tril(A);
    numTermA = nnz(U);

    [row, col, val] = find(U);

    suppIdx = repmat(1:numTermA, 2, 1);
    supportsA = sparse(suppIdx, [row'; col'], 1, numTermA, dimVar);
    coefA = val;
end

if isempty(bb)
    supportsbb = [];
    coefbb = [];
else
    supportsbb = speye(dimVar);
    coefbb = bb(:);
end

if isempty(c)
    supportsc = [];
    coefc = [];
else
    supportsc = zeros(1, dimVar);
    coefc = c;
end

objPoly.supports = [supportsc; supportsbb; supportsA];
objPoly.coef = [coefc; coefbb; coefA];

end


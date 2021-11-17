function [objPoly,ineqPolySys,lbd,ubd] = nondquar(nDim); 

% nondquar(nDim);
% 
% The function in nondquar.mod 
% (http://www.sor.princeton.edu/%7Ervdb/ampl/nlmodels/cute/nondquar.mod)
%
% <Input> 
% nDim: This is the dimension of the function
% Remark: nDim needs to be a multiple of 4. 
%
% <Output>
% objPoly,ineqPolySys,lbd,ubd

if mod(nDim,4) ~= 0 
	error('## nDim needs to be a multiple of 4 ##');
end

%Objective Function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

poly = cell(nDim,1);
poly{1}.typeCone = 1;
poly{1}.sizeCone = 1;
poly{1}.dimVar   = nDim;
poly{1}.degree   = 1;
poly{1}.noTerms  = 2;
poly{1}.supports = sparse(poly{1}.noTerms,poly{1}.dimVar);
poly{1}.supports(1,1)   = 1;
poly{1}.supports(2,2)   = 1;
poly{1}.coef = [1;-1];

for j=1:nDim-2
  poly{j+1}.typeCone = 1;
  poly{j+1}.sizeCone = 1;
  poly{j+1}.dimVar   = nDim;
  poly{j+1}.degree   = 2;
  poly{j+1}.noTerms  = 6;
  poly{j+1}.supports = sparse(poly{j+1}.noTerms,poly{j+1}.dimVar);
  poly{j+1}.supports(1,j)   = 2;
  poly{j+1}.supports(2,j+1) = 2;
  poly{j+1}.supports(3,nDim)= 2;
  poly{j+1}.supports(4,j)   = 1;
  poly{j+1}.supports(4,j+1) = 1;
  poly{j+1}.supports(5,nDim)= 1;
  poly{j+1}.supports(5,j)   = 1;
  poly{j+1}.supports(6,j+1) = 1;
  poly{j+1}.supports(6,nDim)= 1;
  poly{j+1}.coef = [1;1;1;2;2;2];
end
poly{nDim}.typeCone = 1;
poly{nDim}.sizeCone = 1;
poly{nDim}.dimVar   = nDim;
poly{nDim}.degree   = 1;
poly{nDim}.noTerms  = 2;
poly{nDim}.supports = sparse(poly{nDim}.noTerms,poly{nDim}.dimVar);
poly{nDim}.supports(1,nDim-1)   = 1;
poly{nDim}.supports(2,nDim)   = 1;
poly{nDim}.coef = [1;1];

objPoly.typeCone = 1;
objPoly.sizeCone = 1;
objPoly.dimVar   = nDim;
objPoly.degree   = 4;
objPoly.noTerms  = 0;
objPoly.supports = [];
objPoly.coef     = [];
for i=1:nDim
    temp = multiplyPolynomials(poly{i}, poly{i});
    clear poly{i}
    objPoly.noTerms = objPoly.noTerms + temp.noTerms;
    objPoly.supports= [objPoly.supports;temp.supports];
    objPoly.coef    = [objPoly.coef;temp.coef];
end

objPoly = simplifyPolynomial(objPoly);
ineqPolySys = [];
lbd = -1.0e+10*ones(1,nDim);
ubd  = 1.0e+10*ones(1,nDim);

return;

function polyOut = multiplyPolynomials(poly1,poly2, order)
%
% Multiplication of two real valued polynomials 
%       polynomial.nDim --- a positive integer 
%       polynomial.degree --- a positive integer 
%       polynomial.noTerms --- a positive integer 
%       polynomial.supports --- a noTerms \times nDim matrix
%       polynomial.coef --- a noTerms \times 1 vector

%Now, this function can deal with only typeCone = 1 case.
%2005 Jan 9

if nargin < 3
  order = 'grevlex';
end

if poly1.dimVar ~= poly2.dimVar
  error('!!! poly1.nDim ~= poly2.nDim !!!');
end

polyIn.typeCone = poly1.typeCone;
polyIn.sizeCone = poly1.sizeCone;
polyIn.dimVar = poly1.dimVar;
polyIn.noTerms = poly1.noTerms * poly2.noTerms; 
polyIn.supports = sparse(polyIn.noTerms,polyIn.dimVar);
UsedVariable = find(any([poly1.supports;poly2.supports],1));

idx1 = repmat((1:poly1.noTerms),poly2.noTerms,1);
idx1 = idx1(:);
FirSup = poly1.supports(idx1,UsedVariable);
idx2 = (1:poly2.noTerms)';
idx2 = repmat(idx2,poly1.noTerms,1);
SecSup = poly2.supports(idx2,UsedVariable);
polyIn.supports(:,UsedVariable) = FirSup + SecSup;

polyIn.coef = poly2.coef * poly1.coef';
polyIn.coef = polyIn.coef(:);
    
degree = sum(polyIn.supports(:,UsedVariable),2);
polyIn.degree = max(degree);
polyOut = simplifyPolynomial(polyIn,'grevlex');
clear degree polyIn UsedVariable idx1 idx2 FirSup SecSup

return;


% $Header: /home/waki9/CVS_DB/SparsePOPdev/example/POPformat/nondquar.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $

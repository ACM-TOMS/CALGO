function [objPoly,ineqPolySys,lbd,ubd] = BroydenTri(nDim);

% [] = solveBroydenTri(nDim);
% 
% The Broyden Tridiagonal function, which is
% described in "Testing Unconstrained Optimization Software",
% J.J.More et.al, ACM Trans. Math. Soft., 7, p.17-41
%
% <Input> 
% nDim: The dimension of the function
%
% <Output>
% objPoly,ineqPolySys,lbd,ubd

%Objective Function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
poly{1}.typeCone = 1;
poly{1}.sizeCone = 1;
poly{1}.dimVar   = nDim;
poly{1}.degree   = 2;
poly{1}.noTerms  = 4;
poly{1}.supports = sparse(poly{1}.noTerms,poly{1}.dimVar);
poly{1}.supports(1,1) = 1;
poly{1}.supports(2,1) = 2;
poly{1}.supports(3,2) = 1;
poly{1}.coef     = [3;-2;-2;1];

for i=2:nDim-1
  poly{i}.typeCone = 1;
  poly{i}.sizeCone = 1;
  poly{i}.dimVar   = nDim;
  poly{i}.degree   = 2;
  poly{i}.noTerms  = 5;
  poly{i}.supports = sparse(poly{i}.noTerms,poly{i}.dimVar);
  poly{i}.supports(1,i) = 1;
  poly{i}.supports(2,i) = 2;
  poly{i}.supports(3,i-1) = 1;
  poly{i}.supports(4,i+1) = 1;
  poly{i}.coef     = [3;-2;-1;-2;1];
end
poly{nDim}.typeCone = 1;
poly{nDim}.sizeCone = 1;
poly{nDim}.dimVar   = nDim;
poly{nDim}.degree   = 2;
poly{nDim}.noTerms  = 4;
poly{nDim}.supports = sparse(poly{nDim}.noTerms,poly{nDim}.dimVar);
poly{nDim}.supports(1,nDim) = 1;
poly{nDim}.supports(2,nDim) = 2;
poly{nDim}.supports(3,nDim-1) = 1;
poly{nDim}.coef     = [3;-2;-1;1];

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
lbd(1,1) = 0;
ubd  =  1.0e+10*ones(1,nDim);
return;

% $Header: /home/waki9/CVS_DB/SparsePOPdev/example/POPformat/BroydenTri.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $

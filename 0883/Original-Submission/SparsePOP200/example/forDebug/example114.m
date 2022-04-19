function [objPoly,ineqPolySys,lbd,ubd] = example114(nDim);

%
% Linearly dependent coefficient matrix, infeasible
%
% <Input> 
% nDim: The dimension of the function
%
% <Output>
% objPoly,ineqPolySys,lbd,ubd

%Objective Function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rand('state',3201);
objPoly.typeCone = 1;
objPoly.sizeCone = 1;
objPoly.dimVar   = nDim;
objPoly.degree   = 1;
objPoly.noTerms  = nDim;
objPoly.supports = speye(nDim,nDim);
objPoly.coef     = ones(nDim,1);

for i=1:nDim
    ineqPolySys{i}.typeCone = -1;
    ineqPolySys{i}.sizeCone = 1;
    ineqPolySys{i}.dimVar   = nDim;
    ineqPolySys{i}.degree   = 1;
    ineqPolySys{i}.noTerms  = 2;
    ineqPolySys{i}.supports = [sparse(2,nDim)];
    ineqPolySys{i}.supports(2,i) = 1;
    ineqPolySys{i}.coef     = ones(2,1);
end

lbd = zeros(1,nDim);
ubd  =  1.0e+10*ones(1,nDim);
return;

% $Header: /home/waki9/CVS_DB/SparsePOPdev/example/POPformat/BroydenTri.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $

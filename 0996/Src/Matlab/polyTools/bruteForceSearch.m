function [ val,sol ] = bruteForceSearch( objPoly, binaryMode )
%%BRUTEFORCESEARCH    Minimize a polynomial over a binary constraint.
%    objPoly: Polynomial (in sparsePOP format) to minimize.
%    binaryMode: 
%        0 -- minimize objPoly on  {0,1}^n
%       -1 -- minimize objPoly on {-1,1}^n

if ~checkPoly(objPoly)
    error('bruteForceSearch: Input format is not correct');
end

numVar = size(objPoly.supports,2);
allSol = full(enum01MonosUptoDeg(numVar,numVar)');
if binaryMode == -1
    allSol(allSol==0) = -1; 
elseif binaryMode ~= 0
    error('bruteForceSearch: binaryMode must be 0 or -1.');
end
val = Inf; sol = zeros(numVar,1);
for ii = 1:size(allSol,2)
    tmpsol = allSol(:,ii);
    tmpval = evalPoly(objPoly,tmpsol);
    if tmpval < val
        sol = tmpsol;
        val = tmpval;
    end
end

end

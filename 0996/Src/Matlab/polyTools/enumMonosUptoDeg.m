function supports = enumMonosUptoDeg( numVar, maxDegree )
%ENUMMONOSUPTODEG   Enumerate monomials upto maxDegree.
%   Input:
%      numVar:  number of variables
%      maxDegree:  degree of monomials
%   Output:
%      supports: 2 dimensional array whose rows represent 
%                the supports of monomials.

if maxDegree == 0
    supports = sparse(1,numVar);
    return;
end

Col = enumCombWithRep(numVar+1,maxDegree) + 1;
Row = repmat(1:size(Col,2),maxDegree,1);

% The below is faster than sparse(Row(:),Col(:),1);
% because MATLAB uses the column major order
% and thus the second argument should be sorted when constructing a matrix.
supports = sparse(Col(:),Row(:),1)';
supports(:,1) = [];

% basisSupports = grevlexSort(basisSupports); % already sorted

end


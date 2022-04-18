function supports = enum01MonosSameDeg( numVar, degree )
%%ENUM01MONOSSAMEDEG   Enumerate monomials of binary variable with same degree.
%   Input:
%      numVar:  number of variables
%      degree:  degree of monomials
%   Output:
%      supports: 2 dimensional array whose rows represent 
%                the supports of monomials.

if degree == 0
    supports = sparse(1,numVar);
    return;
end

Col = enumCombWithNoRep(numVar,degree) + 1;
Row = repmat(1:size(Col,2),degree,1);

% The below is faster than sparse(Row(:),Col(:),1);
% because MATLAB uses the column major order
% and thus the second argument should be sorted when constructing a matrix.
supports = sparse(Col(:),Row(:),1)';

%supports = grevlexSort(supports); % already sorted

end


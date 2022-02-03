function supports = enum01MonosUptoDeg( numVar, maxDegree )
%%enum01MonosUptoDeg Enumerate monomials of binary variable up to maxDegree.
%   Input:
%      numVar:     number of variables
%      maxDegree:  maximum degree of monomials
%   Output:
%      supports: 2 dimensional array whose rows represent 
%                the supports of monomials.

Col = [];
Row = [];
len = 0;
for ii = 1:maxDegree
    iiCol = enumCombWithNoRep(numVar,ii) + 1;
    iiRow = ones(ii,1) * ( (1:size(iiCol,2)) + len ); % faster than repmat()
    len = len + size(iiCol,2);
    Col = [Col; iiCol(:)];
    Row = [Row; iiRow(:)];
end

% The below is faster than sparse(Row(:),Col(:),1);
% because MATLAB employs the column major order
% and thus the second argument should be sorted when constructing a matrix.
supports = [zeros(1,numVar); sparse(Col(:),Row(:),1)']; 

%supports = grevlexSort(supports); %already sorted

end


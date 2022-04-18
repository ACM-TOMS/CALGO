function X=ctranspose(atual)
%TRANSPOSE (overloaded)
%
% Author: Alexandre Felipe
% 2014, Dec, 8
X = atual;
if ((~isempty(find(X.label == '+',1))) || (~isempty(find(X.label == '-',1))) || (~isempty(find(X.label == '*',1))))
    %The label corresponds not only to a variable, but to an expression
    pre = '(';
    aft = ')';
else
    pre = '';
    aft = '';
end

for i = 1:length(atual.data)
    X.data(i).value = (atual.data(i).value)';
    X.opcode{i} = strcat(pre,atual.opcode{i}, aft, char(39));
end
X.label = [pre,atual.label,aft, char(39)];


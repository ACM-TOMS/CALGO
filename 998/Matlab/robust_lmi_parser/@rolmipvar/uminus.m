function Z = uminus(X, Y)
%MINUS (overloaded)
%
% Author: Alexandre Felipe
% 2014, Dec, 8
Z = X;
for i = 1:length(X.data)
  Z.data(i).value = -(X.data(i).value);
  Z.opcode{i} = strcat('-', X.opcode{i});
end
Z.label = strcat('-', X.label);

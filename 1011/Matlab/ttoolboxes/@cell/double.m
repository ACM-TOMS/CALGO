function C = double(C)
% C = double(C)
% Converts entries of a cell array to double.
%
% E.g.: double({sym(1) {sym(2)}})
%
% Written by tommsch, 2019

for i=1:numel(C);
    C{i}=double(C{i});
end

end
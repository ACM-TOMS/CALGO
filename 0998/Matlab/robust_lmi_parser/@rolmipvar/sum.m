function Z = sum(X)
%SUM (overloaded)
%
% Author: Cristiano M. Agulhari
% 2016, Nov, 18

Z = 0;
for cont=1:length(X.data)
    Z = Z + X.data(cont).value;
end

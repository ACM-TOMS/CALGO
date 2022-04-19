function Y = value(X)
%DOUBLE (overloaded)
%
% Author: Alexandre Felipe
% 2014, Dec, 17
%
% Retrieve the value of the a variable, if the variable has degree 0
% returns its coefficient value, otherwise return a rolmipvar
% with each coefficient evaluated.
%
% This function is intented to convert rolmipvar whose coefficients
% are variables (e.g yalmip sdpvar) to a rolmipvar with constant
% coefficient. This must ease multistage methods.
if(length(X.data) == 1)
    Y = value(X.data(1).value);
else
    Y = X;
    for i = 1:length(Y.data)
        Y.data(i).value = double(X.data(i).value);
    end
end
    
        
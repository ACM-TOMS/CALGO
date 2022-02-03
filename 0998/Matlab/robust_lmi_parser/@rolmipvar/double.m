function S = double(X)
%DOUBLE (overloaded)
%
% Author: Alexandre Felipe
% 2014, Dec, 9
%
% Update: Cristiano M. Agulhari
% 2016, Feb, 16
%
% Retrieve the value of the a variable if it is a polynomial of degree 0,
% or the polynomial of doubles if they were SDPVAR variables.
 
  if length(X.data) == 1
      S = double(X.data(1).value);
  else
      S = X;
      for cont=1:length(S.data)
          S.data(cont).value = double(S.data(cont).value);
      end
  end
return 
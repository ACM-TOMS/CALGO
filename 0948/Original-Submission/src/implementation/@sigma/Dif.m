function y = Dif(y, d)

if nargin == 1, d = 1; end


y.vector = y.vector + d;

end
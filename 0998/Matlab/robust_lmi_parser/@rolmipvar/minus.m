function Z = minus(X, Y)
%MINUS (overloaded)
%
% Author: Alexandre Felipe
% 2014, Dec, 8
h = homogenize(X, Y);
Zf = operation_poly(h{1}, h{2}, '-');
Z = rolmipvar(Zf);


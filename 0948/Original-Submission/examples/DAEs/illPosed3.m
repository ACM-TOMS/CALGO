function f = illPosed3(t, z)
% f = illPosed3(t, z) evaluates the ill-posed problem, which has an over-
% determined part, a well-determined part, and an under-determined part.
% This problem is used in Section 5.2 of R. McKenzie, J. Pryce, 
% N. Nedialkov, G. Tan "DAESA User Guide".
%
% Copyright 2014 N. Nedialkov, J. Pryce, G. Tan

x1 = z(1); x2 = z(2); x3 = z(3);
x4 = z(4); x5 = z(5); x6 = z(6);


f(1) = x1   + x2 +           x5   + x6;
f(2) = x1^3 + x2 +           x5   + x6  ;
f(3) =                       x5   * x6  ;
f(4) =                     - x5^3 + x6^4;
f(5) = x1   + x2 + x3 + x4 + x5   + x6^3;
f(6) =                       x5   + x6  ;

end

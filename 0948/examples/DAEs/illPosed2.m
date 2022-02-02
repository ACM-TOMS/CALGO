function f = illPosed2(t,z,G,L,c)
% f = illPosed2(t,z,G,L,c) evaluates the modified two-pendula problem,
% where equation 3 is missing, and variable 6, mu, is missing.  This
% problem is used in Section 5.2 of R. McKenzie, J. Pryce, N. Nedialkov, 
% G. Tan "DAESA User Guide".
%
% Copyright 2014 N. Nedialkov, J. Pryce, G. Tan

% rename for better readability
x = z(1); y = z(2); la = z(3);
u = z(4); v = z(5); mu = z(6);
% first pendulum
f(1) = Dif(x,2)+x*la;        % $x'' + x\lambda =0$
f(2) = Dif(y,2)+y*la-G;      % $y'' + y\lambda -G=0$
%f(3) = x^2+y^2-L^2;          % $x^2+ y^2  - L^2 = 0$
% modified second pendulum
f(4) = Dif(u,2)  +u;      % $u'' + u\mu =0$
f(5) = Dif(v,3)^2+v-G;    % $(v''')^2 + v\mu - G = 0$

f(6) = u^2+v^2-(L+c*la)^2+Dif(la,2);
% $u^2 + v^2  - (L+c\lambda)^2 + \lambda'' = 0$
end


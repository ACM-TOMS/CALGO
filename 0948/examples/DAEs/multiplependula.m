function f = multiplependula (t,x,p)
% f = multiplependula (t,x,p) evaluates the multiple chained pendula 
% problem, used in Section 5.4 of R. McKenzie, J. Pryce, N. Nedialkov, 
% G. Tan "DAESA User Guide".
%
% Copyright 2014 N. Nedialkov, J. Pryce, G. Tan

G = 9.8; L = 1; c = 0.1;
% first pendulum
f(1) = Dif(x(1),2) + x(1)*x(3);
f(2) = Dif(x(2),2) + x(2)*x(3) -G;
f(3) = x(1)^2 + x(2)^2 - L^2;

% pendulum > 1
for i = 2:p
    xi = 3*i-2; yi = 3*i-1;     % x,y indices
    li = 3*i; li1 = 3*i-3;      % lam_i, lam_{i-1} indices
    f(xi) = Dif(x(xi),2) + x(xi)*x(li); 
    f(yi) = Dif(x(yi),2) + x(yi)*x(li) - G;
    f(li) = x(xi)^2 + x(yi)^2 - (L+c*x(li1))^2;
end
end
function[a,b,c,d,e]=coefficient_periodic_second(p)
% M. Mehra & K. S. Patel
% Inputs
% p=order of accuracy (p=4 or 6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% It computes the coefficents for second derivative approximations for
% interior points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms a b c d e
eqn1=a+b+c-2*d-2*e-1;
eqn2=a+(2^2)*b+(3^2)*c-(factorial(4)/factorial(2))*d-(factorial(4)/factorial(2))*(2^2)*e;
eqn3=a+(2^4)*b+(3^4)*c-(factorial(6)/factorial(4))*d-(factorial(6)/factorial(4))*(2^4)*e;
if p==4
d=1/10;e=0;c=0;
[a,b]=solve(eqn1,eqn2,'a,b');
a=eval(a);b=eval(b);
end
if p==6
d=2/11;e=0;
[a,b,c]=solve(eqn1,eqn2,eqn3,'a,b,c');
 a=eval(a);b=eval(b);c=eval(c);
end
end

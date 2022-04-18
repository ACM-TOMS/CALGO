function[a,b,c,d,e]=coefficient_periodic_second(p)
% M. Mehra & K. S. Patel
% Inputs
% p=order of accuracy (p=4 or 6 or 8 or 10)
% Output
% It computes the coefficents for second derivative approximations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms a b c d e
eqn1=a+b+c-2*d-2*e-1;
eqn2=a+(2^2)*b+(3^2)*c-(factorial(4)/factorial(2))*d-(factorial(4)/factorial(2))*(2^2)*e;
eqn3=a+(2^4)*b+(3^4)*c-(factorial(6)/factorial(4))*d-(factorial(6)/factorial(4))*(2^4)*e;
eqn4=a+(2^6)*b+(3^6)*c-(factorial(8)/factorial(6))*d-(factorial(8)/factorial(6))*(2^6)*e;
eqn5=a+(2^8)*b+(3^8)*c-(factorial(10)/factorial(8))*d-(factorial(10)/factorial(8))*(2^8)*e;
if p==4
d=1/10;e=0;c=0;
[a,b]=solve(eqn1,eqn2,'a,b');a=eval(a);b=eval(b);
end
if p==6
d=2/11;e=0;
[a,b,c]=solve(eqn1,eqn2,eqn3,'a,b,c');a=eval(a);b=eval(b);c=eval(c);
end
if p==8
d=344/1179;
[a,b,c,e]=solve(eqn1,eqn2,eqn3,eqn4,'a,b,c,e');a=eval(a);b=eval(b);c=eval(c);e=eval(e);
end
if p==10
[a,b,c,d,e]=solve(eqn1,eqn2,eqn3,eqn4,eqn5,'a,b,c,d,e');a=eval(a);b=eval(b);c=eval(c);d=eval(d);e=eval(e);
end
end

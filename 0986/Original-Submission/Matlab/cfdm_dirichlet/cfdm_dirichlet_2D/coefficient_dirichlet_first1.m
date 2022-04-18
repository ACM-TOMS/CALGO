function[a1,a2,a3,a4,a5,a6,a7,a8,alp]=coefficient_dirichlet_first1(p)
% M. Mehra & K. S. Patel
% Inputs
% p=order of accuracy (p=4 or 6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% It computes the coefficents for i=1 and i=N (first and last grid point) for first derivative approximations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms a1 a2 a3 a4 a5 a6 a7 a8 alp 
eqn1=a1+a2+a3+a4+a5+a6+a7+a8;
eqn2=a2+(2)*a3+(3)*a4+(4)*a5+(5)*a6+(6)*a7+(7)*a8-alp-1;
eqn3=a2+(2^2)*a3+(3^2)*a4+(4^2)*a5+(5^2)*a6+(6^2)*a7+(7^2)*a8-(factorial(2)/factorial(1))*alp;
eqn4=a2+(2^3)*a3+(3^3)*a4+(4^3)*a5+(5^3)*a6+(6^3)*a7+(7^3)*a8-(factorial(3)/factorial(2))*alp;
eqn5=a2+(2^4)*a3+(3^4)*a4+(4^4)*a5+(5^4)*a6+(6^4)*a7+(7^4)*a8-(factorial(4)/factorial(3))*alp;
eqn6=a2+(2^5)*a3+(3^5)*a4+(4^5)*a5+(5^5)*a6+(6^5)*a7+(7^5)*a8-(factorial(5)/factorial(4))*alp;
eqn7=a2+(2^6)*a3+(3^6)*a4+(4^6)*a5+(5^6)*a6+(6^6)*a7+(7^6)*a8-(factorial(6)/factorial(5))*alp;
if p==4
a5=0;a6=0;a7=0;a8=0;
[a1,a2,a3,a4,alp]=solve(eqn1,eqn2,eqn3,eqn4,eqn5,'a1,a2,a3,a4,alp');
   a1=eval(a1);a2=eval(a2);a3=eval(a3);a4=eval(a4);alp=eval(alp);
end
if p==6;
a7=0;a8=0;
[a1,a2,a3,a4,a5,a6,alp]=solve(eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,'a1,a2,a3,a4,a5,a6,alp');
  a1=eval(a1);a2=eval(a2);a3=eval(a3);a4=eval(a4);a5=eval(a5);a6=eval(a6);alp=eval(alp);
end
end
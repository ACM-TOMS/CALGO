function[b1,b2,b3,b4,b5,b6,b7,b8,b9,bet]=coefficient_dirichlet_second2(p)
% M. Mehra & K. S. Patel
% Inputs
% p=order of accuracy (p=4 or 6)
% Output
% It computes the coefficents for i=2 and i=N-1 (second and second last grid point) for second derivative approximations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms b1 b2 b3 b4 b5 b6 b7 b8 b9 bet
eqn1=b1+b2+b3+b4+b5+b6+b7+b8+b9;
eqn2=b2+(2)*b3+(3)*b4+(4)*b5+(5)*b6+(6)*b7+(7)*b8+(8)*b9;
eqn3=b2+(2^2)*b3+(3^2)*b4+(4^2)*b5+(5^2)*b6+(6^2)*b7+(7^2)*b8+(8^2)*b9-(factorial(2)/factorial(1))*(1+2*bet);
eqn4=b2+(2^3)*b3+(3^3)*b4+(4^3)*b5+(5^3)*b6+(6^3)*b7+(7^3)*b8+(8^3)*b9-(factorial(3)/factorial(1))*(1+(2)*bet);
eqn5=b2+(2^4)*b3+(3^4)*b4+(4^4)*b5+(5^4)*b6+(6^4)*b7+(7^4)*b8+(8^4)*b9-(factorial(4)/factorial(2))*(1+(2^2)*bet);
eqn6=b2+(2^5)*b3+(3^5)*b4+(4^5)*b5+(5^5)*b6+(6^5)*b7+(7^5)*b8+(8^5)*b9-(factorial(5)/factorial(3))*(1+(2^3)*bet);
eqn7=b2+(2^6)*b3+(3^6)*b4+(4^6)*b5+(5^6)*b6+(6^6)*b7+(7^6)*b8+(8^6)*b9-(factorial(6)/factorial(4))*(1+(2^4)*bet);
eqn8=b2+(2^7)*b3+(3^7)*b4+(4^7)*b5+(5^7)*b6+(6^7)*b7+(7^7)*b8+(8^7)*b9-(factorial(7)/factorial(5))*(1+(2^5)*bet);
if p==6;    
b8=0;b9=0;
[b1,b2,b3,b4,b5,b6,b7,bet]=solve(eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,'b1,b2,b3,b4,b5,b6,b7,bet');
  b1=eval(b1);b2=eval(b2);b3=eval(b3);b4=eval(b4);b5=eval(b5);b6=eval(b6);b7=eval(b7);bet=eval(bet);
end
end

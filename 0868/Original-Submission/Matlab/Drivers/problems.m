function y=problems(x,i)
%function y=problems(x,i)
%This code computes the function values
%to a vector input X to problem number I
%where I is an integer between 1 and 23.

%Note: the function in problem 12 may be evaluated at x=0
%if the rules are closed and thus cause a NaN result. 
%This is handled by e-g. da2glob re-defining y: y(isnan(y))=1;
%A user is normally better off to give the correct value in the function
%definition....it is pure luck that da2glob's  choice follows l'Hospitals
%rule in problem 12 (since the value 1 is used always in NaN situations).
%Similarly is +/-Inf replaced by 0 by these codes.
switch i
  case 1,
    y=exp(x);
  case 2,
    y=(x>0.3);
  case 3,
    y=sqrt(x);
  case 4,
    y=(23/25)*cosh(x)-cos(x);
  case 5
    y=1./(x.^4+x.^2+0.9);
  case 6
    y=sqrt(x.^3);
  case 7
    y=1 ./sqrt(x);   
  case 8
    y=1 ./(1+x.^4);
  case 9
    y=2./(2+sin(10*pi.*x));
  case 10
    y=1 ./(1+x);
  case 11
    y=1 ./(1+exp(x));
  case 12
    y=x./(exp(x)-1);
  case 13
    y=sin(100*pi*x)./(pi*x);
  case 14
    y=sqrt(50).*exp(-50.*pi.*x.^2);
  case 15
    y=25*exp(-25*x);
  case 16
    y=50 ./(pi.*(2500.*x.^2+1));
  case 17
    y=50 .*(sin(50 .*pi.*x)./(50 .*pi.*x)).^2;
  case 18
    y=cos(cos(x)+3*sin(x)+2*cos(2*x)+3*sin(2*x)+3*cos(3*x));
  case 19
    y=zeros(1,length(x));index=(x>10^-15);y(index)=log(x(index));
  case 20
    y=1 ./(x.^2+1.005);
  case 21
    y=1 ./cosh(20 .*(x-0.2)) + 1 ./cosh(400 .* (x-0.4)) + 1 ./cosh(8000 .*(x-.6));
  case 22
    y=4*pi^2 .*x.*sin(20*pi*x).*cos(2*pi*x);
  case 23
    y=1 ./(1+(230*x-30).^2);
  otherwise
    y=zeros(1,legth(x));warning('Problem number is not between 1 and 23')
end

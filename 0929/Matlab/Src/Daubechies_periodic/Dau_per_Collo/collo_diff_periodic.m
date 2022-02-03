%function [derivative, error]= collo_diff_periodic(j,D,L,d)

%*******************************************************************
%Programme for differentiating a given function which is periodic on the 
%interval [0,1]
%dependencies
%collo_difmatrix_periodic.m
%********************************************************************
function [derivative, error] = collo_diff_periodic(j,D,L,d)
%j=5;
%D=10;L=1;
N=L*2^j;
 [hk,gk] = wfilters(['db' num2str(D/2)],'r');     % Low and high pass filters coefficients
  hk=hk';gk=gk';
x=(0:(1/N):(L-1/N))';
[D]=collo_difmatrix_periodic(D,j,j,0,1,d); 
syms t;
fun= input('enter the function you want to differntiate in variable t: ');
for i=1:size(x)
   test(i)=subs(fun,x(i));
end
test=test';
c=dst(test,hk);
nd_test=D*c;
derivative=nd_test;
fun_der=diff(fun,t,d);
for i=1:size(x)
   test_d(i)=subs(fun_der,x(i));
end

plot(x,nd_test);
hold on;
plot(x,test_d,'r');
error=norm(nd_test-test_d');

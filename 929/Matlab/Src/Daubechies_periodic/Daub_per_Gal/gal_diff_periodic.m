%function [derivative, error]= gal_diff_periodic(j,d,D,L)

%*******************************************************************
%Programme for differentiating a given function which is periodic on the 
%interval [0,L] using wavelet Galerkin method
% dependencies
%gal_difmatrix_periodic.m
%********************************************************************

%j=5;
%D=10;
function [derivative, error]= gal_diff_periodic(j,d,D,L)
N=L*2^j;
x=(0:(1/N):(L-1/N))';
D_1 = gal_difmatrix_periodic(d,N,L,D);
syms t;
fun= input('enter the function you want to differntiate in variable t: ');
for i=1:size(x)
   fun_n(i)=subs(fun,x(i));
end
fun_n=fun_n';
der_approx=D_1*fun_n;
derivative=der_approx;
fun_der=diff(fun,t,d);
for i=1:size(x)
   der_anal(i)=subs(fun_der,x(i));
end

plot(x,der_approx);
hold on;
plot(x,der_anal,'r');
error=norm(der_approx-der_anal',inf);


%*******************************************************************
%Programme for numerically verifying convergence results with respect to  resolution i.e. J, for wavelet
%galerkin method using daubechies wavelets in periodic case..
%Uses:
%difmatrix
%********************************************************************
clc;clf;
clear all;
D0=8; D1=14;
L=input('Enter the period of the function to be differentiated: ');
j0=2;J=5;
linestyles = cellstr(char('-',':','-.','--','-',':','-.','--','-',':','-',':',...
'-.','--','-',':','-.','--','-',':','-.'));
Markers=['o','x','+','*','s','d','v','^','<','>','p','h','.',...
'+','*','o','x','^','<','h','.','>','p','s','d','v',...
'o','x','+','*','s','d','v','^','<','>','p','h','.'];
d=input('Enter the order of differentiation: ');
syms t;
fun= input('Enter the function in variable t: ');
dfun_1=diff(fun,t,d);
%*****Convergence with respect to level of resolution (collocation points)*****************
for D=D0:2:D1
for j=j0:J
N=L*2^j;
x=(0:(1/N):(L-1/N))';
func=zeros(size(x));
for i=1:length(x)
   func(i)=subs(fun,x(i));
end
[D_1]= gal_difmatrix_periodic(d,N,L,D);
n_funcd_1=D_1*func;
a_funcd_1=zeros(size(x));
for i=1:length(x)
   a_funcd_1(i)=subs(dfun_1,x(i));
   end
error_df1_j(j-j0+1)=norm(n_funcd_1-a_funcd_1,inf);
j_vec(j-j0+1)=j;
end
semilogy(j_vec,error_df1_j,[linestyles{(D-D0+2)/2} Markers((D-D0+2)/2)],'linewidth',3,'markersize',3.5);
hold on;
end
xlabel('J');
p=['E^{(d)}(f,J)'];
ylabel({p});
set(gca,'FontSize',13);
legend('D=8','D=10','D=12','D=14');

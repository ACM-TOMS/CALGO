%*******************************************************************
%Programme for numerically verifying convergence results w.r.t wavelet genus D.
% in wavelet Galerkin method using daubechies wavelets in periodic case
%********************************************************************
clc;clf;
clear all;
D0=8; D1=12;
L=input('Enter the period of the function: ');
j=7;
N=L*2^j;
x=(0:(1/N):(L-1/N))';
func=zeros(size(x));
a_funcd_1=zeros(size(x));
linestyles = cellstr(char('-',':','-.','--','-',':','-.','--','-',':','-',':',...
'-.','--','-',':','-.','--','-',':','-.'));
Markers=['o','x','+','*','s','d','v','^','<','>','p','h','.',...
'+','*','o','x','^','<','h','.','>','p','s','d','v',...
'o','x','+','*','s','d','v','^','<','>','p','h','.'];
syms t;
n=input('Enter the order of differentiation: ');
fun= input('Enter the function in variable t: ');
for i=1:length(x)
   func(i)=subs(fun,x(i));
end
dfun_1=diff(fun,t,n);
for i=1:length(x)
   a_funcd_1(i)=subs(dfun_1,x(i));
   end

%*****Convergence with respect to wavelet genus D****************

for D=D0:2:D1
[D_1]=gal_difmatrix_periodic(n,N,L,D);
n_funcd_1=D_1*func;
error_df1_j((D-D0+2)/2)=norm(n_funcd_1-a_funcd_1,inf);
D_vec((D-D0+2)/2)=D;
end
set(gca,'FontSize',14);
semilogy(D_vec,error_df1_j,[linestyles{n} Markers(n)],'linewidth',3,'markersize',3.5);
xlabel('D');
p=['E^{(d)}(f,J)'];
ylabel({p});


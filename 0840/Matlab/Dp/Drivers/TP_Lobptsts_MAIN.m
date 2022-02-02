% TP_Lobptsts_MAIN.m
% input:      N = number of grid points
%                    c = bandwidth parameter c
% output:   xprol = vector of length(1,N) containing the prolate-Lobatto grid points
%                    wprol = vector of length(1,N) containing the prolate-Lobatto quadrature weights
%                    Deriv(i,j) = d C_{j}/dx(x_i)
%                   Deriv2(i,j) = d**2 C_{j}/dx**2 (x_i)

clf,clear

N=input('The number of interpolation points N=')
cstar=  (pi/2) * (( N-1) + 0.5)
c=input('bandwidth parameter c=')
tic
[xprol,wprol,Deriv,Deriv2d]=TP_LobptswtsPERT(N,c);
CPUtime = toc

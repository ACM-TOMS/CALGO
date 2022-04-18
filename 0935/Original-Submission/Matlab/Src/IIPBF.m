function [sol,reterr,evals]=IIPBF(fx,rho,tau,a,b,abserr,relerr,type)
% Entry point for IIPBF
% [sol,reterr,evals]=IIPBF(f,rhoI,tauI,a,b,abserr,relerr,typeI)

% Reference:
% Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J., and Lucas, 
% S. K. 2013. Algorithm XXX: IIPBF, a MATLAB toolbox for infinite integral 
% of products of two Bessel functions. To appear in ACM Transactions on 
% Mathematical Software. 

% Input
% abserr = required minimum absolute error
% relerr = required minimum relative error
% a = non negative integer
% b = non negative integer
% rho = positive real number
% tau = positive real number
% type = 'JJ', 'JY', or 'YY' (refer to accompanying paper)
%
% Output
% sol = estimated solution to integral
% reterr = predicted error
% evals = number of function evaluations
%
% Author: Sirong Zhang
% Last modified: 10/01/2010 by Jung Hun Kim
% - minor modifications (renaming of variables, function calls etc..)
% - resolve naming conflict of epsilon algorithm
% - parameter check

funcs = Build_funcs(fx, type, a, b, rho, tau);
if abs((rho/tau)-1) < 1e-10
    [zerh1_1,zerh1_2,its1]=first_2_zeros(funcs{4},'1',a,b,rho,tau);its2=0;
    nabserr=abserr/3.0; nrelerr=relerr/3.0;
    [part1,err1,neval1]=dqagea(funcs{1},0,zerh1_1,nabserr,nrelerr);
    [part2,err2,neval2]=dqagiea(funcs{3},zerh1_1,nabserr,nrelerr);
    [part3,err3,neval3]=mWtrans(zerh1_1,zerh1_2,funcs{2},nabserr,nrelerr,a,b,rho,tau,type);
    part4=0; err4=0; neval4=0;
elseif (abs(rho/tau) >1e+2 || abs(tau/rho) >1e+2)
    [zerh1_1,zerh1_2,its1]=first_2_zeros(funcs{4},'1',a,b,rho,tau); its2=0;
    nabserr=abserr/3.0; nrelerr=relerr/3.0;
    [part1,err1,neval1]=dqagea(funcs{1},0,zerh1_1,nabserr,nrelerr);
    [part2,err2,neval2]=epsalg(zerh1_1,zerh1_2,funcs{2},nabserr,nrelerr,rho,tau);
    [part3,err3,neval3]=epsalg(zerh1_1,zerh1_2,funcs{3},nabserr,nrelerr,rho,tau);
    part4=0;
    err4=0;
    neval4=0;
else
    [zerh1_1,zerh1_2,its1]=first_2_zeros(funcs{4},'1',a,b,rho,tau);
    [zerh2_1,zerh2_2,its2]=first_2_zeros(funcs{5},'2',a,b,rho,tau);
    nabserr=abserr/4.0;nrelerr=relerr/4.0;
    [part2,err2,neval2]=mWtrans(zerh1_1,zerh1_2,funcs{2},nabserr,nrelerr,a,b,rho,tau,type);
    [part3,err3,neval3]=epsalg(zerh2_1,zerh2_2,funcs{3},nabserr,nrelerr,rho,tau);
    if (zerh1_1<zerh2_1)
        [part1,err1,neval1]=dqagea(funcs{1},0,zerh1_1,nabserr,nrelerr);
        [part4,err4,neval4]=dqagea(funcs{3},zerh1_1,zerh2_1,nabserr,nrelerr);
    elseif (zerh1_1>=zerh2_1)
        [part1,err1,neval1]=dqagea(funcs{1},0,zerh2_1,nabserr,nrelerr);
        [part4,err4,neval4]=dqagea(funcs{2},zerh2_1,zerh1_1,nabserr,nrelerr);
    end
end

evals=its1+its2+neval1+neval2+neval3+neval4;
sol=part1+part2+part3+part4;

reterr = err1 + err2 + err3 + err4;

%% Plot first two zeros of h1
% saves the h1 function and values of the first two zeros
% (see h1_plotsNACONF2013.m to plot h1)

% zero1 = zerh1_1; zero2 = zerh1_2;
% h1 = funcs{4};
% save('h1plots.mat','h1','zero1','zero2')





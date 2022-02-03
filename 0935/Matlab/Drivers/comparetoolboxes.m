%%%% Comparison Experiments
% Compare Toolboxes  - script runs 10 test cases comparing performances of IIPBF vs.
% Besselint. See reference paper for details.
%   Author: J.T. Ratnanather & Jung Hun Kim
% 
%   Reference: Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J., 
%              and Lucas, S. K. 2012. IIPBF: a MATLAB toolbox for infinite 
%              integrals of product of Bessel functions. Pending review
%              from ACM Transactions on Mathematical Software
% 
%   Revision Date: 09/25/2012 Geoffrey Gunter, gmg@cis.jhu.edu

addpath('..','./BESSELINT')

clear all
close all

format long e

type = 'JJ';

nu = 1;
u = 2;
rho = 1;
tau = 2;
m = 1;


%Case 1

clear all; close all;
format long e;
type = 'JJ';nu = 1;u = 2;rho = 1;tau = 2;m = 1;

%I1 = 2/3;

tstart = tic;
[f1,err1] = besselint([1 1.5],[0 1],0,1.e-13,0);
telapsed11 = toc(tstart);
tstart = tic;
[res1,rel1,evals1]=IIPBF(@(x)1,1,1.5,0,1,1.e-14,1.e-13,type);
telapsed12 = toc(tstart);

err1_1 = err1(1)/(10.^floor(log10(err1(1))));
err1_2 = floor(log10(err1(1)));
rel1_1 = rel1/(10.^floor(log10(rel1)));
rel1_2 = floor(log10(rel1));
diff1 = abs(res1-f1);
diff1_1 = diff1/(10.^floor(log10(diff1)));
diff1_2 = floor(log10(diff1));
if (diff1==0)
    diff1_1 = 0;
    diff1_2 = 0;
end

fprintf('1 & %1.4f & %1.4f & %1.2f %s & %1.2f %s & %1.2f %s \n', telapsed11, telapsed12, err1_1, ['x$10^{', num2str(err1_2),'}$'], rel1_1, ['x$10^{', num2str(rel1_2),'}$'], diff1_1, ['x$10^{', num2str(diff1_2),'}$ \\']);

%Case 2

clear all; close all;
format long e;
type = 'JJ';nu = 1;u = 2;rho = 1;tau = 2;m = 1;

%I2 = 27/4096;
tstart = tic;
[f2,err2] = besselint([1 2],[0 5],-4);
telapsed21 = toc(tstart);
tstart = tic;
[res2,rel2,evals2]=IIPBF(@(x)1./x.^4,1,2,0,5,1.e-14,1.e-13,type);
telapsed22 = toc(tstart);

err2_1 = err2(1)/(10.^floor(log10(err2(1))));
err2_2 = floor(log10(err2(1)));
rel2_1 = rel2/(10.^floor(log10(rel2)));
rel2_2 = floor(log10(rel2));
diff2 = abs(res2-f2);
diff2_1 = diff2/(10.^floor(log10(diff2)));
diff2_2 = floor(log10(diff2));
if (diff2==0)
    diff2_1 = 0;
    diff2_2 = 0;
end

fprintf('2 & %1.4f & %1.4f & %1.2f %s & %1.2f %s & %1.2f %s \n', telapsed21, telapsed22, err2_1, ['x$10^{', num2str(err2_2),'}$'], rel2_1, ['x$10^{', num2str(rel2_2),'}$'], diff2_1, ['x$10^{', num2str(diff2_2),'}$ \\']);

%Case 3

clear all; close all;
format long e;
type = 'JJ';nu = 1;u = 2;rho = 1;tau = 2;m = 1;

%I3 = -6.050747903049e-3;
tstart = tic;
[f3,err3] = besselintr([1 1.1],[0 20],1,1);
telapsed31 = toc(tstart);
tstart = tic;
[res3,rel3,evals3]=IIPBF(@(x)x./(1+x.^2),1,1.1,0,20,1.e-14,1.e-13,type);
telapsed32 = toc(tstart);

err3_1 = err3(1)/(10.^floor(log10(err3(1))));
err3_2 = floor(log10(err3(1)));
rel3_1 = rel3/(10.^floor(log10(rel3)));
rel3_2 = floor(log10(rel3));
diff3 = abs(res3-f3);
diff3_1 = diff3/(10.^floor(log10(diff3)));
diff3_2 = floor(log10(diff3));
if (diff3==0)
    diff3_1 = 0;
    diff3_2 = 0;
end

fprintf('3 & %1.4f & %1.4f & %1.2f %s & %1.2f %s & %1.2f %s \n', telapsed31, telapsed32, err3_1, ['x$10^{', num2str(err3_2),'}$'], rel3_1, ['x$10^{', num2str(rel3_2),'}$'], diff3_1, ['x$10^{', num2str(diff3_2),'}$ \\']);


%Case 4

clear all; close all;
format long e;
type = 'JJ';nu = 1;u = 2;rho = 1;tau = 2;m = 1;

%I4 = 2/pi;
tstart = tic;
[f4,err4] = besselint([1 1],[0 1],-1);
telapsed41 = toc(tstart);
tstart = tic;
[res4,rel4,evals4]=IIPBF(@(x)1./x,1,1,0,1,1.e-14,1.e-13,type);
telapsed42 = toc(tstart);

err4_1 = err4(1)/(10.^floor(log10(err4(1))));
err4_2 = floor(log10(err4(1)));
rel4_1 = rel4/(10.^floor(log10(rel4)));
rel4_2 = floor(log10(rel4));
diff4 = abs(res4-f4);
diff4_1 = diff4/(10.^floor(log10(diff4)));
diff4_2 = floor(log10(diff4));
if (diff4==0)
    diff4_1 = 0;
    diff4_2 = 0;
end

fprintf('4 & %1.4f & %1.4f & %1.2f %s & %1.2f %s & %1.2f %s \n', telapsed41, telapsed42, err4_1, ['x$10^{', num2str(err4_2),'}$'], rel4_1, ['x$10^{', num2str(rel4_2),'}$'], diff4_1, ['x$10^{', num2str(diff4_2),'}$ \\']);

%Case 5

clear all; close all;
format long e;
type = 'JJ';nu = 1;u = 2;rho = 1;tau = 2;m = 1;

%I5 = 4/(3*pi);
tstart = tic;
[f5,err5] = besselint([1 1],[1 1],-2);
telapsed51 = toc(tstart);
tstart = tic;
[res5,rel5,evals5]=IIPBF(@(x)1./x.^2,1,1,1,1,1.e-14,1.e-13,type);
telapsed52 = toc(tstart);

err5_1 = err5(1)/(10.^floor(log10(err5(1))));
err5_2 = floor(log10(err5(1)));
rel5_1 = rel5/(10.^floor(log10(rel5)));
rel5_2 = floor(log10(rel5));
diff5 = abs(res5-f5);
diff5_1 = diff5/(10.^floor(log10(diff5)));
diff5_2 = floor(log10(diff5));
if (diff5==0)
    diff5_1 = 0;
    diff5_2 = 0;
end

fprintf('5 & %1.4f & %1.4f & %1.2f %s & %1.2f %s & %1.2f %s \n', telapsed51, telapsed52, err5_1, ['x$10^{', num2str(err5_2),'}$'], rel5_1, ['x$10^{', num2str(rel5_2),'}$'], diff5_1, ['x$10^{', num2str(diff5_2),'}$ \\']);

%Case 19

clear all; close all;
format long e;
type = 'JJ';nu = 1;u = 2;rho = 1;tau = 2;m = 1;
[K,E] = ellipke(rho.^2/(rho^2+u^2));

%I19 = (K-E)/(2*pi*rho*sqrt(rho^2+u^2));

tstart = tic;
[f19,err19] = besselintc([rho rho],[0 1],1,2*u);
telapsed191 = toc(tstart);
tstart = tic;
[res19,rel19,evals19]=IIPBF(@(x)x.*exp(-2*u.*x),1,1,0,1,1.e-14,1.e-13,type);
telapsed192 = toc(tstart);

err19_1 = err19(1)/(10.^floor(log10(err19(1))));
err19_2 = floor(log10(err19(1)));
rel19_1 = rel19/(10.^floor(log10(rel19)));
rel19_2 = floor(log10(rel19));
diff19 = abs(res19-f19);
diff19_1 = diff19/(10.^floor(log10(diff19)));
diff19_2 = floor(log10(diff19));
if (diff19==0)
    diff19_1 = 0;
    diff19_2 = 0;
end

fprintf('19 & %1.4f & %1.4f & %1.2f %s & %1.2f %s & %1.2f %s \n', telapsed191, telapsed192, err19_1, ['x$10^{', num2str(err19_2),'}$'], rel19_1, ['x$10^{', num2str(rel19_2),'}$'], diff19_1, ['x$10^{', num2str(diff19_2),'}$ \\']);

%Case 20

clear all; close all;
format long e;
type = 'JJ';nu = 1;u = 2;rho = 1;tau = 2;m = 1;
[K,E] = ellipke(rho.^2/(rho^2+u^2));

%I20 = K/(pi*sqrt(rho^2+u^2));
tstart = tic;
[f20,err20] = besselintc([rho rho],[0 0],0,2*u);
telapsed201 = toc(tstart);
tstart = tic;
[res20,rel20,evals20]=IIPBF(@(x)exp(-2*u.*x),1,1,0,0,1.e-14,1.e-13,type);
telapsed202 = toc(tstart);

err20_1 = err20(1)/(10.^floor(log10(err20(1))));
err20_2 = floor(log10(err20(1)));
rel20_1 = rel20/(10.^floor(log10(rel20)));
rel20_2 = floor(log10(rel20));
diff20 = abs(res20-f20);
diff20_1 = diff20/(10.^floor(log10(diff20)));
diff20_2 = floor(log10(diff20));
if (diff20==0)
    diff20_1 = 0;
    diff20_2 = 0;
end

fprintf('20 & %1.4f & %1.4f & %1.2f %s & %1.2f %s & %1.2f %s \n', telapsed201, telapsed202, err20_1, ['x$10^{', num2str(err20_2),'}$'], rel20_1, ['x$10^{', num2str(rel20_2),'}$'], diff20_1, ['x$10^{', num2str(diff20_2),'}$ \\']);

%Case 21
%% not true if (c==a) 

clear all; close all;
format long e;
type = 'JJ';nu = 1;u = 2;rho = 1;tau = 2;m = 1;
[K,E] = ellipke(rho.^2/(rho^2+u^2));

%I21 = ((2*u^2+rho^2)*K-2*(rho^2+u^2)*E)/(pi*rho^2*sqrt(rho^2+u^2));

tstart = tic;
[f21,err21] = besselintc([rho rho],[1 1],0,2*u);
telapsed211 = toc(tstart);
tstart = tic;
[res21,rel21,evals19]=IIPBF(@(x)exp(-2*u.*x),1,1,1,1,1.e-14,1.e-13,type);
telapsed212 = toc(tstart);

err21_1 = err21(1)/(10.^floor(log10(err21(1))));
err21_2 = floor(log10(err21(1)));
rel21_1 = rel21/(10.^floor(log10(rel21)));
rel21_2 = floor(log10(rel21));
diff21 = abs(res21-f21);
diff21_1 = diff21/(10.^floor(log10(diff21)));
diff21_2 = floor(log10(diff21));
if (diff21==0)
    diff21_1 = 0;
    diff21_2 = 0;
end

fprintf('21 & %1.4f & %1.4f & %1.2f %s & %1.2f %s & %1.2f %s \n', telapsed211, telapsed212, err21_1, ['x$10^{', num2str(err21_2),'}$'], rel21_1, ['x$10^{', num2str(rel21_2),'}$'], diff21_1, ['x$10^{', num2str(diff21_2),'}$ \\']);

%Case 22

clear all; close all;
format long e;
type = 'JJ';nu = 1;u = 2;rho = 1;tau = 2;m = 1;

%I22 = 0.75*quadgk(@(z)cos(z).^2./(u^2+cos(z).^2).^2.5,0,pi/2)/pi;

tstart = tic;
[f22,err22] = besselintc([rho rho],[1 1],2,2*u);
telapsed221 = toc(tstart);
tstart = tic;
[res22,rel22,evals22]=IIPBF(@(x)x.^2.*exp(-2*u.*x),1,1,1,1,1.e-14,1.e-13,type);
telapsed222 = toc(tstart);

err22_1 = err22(1)/(10.^floor(log10(err22(1))));
err22_2 = floor(log10(err22(1)));
rel22_1 = rel22/(10.^floor(log10(rel22)));
rel22_2 = floor(log10(rel22));
diff22 = abs(res22-f22);
diff22_1 = diff22/(10.^floor(log10(diff22)));
diff22_2 = floor(log10(diff22));
if (diff22==0)
    diff22_1 = 0;
    diff22_2 = 0;
end

fprintf('22 & %1.4f & %1.4f & %1.2f %s & %1.2f %s & %1.2f %s \n', telapsed221, telapsed222, err22_1, ['x$10^{', num2str(err22_2),'}$'], rel22_1, ['x$10^{', num2str(rel22_2),'}$'], diff22_1, ['x$10^{', num2str(diff22_2),'}$ \\']);

%Case 23

clear all; close all;
format long e;
type = 'JJ';nu = 1;mu=1;u = 2;rho = 1;tau = 2;m = 1;

%I23 = u^(mu-nu)*besseli(nu, rho*u)*besselk(mu, tau*u)

tstart = tic;
[f23, err23] = besselintr([rho,tau], [mu nu], m, u);
telapsed231 = toc(tstart);
tstart = tic;
[res23,rel23,evals23]=IIPBF(@(x)x.^(nu-mu+1)./(u^2+x.^2),rho,tau,mu,nu,1.e-14,1.e-13,type);
telapsed232 = toc(tstart);

err23_1 = err23(1)/(10.^floor(log10(err23(1))));
err23_2 = floor(log10(err23(1)));
rel23_1 = rel23/(10.^floor(log10(rel23)));
rel23_2 = floor(log10(rel23));
diff23 = abs(res23-f23);
diff23_1 = diff23/(10.^floor(log10(diff23)));
diff23_2 = floor(log10(diff23));
if (diff23==0)
    diff23_1 = 0;
    diff23_2 = 0;
end

fprintf('23 & %1.4f & %1.4f & %1.2f %s & %1.2f %s & %1.2f %s \n', telapsed231, telapsed232, err23_1, ['x$10^{', num2str(err23_2) '}$'], rel23_1, ['x$10^{', num2str(rel23_2) '}$'], diff23_1, ['x$10^{', num2str(diff23_2),'}$ \\']);



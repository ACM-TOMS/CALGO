% Example script to test AMGKQ in one dimension and many dimensions.
% Copyright (C) 2014 Robert W. Johnson
% Based on TEST_INT (C) 2009 & 2011 John Burkardt
%
% This file is part of AMGKQ.  See AMGKQ.M for details.
% Copyright (C) 2014 Robert W. Johnson, Alphawave Research.
% This is free software see GPLv3 or greater for copying conditions.
% There is ABSOLUTELY NO WARRANTY not even for MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  For details, see GPLv3 or greater.

% preparations
close all
clear all
tabcell = cell(1,9);
ntab = 0;
disp(' ')

fx = @(X) exp(X)
a = 0
b = 1
exact = exp(1)-1
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) 1./(1+X.^4)
a = 0
b = 1
exact = (log(sqrt(2)+2)-log(2-sqrt(2))+pi)/2^(5/2)
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) 1./(1+exp(X))
a = 0
b = 1
exact = -log(exp(1)+1)+log(2)+1
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) X./(exp(X)-1)
a = 0
b = 1
exact = 7.77504634112247572375054e-1
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) X./(exp(X)+1)
a = 0
b = 1
exact = 1.705573495024381713847106e-1
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) 0.92*cosh(X)-cos(X)
a = -1
b = 1
exact = (2*(23*sinh(1)-25*sin(1)))/25
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) exp(X).*cos(X)
a = 0
b = pi
exact = -(exp(pi)+1) / 2
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) 1./(1+X.^2+X.^4)
a = -1
b = 1
exact = (3*log(3)+sqrt(3)*pi)/6
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) 50./pi./(2500*X.^2+1)
a = 0
b = 1
exact = atan(50)/pi
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) sqrt(X)
a = 0
b = 1
exact = 2 / 3
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) sqrt(50)*exp(-50*pi*X.^2)
a = 0
b = 10
exact = erf(25*2^(3/2)*sqrt(pi))/2
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) 25*exp(-25*X)
a = 0
b = 10
exact = 1-exp(-250)
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) 1./sqrt(X)
a = 0
b = 1
exact = 2
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) log(X)
a = 0
b = 1
exact =-1
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) sqrt(abs(X+0.5))
a = -1
b = 1
exact = sqrt(3)/sqrt(2)+1/(3*sqrt(2))
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) log(abs(X-0.7))
a = 0
b = 1
exact = (7*log(7/10)-7)/10+(3*log(3/10)-3)/10
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) 2./(2+sin(10*pi*X))
a = 0
b = 1
exact = 2/sqrt(3)
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) (sin(50*pi*X)).^2
a = 0
b = 1
exact = 1/2
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) exp(cos(X))
a = 0
b = 2*pi
exact = 2*pi*besseli(0,1)
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) 1./(X.^(1/2)+X.^(1/3))
a = 0
b = 1
exact = 5-6*log(2)
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) exp(- X).*sin(50*X)
a = 0
b = 2*pi
exact =(50 / 2501) *(1-exp(- 2*pi))
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) (X <= exp(1)-2)./(X+2)
a = 0
b = 1
exact = 1-log(2)
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) 1./(1+X.^2)
a = -4
b = 4
exact = 2*atan(4)
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) sqrt(-log(X))
a = 0
b = 1
exact = sqrt(pi) / 2
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) (10*X-1).*(10*X-1.1).*(10*X-1.2).*(10*X-1.3)
a = 0
b = 1
exact = 1627879 / 1500
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) log(X).*sqrt(X)
a = 0
b = 1
exact = -4 / 9
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) log(X)./sqrt(X)
a = 0
b = 1
exact = -4
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) (0.3 <= X)
a = 0
b = 1
exact = 7 / 10
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) (sech(10 *(X-0.2))).^2+(sech(100 *(X-0.4))).^4+(sech(1000 *(X-0.6))).^6
%k = [1:3]'
%fx = @(X) sum((sech(10.^k .*(X-k/5))).^(2*k),1)
a = 0
b = 1
exact = 2.108027355005492875505979e-1
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

k=[1:40]';
fx = @(X) sum(bsxfun(@rdivide, cos (7.^k*X*pi/2), 2.^k), 1)
%fx = @(X) sum(cos (7.^k*X*pi/2)./2.^k, 1)
a = 0
b = 1
exact = sum(2*sin(7.^k*pi/2)./(2*7).^k/pi)
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(X) (1./X).*sin(1./X)
a = 0
b = 1
% exact = pi / 2 - sinint(1)
exact = 6.24713256427713581331318e-1
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

% redo difficult integrands
numfx = ntab
redofx = find(abs(cell2mat(tabcell(:,9))) > 1.5e-8)
numredo = length(redofx)
for kfx = redofx'
    [fx, a, b, exact] = deal(tabcell{kfx,1:4})
    [RES, ERR, NSUB, FL] = amgkq(fx,a,b,[],0,1e3), ACC = RES-exact, disp(' ')
    ntab = ntab + 1;
    tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};
end

% format table for display
fnos = [1:ntab]';
fnos = fnos - floor((fnos - 1)/numfx)*numredo;
headstr = 'No.       A        B    EXACT      RES      ERR     NSUB       FL      ACC';
tabshow = strvcat(headstr,num2str([fnos,cell2mat(tabcell(:,2:end))],'%8.2g '))

% begin ND test
fprintf('\n Multivariate test might take awhile...proceed?  Ctrl-C to abort.\n')
fprintf(' Reduce maximum number of dimensions in script if necessary.\n')
fprintf('Start in '), for ps = 1:9, fprintf('%i ... ',10-ps), pause(1), end
maxnd = 6;
bigtab = cell(maxnd,5);
for fmax = 1:maxnd
    clear -f amgkq
    funnum = [1:fmax]+0;
    numfun = length(funnum)
    bigfun = tabcell(funnum,1);
    bigfunf{fmax} = func2str(bigfun{fmax});
    bigfunchar = char(bigfunf)
    biga = cell2mat(tabcell(funnum,2))
    bigb = cell2mat(tabcell(funnum,3))
    bigexact = prod(cell2mat(tabcell(funnum,4)))
    bf{fmax} = ['bigfun{' num2str(fmax) '}(X(' num2str(fmax) ',:))'];
    if fmax > 1, bfrows = strjoin(bf,'.*'); else bfrows = bf{fmax}; end
    bfstr = ['bigF = @(X) ' bfrows]; eval(bfstr);
    longarg{8} = 2;     % set iteration verbosity
    t = cputime;
    [RES, ERR, NSUB, FL] = amgkq(bigF,biga,bigb,longarg{:}), 
    inttime(fmax) = cputime - t;
    ACC = RES-bigexact, disp(' ')
    [bigtab{fmax,:}] = deal(RES, ERR, NSUB, FL, ACC);
    %save -b temp.mat
end
%load temp.mat

% format table for display
for f = 1:fmax
    nG = 7; nD = f; nX = (2*nG+1)^nD; nDnX = nD*nX; ncnt(f,:) = [nD, nX, nDnX];
end
fnos = [1:ntab]';
fnos = fnos - floor((fnos - 1)/numfx)*numredo;
headstr = 'ND        NX       NDNX      RES      ERR     NSUB       FL      ACC';
tabshow = strvcat(headstr,num2str([ncnt, cell2mat(bigtab)],'%i %10i %10i %8.2g %8.2g %8i %8i %8.2g \n'))

% make figure
pappos = [1 1 6 2.5];
figure(1), clf
set(gcf,'position',pappos*90,'paperposition',pappos)
for f = 1:3
    sb{f} = subplot(1,3,f); set(sb{f},'linewidth',2,'units','normalized')
    switch f
    case 1, plot(1:fmax,log10(ncnt(:,2)),'s-','linewidth',5); ylim([1 8])
    case 2, plot(1:fmax,log10(ncnt(:,3)),'s-','linewidth',5); ylim([1 8])
    case 3, plot(1:fmax,log10(inttime),'s-','linewidth',5); ylim([-1 4]), set(sb{f},'ytick',-1:4)
    end
    text(.5,-.25,'N_D','units','normalized')
    ylbl = {'LOG_{10} N_X','LOG_{10} N_D N_X','LOG_{10} TIME'};
    text(-.35,.5,ylbl{f},'units','normalized','rotation',90,'horizontalalignment','center')
    text(1.05,1.0,['(' char(96+f) ')'],'units','normalized','verticalalignment','bottom')
end
for f = 1:3
    apos = get(sb{f},'position'); 
    apos(1) = apos(1) + (f-2)*.05 + .01;
    apos(3) = apos(3) - .01;
    apos(2) = apos(2) + .12;
    apos(4) = apos(4) - .15;
    set(sb{f},'position',apos,'xtick',1:fmax)
end


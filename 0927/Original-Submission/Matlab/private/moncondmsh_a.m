function [nptcond,cond_epsi,fatt_r1r2,fatt_r2] = moncondmsh_a(nmsh,omega,h,a,b,moncondmsh_par)
%
%   Private function for twpbvpc
%
%  
%
%       Authors:
%
%       Jeff R. Cash 
%            (Department of Mathematics, Imperial College,  London, England.)
%       Davy  Hollevoet 
%            (Vakgroep Toegepaste Wiskunde en Informatica, Universiteit Gent, Belgium.)
%       Francesca Mazzia  
%            (Dipartimento di Matematica, Universita' di Bari, Italy)
%       Abdelhameed Nagy Abdo
%            (Dipartimento di Matematica, Universit\`a di Bari, Italy)
%            (Dept. of Mathematics, Faculty of Sciences, Benha  University,Egypt)
%            
%
sfatt_alpha = moncondmsh_par.sfatt_alpha;


alfa2 =(sfatt_alpha/(b-a))*sum(abs(omega(1:nmsh-1)- omega(2:nmsh)));
epsi_gamma2 = (abs(omega(1:nmsh-1)- omega(2:nmsh))  + h.*alfa2);
cond_epsi = epsi_gamma2;


r1 = max(cond_epsi);

cond_epsi=cond_epsi/r1;

r2 = sum(cond_epsi)/(nmsh-1);
% to plot the monitor functions
% figure(1)
% semilogy(1:length(h),epsi_gamma1/(max(epsi_gamma1)),'r+',1:length(h),r2*1e-3,'-',1:length(h), epsi_gamma2/max(epsi_gamma2),'og')
% figure(2)
% semilogy(1:length(h),epsi_gamma1./h,'+',1:length(h), epsi_gamma2./h,'o')
% pause

sfatt_r2 = moncondmsh_par.sfatt_r2;
fatt_r2  = r2*sfatt_r2;
sfatt_r1r2= moncondmsh_par.sfatt_r1r2;
fatt_r1r2=max(sfatt_r1r2,r2) ;  % nptm symbole to intervals
nptm = sum(cond_epsi >= fatt_r1r2);


if (nptm <= 1)
    nptcond =  14;                       % nptcond symbole the number of adding points
elseif (nptm <= 2)
    nptcond =  10;
elseif (nptm <= 4)
    nptcond =  8;
elseif (nptm <= 8)
    nptcond = 6;
elseif (nptm <= (nmsh)/20 )
    nptcond = 4;
else
    nptcond = 2;
end


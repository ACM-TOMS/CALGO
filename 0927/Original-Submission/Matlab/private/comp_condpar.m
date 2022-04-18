function [omega,condpar,stabcond] = comp_condpar(A,m,nmsh,h,a,b,condpar_old,linear,compcond_par)
% function for the computation of the conditioning parameters
%
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
flmax = realmax;
fatt_nl_sigma = compcond_par.fatt_nl_sigma;
fatt_lin_sigma = compcond_par.fatt_lin_sigma;

nblk = nmsh;                %number of blocks
mnblk = m*(nblk);
omega = zeros(1,nblk);        % omega is row of dimenison (1*N+1) 
ind = 1:m;
bv = zeros(mnblk,m);           % rigth hand side of linear system of dimension (((N+1)*m)*m)
bv(1:m,1:m) = speye(m);
bomega = zeros(nblk,m);
yy = A.Q*(A.U\(A.L\(A.P*bv)));         % solution of linear system of dimension (((N+1)*m)*m)
for i = 1:nblk              
     omega(i) = norm(yy(ind,:) ,'inf');     % compute omega: block by block i.e., omega(1) is first block of yy and so on                 
     for j=1:m
         bomega(i,j) = max(abs(yy(ind,j))); 
     end
     ind = ind + m;
end
[~,indmaxom] = max( sum(abs( yy ),2));


gamma1 = (1/(b-a))* sum(h .* max([omega(2:nmsh);omega(1:nmsh-1)])); 
kappa1 = max(omega(1:nmsh));


bomega_max=zeros(nmsh-1,m);
for j=1:m
    bomega_max(:,j)  =  max( [bomega(2:nmsh,j)' ; bomega(1:nmsh-1,j)' ] )';
end
bgammai=zeros(1,m);
for j=1:m
    bgammai(j) = sum(h'.*bomega_max(:,j))/(b-a);
end


sigma = max(max(abs(bomega))./bgammai); 
% sigma = kappa1/gamma1;  %old value of sigma
% T is the initial value for the Higham algorithm
X0 = [ones(mnblk,1)/mnblk, zeros(mnblk,1)];
X0(indmaxom,2)=1;

D=kron(spdiags([1;h(:)],0,nmsh,nmsh),speye(m));
[kappa, ~, W]  = normest1(@condaux,2,X0,A.L,A.U,A.P,A.Q,D);             % compute the infinity norm of Matrix by using one norm '' Higham algorithm''
 
kappa2=norm(W(m+1:end),1);


kappaold=condpar_old.kappa;
kappa1old=condpar_old.kappa1;
gamma1old=condpar_old.gamma1;
sigmaold=condpar_old.sigma;


stab_kappa = (abs(kappaold-kappa)/(kappa)  < 5e-2)  & kappa < flmax;
stab_kappa1 = (abs(kappa1old-kappa1)/(kappa1) < 5e-2)  & kappa1 < flmax & gamma1 < flmax;
stab_gamma1 = (abs(gamma1old-gamma1)/(gamma1) < 5e-2)  & gamma1 < flmax & kappa1 <  flmax;
stab_all = (stab_kappa & stab_kappa1 & stab_gamma1);
if linear
   stiff_cond =  ( sigma >=  fatt_lin_sigma  );
else
   stiff_cond =  ( sigma >=  fatt_nl_sigma  ); 
end
ill_cond   =  ( kappa >  1e14  ) ;

condpar.kappa = kappa;
condpar.kappa1 = kappa1;
condpar.kappa2 = kappa2;
condpar.gamma1 = gamma1;
condpar.sigma = sigma;

 
stabcond.kappa = stab_kappa;
stabcond.kappa1 = stab_kappa1;
stabcond.gamma1 = stab_gamma1;
stabcond.all = stab_all;
stabcond.stiff = stiff_cond;
stabcond.ill = ill_cond;


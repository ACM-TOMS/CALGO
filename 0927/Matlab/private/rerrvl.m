function [rerr,remax,itlmx]=rerrvl(u,usvrex,ltol,tol,adjerr)
%
%   Private function for twpbvpc
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
    ntol=size(tol,1);
    nmold=size(usvrex,2);
    
    itlmx=1;
    remax=0;
    rerr=zeros(size(usvrex));
    
    for it=1:ntol
        icmp=ltol(it);
        rerr(icmp,:)=abs(usvrex(icmp,:)-u(icmp,1:2:end));
        denom=abs(usvrex(icmp,:));
        denom(denom<1)=1;
        rerel=rerr(icmp,:)./denom;
        rerelmax=max(rerel);
        if rerelmax>remax
            remax=rerelmax;
            itlmx=it;
        end
    end
    
    if adjerr
        rerr=max(rerr(ltol,1:end-1),rerr(ltol,2:end));
    end
    
end
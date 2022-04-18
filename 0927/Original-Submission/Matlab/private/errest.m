function [errok,errmax]=errest(ltol,tol,etest,y,yold)
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
    errok=1;
    errmax=0;
    for it=1:size(tol,1)
        icmp=ltol(it);

        er=y(icmp,:)-yold(icmp,:);
        denum=yold(icmp,:);
        denum(abs(denum)<1)=1;
        errel=abs(er./(tol(it)*denum));
        errmax = max(errmax,  max(errel));
        errok=errok & all(errel<=etest(it));
    end
    
end
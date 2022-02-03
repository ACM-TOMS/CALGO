function [defexp,problem]=dfexcl_m(problem,h,t,y,fty)
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
    hb=repmat(h,problem.ncomp,1);
    % todo: pass hb around?
    
    yn=y(:,1:end-1);
    ynpo=y(:,2:end);
    fn=fty(:,1:end-1);
    fnpo=fty(:,2:end);
    
    tn=t(1,1:end-1);
    tnpof=tn+h./4;
    tnpoh=tn+h./2;
    tnptf=tn+3*h./4;
    
    tmp1=(5*ynpo+27*yn)./32+ hb.*(9*fn-3*fnpo)./64;
    tmp2=(27*ynpo+5*yn)./32+ hb.*(3*fn-9*fnpo)./64;
    tmp3=problem.f(tnpof,tmp1);
    tmp4=problem.f(tnptf,tmp2);
     
    tmp1=(ynpo+yn)./2 + 5*hb.*(fnpo-fn)./24 + 2*hb.*(tmp3-tmp4)./3;
    tmp2=problem.f(tnpoh,tmp1);
    
    defexp=(yn-ynpo)+hb.*(7*(fnpo+fn)./90+16*(tmp3+tmp4)./45+2*tmp2./15);
    if ~problem.vectorized 
         problem.NFUN = problem.NFUN + length(tnpof)+ length(tnptf) + length(tnpoh);
       else
         problem.NFUN = problem.NFUN + 3;
    end
    
end
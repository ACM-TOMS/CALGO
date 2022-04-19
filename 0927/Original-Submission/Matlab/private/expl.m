function [defexp,problem] = expl(problem,h,t,y,fty)
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
    yn=y(:,1:end-1);
    ynpo=y(:,2:end);
    fn=fty(:,1:end-1);
    fnpo=fty(:,2:end);
 
    t1 = (5*ynpo + 27*yn)./32+ hb.*(9*fn - 3*fnpo)./64;
    t2 = (27*ynpo + 5*yn)./32+ hb.*( 3*fn -9*fnpo )./64;
    
    tn=t(1,1:end-1);
    tnpof=tn+h./4;
    tnptf=tn+3*h./4;
    tnpoh=tn+h./2;
    
    t3=problem.f(tnpof,t1);
    t4=problem.f(tnptf,t2);

    t1 = (ynpo + yn)./2 + 5*hb.*(fnpo - fn)./24 - 2*hb.*(t4 - t3)./3;
    t2 = problem.f(tnpoh,t1);
    defexp = hb.*(7*(fnpo + fn)./90  + 16*(t3 + t4)./45 + 2*t2./15) - ynpo + yn;
    if ~problem.vectorized 
         problem.NFUN = problem.NFUN + length(tnpof)+ length(tnptf) + length(tnpoh);
       else
         problem.NFUN = problem.NFUN + 3;
    end
  end

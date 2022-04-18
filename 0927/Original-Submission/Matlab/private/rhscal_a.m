function [rhs,rnsq,fty,problem]=rhscal_a(problem,h,t,y,lambda,dc)
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
    nmsh=size(y,2);
    
    yn=y(:,1:end-1);
    ynpo=y(:,2:end);
    tn=t(1,1:end-1);
    tnpo=t(1,2:end);
    fty=problem.f(t,y,lambda);
    fn=fty(:,1:end-1);
    fnpo=fty(:,2:end); 
    
    hb=repmat(h,problem.ncomp,1);
    ynh=(yn+ynpo)./2-hb.*(fnpo-fn)./8;
    tnh=(tnpo+tn)./2;
    fnpoh=problem.f(tnh,ynh,lambda);

    nfrows=problem.ncomp*(nmsh-1);
    rhsi=reshape(yn-ynpo+dc+hb.*(fn+4.*fnpoh+fnpo)./6,nfrows,1);
    rhsbc=-problem.g(y(:,1),y(:,end),lambda);
    rhs=cat(1,rhsbc,rhsi);
    
    rnsq=sum(rhs.^2);
    if ~problem.vectorized 
     problem.NFUN = problem.NFUN + length(t)+ length(tnh);
  else
     problem.NFUN = problem.NFUN + 2;
    end
 problem.NBC = problem.NBC + 1;

    %if problem.debug, disp(sprintf('rnsq %g',rnsq)); end
    %if problem.debug, disp(rhs'); end
end
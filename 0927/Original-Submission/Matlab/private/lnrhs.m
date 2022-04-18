function  [rhs,f,problem] = lnrhs(problem,h,t,nmsh,y)
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
f = problem.f(t,y);
h = t(2:nmsh) - t(1:nmsh-1);
hb=repmat(h,problem.ncomp,1);
uint = (1/2)*(y(:,1:nmsh-1) + y(:,2:nmsh)) - (1/8)* hb.*(f(:,2:nmsh) - f(:,1:nmsh-1));
thalf = (t(1:nmsh-1) + t(2:nmsh))/2;
fhalf = problem.f(thalf,uint);
rhsi(:,2:nmsh) = - y(:,2:nmsh) + y(:,1:nmsh-1) + (1/6)*hb.*(f(:,1:nmsh-1) + f(:,2:nmsh) + 4*fhalf(:,1:nmsh-1) );
rhsi(:,1) = -problem.g(y(:,1),y(:,end));
rhs = reshape(rhsi,problem.ncomp*nmsh,1);
if ~problem.vectorized 
     problem.NFUN = problem.NFUN + length(t)+ length(thalf);
  else
     problem.NFUN = problem.NFUN + 2;
end

problem.NBC = problem.NBC + 1;

return
end
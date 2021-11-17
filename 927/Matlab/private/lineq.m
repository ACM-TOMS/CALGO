function [y,iflag,problem,rhs,fty,J,bhold,chold] = lineq(problem,h,t,nmsh,y,dc,ludone,J)
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

iflag = 0;

ncomp=problem.ncomp;
if (~ludone) 
  [rhs,fty,problem] = lnrhs(problem,h,t,nmsh,y);            %rigth hand side
  [J,bhold,chold]=jaccal(problem,h,t,y,fty);        % matrix 
  [dy,J]=solvesls(J,rhs);                           % solution
  ludone = 1; 
else
    delu = reshape([zeros(ncomp,1) dc],ncomp*nmsh,1);
   [dy,J]=solvesls(J,delu); 
   
end 
if checksingular()
   iflag=-1;
   if problem.debug, disp('singular Jacobian'); end
 
end

delu = reshape(dy,ncomp,nmsh);
y=y+delu;

return

end
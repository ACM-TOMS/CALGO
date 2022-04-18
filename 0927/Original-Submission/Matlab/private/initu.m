function y=initu(problem,t,ncomp)
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
    if (problem.debug) fprintf('initu /n') ; end
    if isscalar(problem.y0)
    	y=repmat(problem.y0,ncomp,size(t,2));
    elseif size(problem.y0,2)==1
        y=repmat(problem.y0,1,size(t,2));
    else
        y=interp1q(problem.x0',problem.y0',t')';
    end
end
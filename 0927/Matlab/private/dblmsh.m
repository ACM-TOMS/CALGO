function [t,h,nmsh,maxmsh,told]=dblmsh(problem,t,nmax)

    % doubles the number of intervals in mesh t if the new number of
    % meshpoints does not exceed nmax
    %
    % returns t: new meshpoints
    %         h: new intermesh distances
    %         nmsh: new number of meshpoints
    %         maxmsh: if nmax would be exceeded by this action
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
    
    told=t;
    nmsh=size(t,2);
    
    if 2*nmsh-1>nmax
        maxmsh=1;
    else
        maxmsh=0;
        nmsh=2*nmsh-1;
        
        tp=t(1:end-1);
        tn=(tp+t(2:end))./2;
        t=cat(2,reshape(cat(1,tp,tn),1,[]),t(end));
    end
    
    h=t(2:end)-t(1:end-1);
    if ~maxmsh
       if problem.debug, disp(sprintf('dblmsh: new mesh has %g points',nmsh)); end
    else
       if problem.debug, disp(sprintf('dblmsh: maximum mesh exceeded, nmnew=%g',2*nmsh-1)); end
    end

end
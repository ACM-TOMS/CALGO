function [t,h,nmsh,maxmsh]=smpmsh(problem,t,intref,numadd,nmax)
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
    % inserts points into mesh t
    % if intref=1, numadd points are inserted in the first interval
    % if intref=nmsh-1, numadd points are inserted in the last interval
    % else numadd points are added in interval [intref-1,intref,intref+1]
    % only if the new number of meshpoints does not exceed nmax
    %
    % returns t: new meshpoints
    %         h: new intermesh distances
    %         nmsh: new number of meshpoints
    %         maxmsh: if nmax would be exceeded by this action

    maxmsh=0;
    nmsh=size(t,2);
    numadd=min(49,max(numadd,4));
    
    if problem.debug, disp(sprintf('nmsh %g, intref %g, numadd %g',nmsh,intref,numadd)); end
    
    if intref==1
        if nmsh+numadd>nmax
            maxmsh=1;
        else
            nmsh=nmsh+numadd;
            t=cat(2,linspace(t(1),t(2),numadd+2),t(3:end));
        end
    elseif intref==nmsh-1
        if nmsh+numadd>nmax
            maxmsh=1;
            nmsh=nmsh+numadd;
        else
            t=cat(2,t(1:end-2),linspace(t(end-1),t(end),numadd+2));
            nmsh=nmsh+numadd;
        end
    else
        if nmsh+3*numadd>nmax
            maxmsh=1;
        else
            numadd=min(numadd,9);
            nmsh=nmsh+numadd*3;
            t1=linspace(t(intref-1),t(intref),numadd+2);
            t2=linspace(t(intref),t(intref+1),numadd+2);
            t3=linspace(t(intref+1),t(intref+2),numadd+2);
            t=[t(1:intref-2) t1(1:end-1) t2(1:end-1) t3 t(intref+3:end)];
        end
    end
    h=t(2:end)-t(1:end-1);
    
    if problem.debug, disp(sprintf('smpmsh: new mesh: %g points',nmsh)); end
end
function ynew=interpu(told,yold,tnew)
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
    % interpolates yold in points tnew for every component
    %
    % returns ynew: interpolated values

    ncomp=size(yold,1);
    ynew=zeros(ncomp,size(tnew,2));
    
    for ic=1:ncomp
        %interp1q to avoid input checks...
        ynew(ic,:)=interp1q(told',yold(ic,:)',tnew')';
    end
end
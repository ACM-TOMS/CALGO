function [summod,ebigst,index,esecnd]=stats(elem)
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
    % gives stats from vector elem
    %
    % returns summod: sum of absolute values of elem
    %         ebigst: largest element in magnitude
    %         index: index of element ebigst
    %         escnd: magnitude of second largest (strictly) element

    elem=abs(elem);
    summod=sum(elem);
    [ebigst index]=max(elem);
    esecnd=max(elem(elem<ebigst));
end
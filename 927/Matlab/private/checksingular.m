function b=checksingular()
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
    [wn, wid]=lastwarn();
    %b=strcmp(wid,'MATLAB:singularMatrix')||strcmp(wid,'MATLAB:nearlySingularMatrix'); % new modify ||strcmp(wid,'MATLAB:nearlySingularMatrix')
    b=strcmp(wid,'MATLAB:singularMatrix');
   lastwarn('');
end
function f = condaux(flag,X,L,U,P,Q, D)
% ''Higham algorithm''
% CONDAUX  Auxiliary function for estimation of condition.
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
switch flag
   case 'dim'
       f = max(size(L));
   case 'real'
       f = 1;
   case 'notransp'
       f = P'*(L'\(U' \(Q'*X)));
       f = D*f;
   case 'transp'
       f=D*X;
       f =Q* (U\(L\ (P*f))); %
end

return
end

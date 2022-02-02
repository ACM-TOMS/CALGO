function [x,A]=solvesls(A,b)
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
    if ~isstruct(A)
        [L,U,P,Q]=lu(A);
        A=struct('L',L,'U',U,'P',P,'Q',Q);
    end
    
    Pb=A.P*b;
    xp=A.L\Pb;
    xpp=A.U\xp;
    x=A.Q*xpp;
end
function dfrat=ratcor(h,bhold,defimp)
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
    ninter=size(h,2);
    ncomp=size(defimp,1);
    Ni=ncomp*ninter;
    
    hb=repmat(h/2,ncomp,1); 
    defimpc=reshape(defimp,Ni,1);
    dfrat=reshape(defimpc-hb(:).*(bhold*defimpc),ncomp,ninter);
end
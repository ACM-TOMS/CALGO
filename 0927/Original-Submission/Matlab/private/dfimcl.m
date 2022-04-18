function defimp=dfimcl(chold,defexpl)
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
    nmsh=size(chold,3);
    defimp=zeros(size(chold,1),nmsh);

    ws_sing=warning('off','MATLAB:singularMatrix');
    % new change add only the next pragraph
    % maybe we can apply the same trick as in jaccal and ratcor here:
    % rewrite this stack of matrices chold as one matrix, defexpl as one
    % vector and then solve for the different defimp in one pass
    % but what happens when one of those matrices is (nearly) singular?
    for im=1:nmsh
        defimp(:,im)=chold(:,:,im)\defexpl(:,im);
    end
    
    if checksingular()
        % do call checksingular to clear the lasterror!
    end
    
    warning(ws_sing.state,'MATLAB:singularMatrix');
    
end
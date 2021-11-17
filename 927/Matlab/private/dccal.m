function [ratdc,dfexmx,incmp,inmsh,intol,derivm,dfimmx,rat1,rat2]=dccal(problem,t,y,fty,defexp,defimp,ltol,dfctol)
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
    tstrat=0.1;

    ninter=size(y,2)-1;
    ntol=size(ltol,1);
    
    [compmaxv,compmaxi]=max(abs(defexp(ltol,:)),[],2);
    [dfexmx,intol]=max(flipdim(compmaxv,1)); %flipdim here to mimic the
    intol=ntol-intol+1;                     %behaviour of the orignal code
    inmsh=compmaxi(intol);
    incmp=ltol(intol);
    
    derivm=max(abs(fty(incmp,:)));
    
    smtest=tstrat*dfexmx;
            
    if problem.debug, disp('defexp(incmp)'); end
    %if problem.debug, disp(sprintf('%g ',defexp(incmp,:))); end
    
    if problem.debug, disp('defimp(incmp)'); end
    %if problem.debug, disp(sprintf('%g ',defimp(incmp,:))); end

    texp=defexp(incmp,:);
    timp=defimp(incmp,:);
    dfimmx=max(abs(timp));
    abtexp=abs(texp);
    
    timp(abs(timp)<dfctol)=dfctol;
    ratdc=texp./timp;
    abrat=abs(ratdc);
    %if problem.debug, disp(['abrat: ' sprintf('%g ',abrat)]); end
    
    rat1=0;
    rat2=0;
    
    for im=1:ninter
        if abtexp(im) <= dfctol
            ratdc(im)=1;
        else
            rat2=max(abrat(im),rat2);
            if abtexp(im)>=smtest && abrat(im)>=rat1
                rat1=abrat(im);
            end
        end
    end
    
end
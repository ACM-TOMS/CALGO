function [t,h]=unimsh(a,b,nmsh,fixpnt)
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
    % creates a uniform mesh with 'nmsh' meshpoints, including a and b
    %
    % returns t: meshpoints
    %         h: intermesh distances


    nfxpnt=numel(fixpnt);
    if nfxpnt==0
        h=(b-a)/(nmsh-1);
        t=linspace(a,b,nmsh);
        h=repmat(h,1,nmsh-1);
    else
        if nfxpnt>0

            if nmsh<nfxpnt+2
                nmsh=nfxpnt+2;
            end

            ninter=nmsh-1;
            xx=NaN*ones(1,nmsh);
            xx(1)=a;
            xx(end)=b;

            ileft=1;
            xleft=a;

            totint=b-a;
            ndif=ninter-nfxpnt;

            for j=1:nfxpnt+1
                if j<=nfxpnt
                    xright=fixpnt(j);
                    nmin=ninter*(xright-a)/totint+1.5;
                    if nmin>ndif+j
                        nmin=ndif+j;
                    end
                    iright=floor(max(ileft+1,nmin));
                else
                    xright=b;
                    iright=nmsh;
                end

                xx(iright)=xright;
                npt=iright-ileft-1;
                dx=(xright-xleft)/(npt+1);
                for i=1:npt
                    xx(ileft+i)=xleft+i*dx;
                end
                ileft=iright;
                xleft=xright;
            end

            t=xx;
            h=xx(1,2:end)-xx(1,1:end-1);
        end
    end
end

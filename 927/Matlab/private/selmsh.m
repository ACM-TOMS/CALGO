function [t,h,nmsh,maxmsh,double,told]=selmsh(problem,ltol,tol,fixpnt,ipow,nmax,told,u,ermeas)
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

    if problem.debug, disp(sprintf('selmsh, ipow=%g',ipow)); end

    erdcid=5;
    rlndec=log(erdcid);
    frcpow=1/ipow;
    ithres=0;
    thres=1;
    phitst=0.1;
    
    nmsh=size(u,2);
    ninter=nmsh-1;
    nmold = nmsh;
    maxmsh=0;
    double=0;
   
    denom=abs(u(ltol,1:end-1));
    denom(denom<1)=1;
    denom=denom.*repmat(tol,1,ninter);
    ermeas(ltol,:)=abs(ermeas(ltol,:))./denom;
    
    [ermx,ihcomp]=max(ermeas,[],1);
    errmax=max(ermx);
    
    if problem.debug, disp(sprintf('errmax %g',errmax)); end
    
    if errmax>0 && errmax <= erdcid
        if errmax>1
            ii=1;
            decii=erdcid;
        else
            ilg=floor(-log(errmax)/rlndec);
            ii=2+ilg;
            decii=erdcid^ii;
        end
        
        errmax=decii*errmax;
        ermx=decii*ermx;
        ermeas=decii*ermeas;
    end
    
    if problem.debug, disp('ermx'); end
    if problem.debug, disp(sprintf('%g ',ermx)); end
    
    selmshdone=0;
    toolargemsh=0;
    
    while ~selmshdone
        
        nmest=nmsh;
        irefin=floor(ermx.^frcpow)+1;
        ermxthr=ermx<thres;
        irefin(ermxthr)=1;
        nmest=nmest+sum(irefin(~ermxthr))-size(ermx,2);
        
        if problem.debug, disp(sprintf('nmest %g, irefin:',nmest)); end
        %if problem.debug, disp(sprintf('%g ',irefin)); end
        
        if nmest>nmax
            toolargemsh=1;
        elseif nmest-1 >= 3*ninter
            double=1;
            [t,h,nmsh,maxmsh,told]=dblmsh(problem,told,nmax);
            return;
        else
            %let's refine the mesh

            new=1;
            t=NaN*ones(1,nmest+sum(irefin==1));
            t(1)=told(1);

            rlen=told(2)-told(1);
            slen=rlen;
            if irefin(1)>1
                dx=rlen/irefin(1);
                s=2:irefin(1);
                t(s)=t(1)+dx*(s-1);
                new=new+irefin(1)-1;
            end

            ifxcnt=1;
            if numel(fixpnt)==0
                fxnext=1.1*abs(told(end));
            else
                fxnext=fixpnt(1,ifxcnt);
            end

            jtkout=0;

            for im=2:ninter
                rlold=rlen;
                rlen=told(im+1)-told(im);

                remove=0;

                if told(im)==fxnext
                    ifxcnt=ifxcnt+1;
                    if ifxcnt>size(fixpnt,2)
                        fxnext=1.1*told(end);
                    else
                        fxnext=fixpnt(ifxcnt);
                    end
                elseif irefin(im)==1
                    slen=slen+rlen;

                    if jtkout==0
                        ind1=ihcomp(im-1);
                        phihat=ermeas(ind1,im-1)/(rlold^ipow);
                    end
                    phihat=max(phihat,ermeas(ihcomp(im),im)/(rlen^ipow));
                    val1=phihat*(slen^ipow);
                    if val1<=phitst && jtkout<4
                        jtkout=jtkout+1;
                        remove=1;
                    end
                end
                if ~remove
                    jtkout=0;

                    new=new+1;
                    t(new)=told(im);
                    if irefin(im)>1
                        dx=rlen/irefin(im);
                        s=(1:irefin(im)-1);
                        t(new+s)=told(im)+s*dx;
                        new=new+irefin(im)-1;
                    end
                    slen=rlen;

                    if new>nmax
                        toolargemsh=1;
                        break;
                    elseif new>3*ninter
                        double=1;
                        [t,h,nmsh,maxmsh,told]=dblmsh(problem,told,nmax);
                        return;
                    end
                end
            end
        end
        if toolargemsh
            toolargemsh=0;
            if 2*nmsh-1<nmax
                double=1;
                [t,h,nmsh,maxmsh,told]=dblmsh(problem,told,nmax);
                return;
            elseif thres<errmax && ithres<3
                ithres=ithres+1;
                thres=erdcid*thres;
                thres=min(errmax,thres);
                t=told;
            else
                nmsh=nmold;
                t=told;
                maxmsh=1;
                selmshdone=1;
            end
        else
            new=new+1;
            t(new)=told(end);
            nmsh=new;
            maxmsh=0;
            
            t=t(1:new);
            
            if problem.debug, disp(sprintf('selmsh: nmsh: %g new mesh:',nmsh)); end
            if problem.debug, disp(sprintf('%g ',t)); end
            
            selmshdone=1;
        end
    end
    h=t(2:end)-t(1:end-1);
end
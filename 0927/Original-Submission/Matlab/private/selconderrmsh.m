function [t,h,nmsh,maxmsh,double,told]=selconderrmsh(problem,ltol,tol,fixpnt,ipow,nmax,told,u,ermeas,omega,stabcond,moncondmsh_par)
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

    if problem.debug, disp(sprintf('selconderrmsh, ipow=%g',ipow)); end

    erdcid=5;
    rlndec=log(erdcid);
    frcpow=1/ipow;
    ithres=0;
    thres=1;
    phitst=0.1;
    
    nmsh=size(u,2);
    ninter=nmsh-1;
    
    maxmsh=0;
    double=0;
    ermeas_old=ermeas;
    for im=1:ninter
        denom=abs(u(ltol,im));
        denom(denom<1)=1;
        denom=denom.*tol;
        ermeas(ltol,im)=abs(ermeas(ltol,im))./denom;    
    end

    [ermx,ihcomp]=max(ermeas,[],1);
    errmax=max(ermx);
   
    if problem.debug, disp(sprintf('errmax %g',errmax)); end
    h=told(2:end)-told(1:end-1);
    if (  ~stabcond.all && errmax >= 1e20 ) 
%   only the conditioning      
        [t,h,nmsh,maxmsh]=selcondmsh(problem,told,fixpnt,nmax,nmsh,omega,h,problem.a,problem.b,moncondmsh_par);                       
        
    else
    
     [nptcond,cond_epsi,fatt_r1r2,fatt_r2] = moncondmsh(nmsh,omega,h,problem.a,problem.b,moncondmsh_par);
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
        irefin=fix(ermx.^frcpow)+1;
        ermxthr=ermx<thres;
        irefin(ermxthr)=1;
        nmest=nmest+sum(irefin(~ermxthr))-size(ermx,2);
        
        if (nptcond >= 4)          
            indrefin = find( cond_epsi(1:end)>=fatt_r1r2);
            nmest = nmest + sum(max(0,nptcond-irefin(indrefin)-1));
            irefin(indrefin)=max(nptcond, irefin(indrefin) );
        end
        
        
        if problem.debug, disp(sprintf('nmest %g, irefin:',nmest)); end
        if problem.debug, disp(sprintf('%g ',irefin)); end
        
        if nmest>nmax
            toolargemsh=1;
        elseif nmest-1 >= 3*ninter
            double=1;
            [t,h,nmsh,maxmsh,told]=dblmsh(problem,told,nmax);
            return;
        else
            %let's refine the mesh

            new=1;
            t=told(1);

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
                fxnext=1.1*told(end);
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
                    if val1<=phitst && jtkout<4 && cond_epsi(im) < fatt_r1r2
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
                nmsh=length(told);
                maxmsh=1;
                t=told;
                selmshdone=1;
            end
        else
            new=new+1;
            t(new)=told(end);
            nmsh=new;
            maxmsh=0;
            selmshdone=1;
        end
     end
    
    end 
    h=t(2:end)-t(1:end-1);
    
    if problem.debug, disp(sprintf('selconderrmsh: nmsh: %g new mesh:',nmsh)); end
    if problem.debug, disp(sprintf('%g ',t)); end
end
function [t,h,nmsh,maxmsh,double,told,nvold,voldmsh,phiold,problem]=selmsh_a(problem,nvold,voldmsh,phiold,ltol,tol,fixpnt,ipow,nmax,t,told,u,ermeas)
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

if problem.debug, fprintf('selmsh, ipow=%g\n',ipow); end

erdcid=5;
rlndec=log(erdcid);
frcpow=1/ipow;
ithres=0;
thres=1;
phitst=0.1;
nmsh=length(t);
ninter=nmsh-1;
maxmsh=0;
double=0;
nmold = nmsh;
problem.iprec = min(problem.iprec,1);
if (problem.iatt == -1), nvold = nmsh; end
if (problem.iatt == 0), voldmsh = told; end
told = t;
denom = min(abs(u(ltol,1:end-1)),abs(u(ltol,2:end)));
denom(denom<1)=1;
denom=denom.*repmat(tol,1,ninter);
ermeas(ltol,:)=abs(ermeas(ltol,:))./denom;

[ermx,ihcomp]=max(ermeas,[],1);
errmax=max(ermx);

philrg = 0;
esum = 0;
mshchng = 0;
selmshdone=0;
toolargemsh=0;

%default value in case selmsh_a is ended abruptly
h=t(2:end)-t(1:end-1);

while ~selmshdone
    
    nmest = nmsh;
    errim = ermx.^frcpow;
    irefin=floor(errim)+1;
    ermxthr=ermx<thres;
    irefin(ermxthr)=1;
    
    nmest=nmest+sum(irefin(~ermxthr))-size(ermx,2);
    mshchng =max(irefin) > 1;
        
    him = told(2:end)-told(1:end-1);
    phiim = errim./him;
    esum = sum(errim);
    if (problem.iatt == -1)
        phiold = phiim;
    else
        [philrg,imreg] = max(phiim);
        hordlrg = errim(imreg);
    end
    
    if (problem.iatt == 0)
        problem.hord(2) = hordlrg;
        problem.pmax(2) = philrg;
        xloc1 = told(imreg);
        xloc2 = told(imreg+1);
        ichkpt = floor(nvold/2);
        if (xloc1 < voldmsh(ichkpt)), ichkpt = 1; end
        xb = voldmsh(ichkpt);
        for j = ichkpt:nvold-1
            xa = xb;
            xb = voldmsh(j+1);
            if (xloc1 >= xa && xloc1 < xb)
                if (xloc2-xb < (xloc2-xloc1)/2.0)
                    problem.pmax(1) = phiold(j);
                    problem.hord(1) = (xb-xa)*problem.pmax(1);
                else
                    problem.pmax(1) = phiold(j+1);
                    problem.hord(1) = (voldmsh(j+2)-xb)*problem.pmax(1);
                end
            end
        end
        problem.hsml = (xloc2-xloc1)/irefin(imreg);
        problem.npr(2) = floor(esum);
    end
    
    if (problem.iatt == -1), problem.npr(1) = floor(esum); end
    if (problem.ifinal == 1 && mshchng == 0), problem.iatt = 0; end
    
    
    
    
    if (nmest > nmax)
        toolargemsh=1;
    else
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
            if (t(s) == t(s-1))
                problem.iprec = 2;
                return
            end
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
                if (t(new) == t(new-1))
                    problem.iprec = 2;
                    return
                end
               if irefin(im)>1
                    dx=rlen/irefin(im);
                    s=(1:irefin(im)-1);
                    t(new+s)=told(im)+s*dx;
                    new=new+irefin(im)-1;
                    if (t(new) == t(new-1))
                        problem.iprec = 2;
                        return
                    end
                end
                slen=rlen;
                
                if new>nmax-1
                    toolargemsh=1;
                    break;
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
        if (t(new) == t(new-1))
            problem.iprec = 2;
            return
        end
        nmsh=new;
        maxmsh=0;
        t=t(1:new);
        
        if problem.debug, fprintf('selmsh: nmsh: %g\n ',nmsh); end % new mesh:
    %    if problem.debug, fprintf('%g',t); end
        selmshdone=1;
    end
end
h=t(2:end)-t(1:end-1);
end
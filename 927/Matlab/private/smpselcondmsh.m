function [t,h,nmsh,maxmsh,ddouble,told]=smpselcondmsh(problem,t,intref,numadd,nmax,nmsh,omega,h,a,b,moncondmsh_par)
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
% inserts points into mesh t
% if intref=1, numadd points are inserted in the first interval
% if intref=nmsh-1, numadd points are inserted in the last interval
% else numadd points are added in interval [intref-1,intref,intref+1]
% only if the new number of meshpoints does not exceed nmax
%
% returns t: new meshpoints
%         h: new intermesh distances
%         nmsh: new number of meshpoints
%         maxmsh: if nmax would be exceeded by this action
if problem.debug, disp(sprintf('smpselcondmsh:  mesh init: %g points',nmsh)); end
fixpnt=problem.fixpnt;
maxmsh=0;
ddouble=0;
nmsh=size(t,2);
ninter=nmsh-1;
numadd=min(49,max(numadd,4));
nmold = nmsh;
told=t;
nmest=nmsh;
toolargemsh=0;

[nptcond,cond_epsi,fatt_r1r2,fatt_r2] = moncondmsh(nmsh,omega,h,a,b,moncondmsh_par);


irefin = ones(nmsh-1,1);
indrefin = find( cond_epsi(1:end)> fatt_r1r2);
irefin(indrefin)=nptcond;
nmest=nmest+sum(irefin(indrefin))-length(irefin(indrefin));

if intref==1
    irefin(1)=max(numadd,irefin(1));
    nmest=nmest+numadd-1;
elseif intref==nmsh-1
    irefin(nmsh-1)=max(numadd,irefin(nmsh-1));
    nmest=nmest+numadd-1;
else
    numadd=max(min(numadd,9),4);
    irefin(intref-1)=max(numadd, irefin(intref-1) );
    irefin(intref)=max(numadd,irefin(intref));
    irefin(intref+1)=max(numadd,irefin(intref+1));
    nmest=nmest+3*numadd-1;
end



if problem.debug, disp(sprintf('nmest %g, irefin:',nmest)); end
if problem.debug, disp(sprintf('%g ',irefin)); end
if problem.debug, disp(sprintf('intref %g',intref)); end
if nmest>nmax
    toolargemsh=1;
  
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
            
            
            if  jtkout<2&& cond_epsi(im) < 5e-1*fatt_r2
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
        end
        
        if new>nmax
            toolargemsh=1;
            break;
        elseif new>3*ninter
             ddouble=1;
             [t,h,nmsh,maxmsh,told]=dblmsh(problem,told,nmax);
             return;
        end
    end
end

if toolargemsh
    if 2*nmsh-1<nmax
        ddouble=1;
        [t,h,nmsh,maxmsh,told]=dblmsh(problem,told,nmax);
        return;
    else
        nmsh=nmold;
        maxmsh=1;
        t = told;
    end
else
    new=new+1;
    t(new)=told(end);
    nmsh=new;
    maxmsh=0;
    
end

h=t(2:end)-t(1:end-1);


if problem.debug, disp(sprintf('smpselcondmsh: new mesh: %g points',nmsh)); end
end
function [J,bhold,chold]=jaccal(problem,h,t,y,fty)
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
    nmsh=size(y,2);
    ninter=nmsh-1;
    N=nmsh*problem.ncomp;
    Ni=problem.ncomp*ninter;
    
    % some index vectors to be reused by makedsq and then (mm) buildspjac
    [m2j,m2ib]=meshgrid(1:Ni,1:problem.ncomp);
    m2ic=reshape(meshgrid(0:problem.ncomp:Ni-1,1:problem.ncomp*problem.ncomp),problem.ncomp,Ni);
    m2i=m2ib(:)+m2ic(:);
    m2j=m2j(:);
    %%

    yn=y(:,1:end-1);
    ynpo=y(:,2:end);
    tn=t(1,1:end-1);
    tnpo=t(1,2:end);
    fn=fty(:,1:end-1);
    fnpo=fty(:,2:end);
    dfty=problem.df(t,y);
    dfn=dfty(:,:,1:end-1);
    dfnpo=dfty(:,:,2:end);
    
    hb=repmat(h,problem.ncomp,1);
    
    ynh=(yn+ynpo)./2+hb.*(fn-fnpo)./8;
    tnh=(tnpo+tn)./2;
    
    dftm1=dfn;
    dftm2=problem.df(tnh,ynh);
    dsq=makedsq(dftm2,dftm1);
    
    I=repmat(eye(problem.ncomp,problem.ncomp),[1,1,ninter]);
    
    h=repmat(reshape(h,1,1,nmsh-1),[problem.ncomp,problem.ncomp,1]);
    
    ajac1=-I-h.*(dftm1./2+dftm2+h.*dsq./4)./3;
    
    dftm1=dfnpo;
    [dsq,bhold]=makedsq(dftm2,dftm1);
    chold=I-h.*(dftm1./2+dftm2-h.*dsq./4)./3;
    %ajac2=chold;
    
    [dga,dgb]=problem.dg(y(:,1),y(:,end));
    J=buildspjac(ajac1,chold,dga,dgb);
    
    function [idsq,m2s]=makedsq(m1,m2)
        
        m2s=sparse(m2i,m2j,m2(:),Ni,Ni);
        idsq=reshape(reshape(m1,problem.ncomp,Ni)*m2s,problem.ncomp,problem.ncomp,ninter);
    end

    function spJ=buildspjac(ajac1,ajac2,bca,bcb)
        
        [ljm,lim]=meshgrid([1:problem.ncomp Ni+1:N],1:problem.ncomp);
        im=lim(:);
        jm=ljm(:);
        vm=reshape(cat(2,bca,bcb),problem.ncomp*problem.ncomp*2,1);
        
        m2i=m2i+problem.ncomp;
        im=cat(1,im,cat(1,m2i,m2i));
        jm=cat(1,jm,cat(1,m2j,m2j+problem.ncomp));
        vm=cat(1,vm,cat(1,ajac1(:),ajac2(:)));
        
        spJ=sparse(im,jm,vm,N,N);
    end
end
function [y,iflag,problem,rhs,rnsq,fty]=fixjac(problem,iorder,ltol,tol,t,y,fty,defcor,defnew,rhs,J)
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
    lmtfrz=8;
    huge=1e30;
    rngrow=16;
    xlarge=1e6;
    tolfct=0.1;
    
    epsmch=1.11022302E-16;
    
    if problem.debug, disp('fixed Jacobian iterations'); end
    
    h=t(2:end)-t(1:end-1);
    
    rhstri=rhs;
    rhstri(problem.ncomp+1:end)=rhstri(problem.ncomp+1:end)+reshape(defnew,[],1);
    rnsq=sum(rhstri.^2);
    
    if rnsq>huge || ( iorder==8 && rnsq>xlarge )
        if problem.debug, (fprintf('Large residual, rnsq = %g',rnsq)); end
        iflag=-2;
        return
    end
    rhs=rhstri;
    
    iter=0;
    while true
        
        if problem.debug, (fprintf('iter %d rnsq %g',iter,rnsq)); end
        if rnsq<=epsmch
            iflag=0;
            return;
        end
        iter=iter+1;
        
        %dy=J\rhs;
        dy=solvesls(J,rhs);
        delu=reshape(dy,problem.ncomp,[]);
        ytrial=y+delu;
        
        rnold=rnsq;
        [rhstri,rnsq,fty,problem]=rhscal(problem,h,t,ytrial,defcor);
        
        better=false;
        if rnsq<rnold
            better=true;
            y=ytrial;
            rhs=rhstri;
        end
        
        if iter>=lmtfrz || rnsq>rnold/rngrow
            if better
                iflag=-3;
            else
                iflag=-2;
            end
            if problem.debug, (fprintf('failure of fixed jacobian, iflag = %d',iflag)); end
            return
        end
        
        errok=true;
        for it=1:size(ltol,1)
            itol=ltol(it);
            erdenom=abs(y(itol,:));
            erdenom(erdenom<1)=1;
            er=abs(delu(itol,:))./erdenom;
            % implemented as a negation of a > because er can be NaN
            % the original code accepts a NaN as lower than the tolerance
            errok=errok & ~any(er>tolfct*tol(it));
        end
        
        if errok
            if problem.debug, (fprintf('fixed jacobian convergence, iter: %d, rnsq: %g',iter,rnsq)); end
            iflag=0;
            return;
        end
    end
end
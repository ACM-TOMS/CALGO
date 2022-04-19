function [y,iflag,rhs,problem,fty,it,J,bhold,chold]=newteq(problem,h,t,y,dc,ltol,tol,newton_par)
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
    it=-1;
    iflag=0;
    lmtnwt=39;
    ncomp=size(y,1);
    nmsh=size(y,2);
    
    imerit=1;
    alfsml=1e-4;
    alfmax=1.1;
    stpfct=2.0;
    cnvfct=0.1;
    
    alfold=1;
    grfct=newton_par.grfct;   
    greps=newton_par.greps;
    itwtmx=8;
    itonew=5;
    itwtch=0;
    wmprev=realmax;
    wmbest=realmax;
    
   
    epsmch=eps/2;
    
    Ni=ncomp*nmsh;
   
    
    ws_sing=warning('off','MATLAB:singularMatrix');
    ws_nsing=warning('off','MATLAB:nearlySingularMatrix');
    
    [rhs,rnsq,fty,problem]=rhscal(problem,h,t,y,dc);
    
    while 1
        it=it+1;
        if problem.debug, disp(sprintf('newton iteration %g',it)); end
        
        if it>=lmtnwt
            %error('Too many Newton iterations, iter=%g',it)
            if problem.debug, disp('too many newton iterations'); end
            iflag=-2;
            break;
        end
        
        if rnsq<=greps*epsmch %todo: make this eps
            if problem.debug, disp(sprintf('rnsq small: %g',rnsq)); end
            [J,bhold,chold]=jaccal(problem,h,t,y,fty);
            [dy,J]=solvesls(J,rhs);
            if checksingular()
                iflag=-1;
                if problem.debug, disp('singular Jacobian'); end
            
                break;
            end
            delu = reshape(dy,ncomp,nmsh);
            %if problem.debug, disp(rhs); end
            break;
        end
        
        if wtchdg(it,rnsq)~=0
            %error('Watchdog tests fail, iter=%g',it);
            if problem.debug, disp('watchdog abort'); end
            iflag=-3;
            break;
        end
        [J,bhold,chold]=jaccal(problem,h,t,y,fty);
        
        %dy=J\rhs;
        [dy,J]=solvesls(J,rhs);
        if checksingular()
            iflag=-1;
            if problem.debug, disp('singular Jacobian'); end
            break;
        end
        delu = reshape(dy,ncomp,nmsh);
        
        fmtry=sum(dy.^2);
        
        alfa=1;
        if stpfct*alfold < 1
            alfa=stpfct*alfold;
        end
        
        if alfa < alfsml
            alfa=alfsml;
        end
        
        if problem.debug, disp(sprintf('rnsq %g',rnsq)); end
        [alfa,iflag,alfold]=getptqm(problem,@meritfunct,alfa,fmtry,alfold,rnsq);
        if iflag~=0
            if problem.debug, disp(sprintf('getptq error %d',iflag)); end
            %error('getptq');
            break;
        end
        
        y=y + alfa*delu;
      
        wmprev=rnsq;
        [rhs,rnsq,fty,problem]=rhscal(problem,h,t,y,dc);
        
        conv=true;
        for itc=1:size(tol,1)
            icmp=ltol(itc);
            erdenom=abs(y(icmp,:));
            erdenom(erdenom<1)=1;
            er=abs(delu(icmp,:))./erdenom;
            conv=conv&all(er<=cnvfct*tol(itc));
        end
        if conv
            if problem.debug, disp(sprintf('Convergence, iter = %g , rnsq= %g',it+1,rnsq)); end
            break;
        end
    end
    
    warning(ws_nsing.state,'MATLAB:nearlySingularMatrix');
    warning(ws_sing.state,'MATLAB:singularMatrix');
    
    function iflag=wtchdg(iter,wmerit)
        if problem.debug, disp(sprintf(' iter %g wmerit %g wmbest %g wmprev %g', iter,wmerit,wmbest,wmprev)); end
        if problem.debug, disp(sprintf(' itwtch %g alfold %g',itwtch,alfold)); end
        iflag=0;
        if wmerit<=wmbest
            itwtch=0;
            wmbest=wmerit;
        else
            itwtch=itwtch+1;
            if alfold>=1/2
                iflag=0;
                return
            end
            if wmerit<=wmprev && itwtch <= 2*itwtmx 
                iflag=0;
            elseif itwtch >= itwtmx 
                iflag=-1;
            elseif iter>=itonew && wmerit>grfct*wmbest
                iflag=-2;
            end
        end
    end

    function  [fmtry,rnsqtr]=meritfunct(x)
        utrial= x*delu+y;
        [rhstri,rnsqtr]=rhscal(problem,h,t,utrial,dc);
        if imerit
           %xmerit=J\rhstri;
           xmerit=solvesls(J,rhstri);
           fmtry=sum(xmerit.^2);
        else
           fmtry=rnsqtr;
        end
        if problem.debug, disp(sprintf('alpha %g merit %g',x,fmtry)); end
    end
end

 
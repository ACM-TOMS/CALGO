function [y,t,err,succes,maxmsh,stabcond,condpar,problem,iseries]=twpbvp_m(problem,nmsh,nmax,ltol,tol,liseries)
%
%   Private function implementing the twpbvp_m and twpbvpc_m solvers 
%   based on the fortran solver twpbvpc
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
pdebug = problem.debug;
iseries = zeros(1,liseries);
linear=problem.linear;
use_c=problem.conditioning;
comp_c=problem.conditioning;
if ~comp_c
    use_c = 0;
    stabcond.stiff =0;
    stabcond.all=0;
end
itcondmax=10;
itcond=0;
condpar.kappa = realmax;
condpar.gamma1 = realmax;
condpar.kappa1 = realmax;
condpar.sigma=0;
double=0;
dfold=0;
oldrt1=realmax;
small=1e-2;
bigfac=10;
quan6 = 0.1;
newton_par.grfct = 100;
newton_par.greps = 1;

moncondmsh_par.sfatt_alpha = 0.08;
moncondmsh_par.sfatt_r1r2 = 0.5;
moncondmsh_par.sfatt_r2 = 1e-5;
compcond_par.fatt_nl_sigma = 5;
compcond_par.fatt_lin_sigma = 10;

%  The routine stcons calculates integration constants
[c_alp, c_bet,c_abc]=stcons();
% etest6
if linear
    etest6 = tol.^(1/3);
    etest6(etest6<quan6)=quan6;
    etest6=1./etest6;
else
    etest6=ones(size(tol));
end
etest6=10*ones(size(tol));

%mshref
numbig=0;
nummed=0;
numpt=14;
bigfac=10;
% first 4
first4 = 1;


smooth=0;
reaft6=0;
strctr=0;
trst6=1;
maxmsh=0;
er6old=realmax;
er8old=realmax;
quan8=0.025;
efact=100;
fxfct=10;
epsmch=1.11022302E-16;

if isscalar(nmsh)
    nmold=nmsh;
    [t,h]=unimsh(problem.a,problem.b,nmsh,problem.fixpnt);
else
    t=nmsh;
    h=t(2:end)-t(1:end-1);
    nmsh=size(t,2);
    nmold=nmsh;
end
y=initu(problem,t,problem.ncomp);
done=0;
succes=0;
err=0;
%in case rhs is small enough in the first iteration
ncomp=problem.ncomp;
nfail4 = 0; % number of newton failure starting from the initial guess

indnmesh = 0;

while ~done
    if pdebug
        (fprintf('========== nmsh = %g:', nmsh));
    end
    
    indnmesh = indnmesh+1;
    iseries(indnmesh) = nmsh;
    if indnmesh > liseries
        succes = 2;
        done =1;
        return
    end
    ninter=nmsh-1;
    dc=zeros(problem.ncomp,ninter);
    if (linear)
        ludone = 0;
        [y,iflag,problem,rhs,fty,J,bhold,chold] = lineq(problem,h,t,nmsh,y,dc,ludone);
        fty = problem.f(t,y);
        if ~problem.vectorized
            problem.NFUN = problem.NFUN + length(t);
        else
            problem.NFUN = problem.NFUN + 1;
        end
    else
        [y,iflag,rhs,problem,fty,itnwt,J,bhold,chold]=newteq(problem,h,t,y,dc,ltol,tol,newton_par);
    end
    if iflag==0
        if comp_c
            condpar_old=condpar;
            [omega,condpar,stabcond] = comp_condpar(J,problem.ncomp,nmsh,h,problem.a,problem.b,condpar_old,linear,compcond_par);
            if pdebug
                (fprintf('sigma = %g gamma1 = %g kappa1 = %g kappa = %g',condpar.sigma,condpar.gamma1,condpar.kappa1,condpar.kappa));
                (fprintf('stab_kappa = %g stab_gamma1 = %g stab_kappa1 = %g stiff_cond = %g ill_cond = %g',stabcond.kappa,stabcond.gamma1,stabcond.kappa1,stabcond.stiff,stabcond.ill));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%% conv4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if pdebug, disp('conv4'); end
        
        if first4
            dfold =0;
            oldrt1=realmax;
            savedu=0;
            reposs=0;
            first4=0;
        end
        succes=0;
        
        [defexp,problem]=dfexcl_m(problem,h,t,y,fty);
        def=defexp;
        
        if smooth || reaft6
            onto6=1;
            
            if pdebug, (fprintf('proceed to order 6 smooth %d, reaft6 %d',smooth,reaft6)); end
        else
            defimp=dfimcl(chold,defexp);
         
            dfctol=1e5*1.110E-16; %to enable matching of results with original code
            [ratdc,dfexmx,incmp,inmsh,intol,derivm,dfimmx,rat1,rat2]=dccal(problem,t,y,fty,defexp,defimp,ltol,dfctol);
            %if pdebug, (fprintf('dfexmx %g incmp %g inmsh %g intol %g derivm %g /n',dfexmx,incmp,inmsh,intol,derivm)); end
            %if pdebug, (fprintf('rat1 %g rat2 %g dfimmx %g /n',rat1,rat2,dfimmx)); end
            
            tolval=tol(intol);
            [onto6,smooth,callrt,strctr,oscchk,double,reposs]=decid4(problem,linear,rat1,rat2,dfexmx,dfimmx,derivm,dfold,tolval,oldrt1); % todo: fix linear
            %if pdebug, (fprintf('rat1 %g oldrt1 %g onto6 %g strctr %g /n',rat1,oldrt1,onto6,strctr)); end
            oldrt1=rat1;
            dfold=dfexmx;
            
            if callrt
                if pdebug, sprintf('callrt /n'); end
                if pdebug, sprintf('called ratcor /n'); end
                def=ratcor(h,bhold,defimp);
                %if pdebug, disp(['ratcor ',sprintf('%g ',def)]); end
            elseif linear
                if oscchk
                    if pdebug, disp('oscchk'); end
                    [double, onto6, trst6,smooth,inmsh] = osc(problem,nmsh,incmp,defexp,dfexmx,ratdc,trst6,smooth,double,onto6,inmsh); %call osc
                elseif reposs
                    if pdebug, disp('reposs'); end
                    if savedu
                        adjerr = 0;
                        [rerr,remax,itlmx]=rerrvl(y,ysavedu,ltol,tol,adjerr); %richardson interpolation
                       
                         if (remax < 1.5*tol(itlmx))
                             err =remax;
                             succes = 1;
                             done=1;
                             disp('ri 0 ')
                             if pdebug, disp('done '); end
                         end
                    end
                    ysavedu=y;
                    double= 1;
                    savedu = 1;
                    onto6 = 0;
                end
            end
        end
        if ~onto6 && ~done
            if (double)
                nmold=nmsh;
                
                
                if pdebug, disp('double ~onto 6'); end
                [t,h,nmsh,maxmsh,told]=dblmsh(problem,t,nmax);
                
            else
                
                drat=dfexmx/(max(1,abs(y(incmp,inmsh)))*tol(intol));
                if pdebug, (fprintf('drat=%g, u(incmp,inmsh)=%g, tol(intol)=%g',drat,y(incmp,inmsh),tol(intol))); end
                numadd=fix(drat^(1.0/6.0));
                if pdebug, (fprintf('numadd %g',numadd)); end
                nmold=nmsh;
                if use_c && stabcond.stiff
                    [t,h,nmsh,maxmsh,double,told]=smpselcondmsh(problem,t,inmsh,numadd,nmax,nmsh,omega,h,problem.a,problem.b,moncondmsh_par);
                else
                    [t,h,nmsh,maxmsh]=smpmsh(problem,t,inmsh,numadd,nmax);
                end
                
            end
            if ~maxmsh
                y=initu(problem,t,problem.ncomp);
            else
                if pdebug, disp('too many meshpoints'); end
                onto6=0;
                onto8=0;
                done=1;
            end
        end
    elseif iflag==-1
        %%%%%%%%%%%%%%%%%%% fail4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if pdebug, disp('fail4'); end
        onto6=0;
        
        if pdebug, (fprintf('singulier, double mesh %g',nmsh*2-1)); end
        nmold=nmsh;
        [t,h,nmsh,maxmsh,told]=dblmsh(problem,t,nmax);
        if maxmsh
            if pdebug, disp('too many meshpoints'); end
            onto6=0;
            onto8=0;
            done=1;
        end
        ysave=y;
        y=initu(problem,t,problem.ncomp);
        nfail4 =0;
    else % iflag~=-1
        if comp_c
            condpar_old=condpar;
            [omega,condpar,stabcond] = comp_condpar(J,problem.ncomp,nmsh,h,problem.a,problem.b,condpar_old,linear,compcond_par);
            if pdebug
                (fprintf('sigma = %g gamma1 = %g kappa1 = %g kappa = %g',condpar.sigma,condpar.gamma1,condpar.kappa1,condpar.kappa));
                (fprintf('stab_kappa = %g stab_gamma1 = %g stab_kappa1 = %g stiff_cond = %g ill_cond = %g',stabcond.kappa,stabcond.gamma1,stabcond.kappa1,stabcond.stiff,stabcond.ill));
            end
        end
        %%%%%%%%%%%%%%%%%%% fail4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if pdebug, disp('fail4'); end
        onto6=0;
        nfail4=nfail4+1;
        told=t;
        maxmsh=mshref(4);
        if ~maxmsh
            ysave=y;
            if linear || itnwt==0 || nfail4 >=3
                if pdebug, disp('initu'); end
                y=initu(problem,t,problem.ncomp);
                nfail4=0;
            else
                if pdebug, disp('interp'); end
                y=interpu(told,y,t);
            end
        else
            if pdebug, disp('too many meshpoints'); end
            onto6=0;
            onto8=0;
            done=1;
        end
    end
    
    if succes
        done=1;
        onto6=0;
    end
    
    
    if onto6
        
        savedu=0;
        if pdebug, disp('start 6th order'); end
        %%***** logic for trying to calculate 6th order solution *****
        
        yold=y;
        if (linear)
            ludone = 1;
            [y,iflag,problem] = lineq(problem,h,t,nmsh,y,def,ludone,J);
        else
            [y,iflag,problem,rhs,rnsq,fty]=fixjac(problem,6,ltol,tol,t,y,fty,def,def,rhs,J);
            
            if iflag==-3 && rnsq<fxfct*epsmch
                [y,iflag,rhs,problem,fty]=newteq(problem,h,t,y,def,ltol,tol,newton_par);
            end
        end
        
        if iflag==0
            %%%%%%%%%%%%%%%%%%%%%% conv6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if pdebug, disp('conv6'); end
            
            succes=0;
            reaft6=0;
            
            
            [errok,errsum6]=errest(ltol,tol,etest6,y,yold);
            if pdebug, (fprintf('errsum6 %g (errok %g)',errsum6,errok)); end
            if errok && trst6 % todo: should depend on oscchk
                if pdebug, disp('succes!'); end
                err=errsum6;
                succes=1;
            else
                onto8=1;
            end
        else
            %%%%%%%%%%%%%%%%%%%%%% fail6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if pdebug, disp('fail6'); end
            succes=0;
            onto8=0;
            y=yold;
            
            if pdebug, (fprintf('fail6: double %d',double)); end
            if pdebug, (fprintf('fail6: reaft6 %d',reaft6)); end
            
            if ~reaft6 || ~double
                reaft6=1;
                ysave=y;
                told=t;
                
                maxmsh=mshref(6);
                if ~maxmsh
                    y=interpu(told,y,t);
                end
            else
              
                adjrer=1;
                [rerr,remax,itlmx]=rerrvl(y,ysave,ltol,tol,adjrer);
                if pdebug, (fprintf('remax %g 8*tol %g',remax,8*tol(itlmx))); end
             
                    if use_c && stabcond.stiff
                        
                        [t,h,nmsh,maxmsh,double,told]=selconderrmsh(problem,ltol,tol,problem.fixpnt,4,nmax,told,ysave,rerr,omega(1:2:end),stabcond,moncondmsh_par);
                    else
                        [t,h,nmsh,maxmsh,double,told]=selmsh(problem,ltol,tol,problem.fixpnt,4,nmax,told,ysave,rerr);
                    end
                    if ~maxmsh
                        if ~double
                            y=interpu(told,ysave,t);
                        else
                            reaft6=1;
                            ysave=y;
                            told=t;
                            
                            maxmsh=mshref(6);
                            if ~maxmsh
                                y=interpu(told,y,t);
                            end
                        end
                    end
              
            end
        end
        
        if succes
            done=1;
            onto8=0;
        elseif maxmsh
            if pdebug, disp('too many meshpoints'); end
            done=1;
            onto8=0;
        end
        
        if onto8
            
            if pdebug, disp('start 8th order'); end
            %%***** logic for trying to calculate 8th order solution *****
            if linear, fty = problem.f(t,y);
                if ~problem.vectorized
                    problem.NFUN = problem.NFUN + length(t);
                else
                    problem.NFUN = problem.NFUN + 1;
                end
            end
            [def8,problem]=df8cal_m(problem,h,t,y,fty,c_alp,c_bet,c_abc);
            def6 = def;
            if linear
                def=def8;
            else
                def =def + def8;
            end
            told=t;
            yold=y;
            if (linear)
                ludone = 1;
                [y,iflag,problem] = lineq(problem,h,t,nmsh,y,def,ludone,J);
            else
                [y,iflag,problem,rhs,rnsq]=fixjac(problem,8,ltol,tol,t,y,fty,def,def8,rhs,J);
                
                if iflag==-3 && rnsq<fxfct*epsmch
                    [y,iflag,rhs,problem]=newteq(problem,h,t,y,def,ltol,tol,newton_par);
                end
            end
            if iflag==0
                %%%%%%%%%%%%%%%%%% conv8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if pdebug, disp('conv8'); end
                if linear
                    etest8=tol.^(1/4);
                    etest8(etest8<quan8)=quan8;
                    etest8=1./etest8;
                else
                    etest8=ones(size(tol));
                end
                
                
                
                if (linear && strctr && nmsh < nmold), etest8=ones(size(tol)); end
                etest8=10*ones(size(tol));
                
                [errok,errsum8]=errest(ltol,tol,etest8,y,yold);
                if pdebug, (fprintf('errsum8 %g (errok %g)',errsum8,errok)); end
                if errok
                    succes=1;
                else
                    if pdebug, (fprintf('err6 %g, err8 %g, er6old %g, er8old %g',errsum6,errsum8,er6old,er8old)); end
                    
                    if (nmsh<nmold && errsum6>efact*er6old && errsum8>efact*er8old) || (nmsh < 3*nmold && errsum8 > er8old)
                        nmold=nmsh;
                        [t,h,nmsh,maxmsh,told]=dblmsh(problem,t,nmax);
                        if ~maxmsh
                            er6old=errsum6;
                            er8old=errsum8;
                            
                            if linear
                                y=initu(problem,t,problem.ncomp);
                            else
                                if pdebug, disp('interp'); end
                                y=interpu(told,yold,t);
                            end
                        end
                    else
                        told=t;
                        
                        nmold=nmsh;
                        er6old=errsum6;
                        er8old=errsum8;
                        
                        if errsum8<=errsum6
                            if use_c && stabcond.stiff
                                [t,h,nmsh,maxmsh,double,told]=selconderrmsh(problem,ltol,tol,problem.fixpnt,6,nmax,t,y,def8,omega,stabcond,moncondmsh_par);
                            else
                                [t,h,nmsh,maxmsh,double,told]=selmsh(problem,ltol,tol,problem.fixpnt,6,nmax,t,y,def8);
                            end
                            yold=y;
                        else
                            if linear, etest8=ones(size(tol)); end
                            if use_c && stabcond.stiff
                                [t,h,nmsh,maxmsh,double,told]=selconderrmsh(problem,ltol,tol,problem.fixpnt,4,nmax,t,yold,def6,omega,stabcond,moncondmsh_par);
                            else
                                [t,h,nmsh,maxmsh,double,told]=selmsh(problem,ltol,tol,problem.fixpnt,4,nmax,t,yold,def6);
                            end
                        end
                        if ~maxmsh
                            if linear
                                y=initu(problem,t,problem.ncomp);
                            else
                                if pdebug, disp('interp'); end
                                y=interpu(told,yold,t);
                            end
                        end
                    end
                end
            else
                %%%%%%%%%%%%%%%%%% fail8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if pdebug, disp('fail8'); end
                
                nmold=nmsh;
                if use_c && stabcond.stiff
                    [t,h,nmsh,maxmsh,double,told]=selconderrmsh(problem,ltol,tol,problem.fixpnt,4,nmax,t,yold,def6,omega,stabcond,moncondmsh_par);
                else
                    [t,h,nmsh,maxmsh,double,told]=selmsh(problem,ltol,tol,problem.fixpnt,4,nmax,t,yold,def6);
                end
                if ~maxmsh
                    if pdebug, disp('interp'); end
                    y=interpu(told,yold,t);
                end
            end
            if succes
                done=1;
                err=errsum8;
            elseif maxmsh
                if pdebug, disp('too many meshpoints'); end
                done=1;
            end
        end
    end
end

if succes
    if pdebug, (fprintf('meshpoints: %g',nmsh)); end
end


    function maxmsh=mshref(iorder)
        
        if pdebug, (fprintf('mshref: nummed %g numbig %g',nummed,numbig)); end
        
        nodouble = (iorder == 4) && (use_c && stabcond.stiff && ~stabcond.all);
        forcedouble = 0;
        if (use_c) && (itcond == itcondmax)
            itcond = 0;
            forcedouble = 1;
        end
        
        
        rhscomp=reshape(rhs(problem.ncomp+1:end),problem.ncomp,nmsh-1);
        
        largefound=0;
        
        nup=problem.ncomp;
        if iorder>4
            nup=size(tol,1);
        end
        
        for ic=1:nup 
          
            [sumrhs,rbigst,intref,rsecnd]=stats(rhscomp(ic,:));
            tstval=bigfac*(sumrhs-rbigst)/ninter;
            
            if pdebug, (fprintf('ic %g, tst %g, bigst %g, second %g',ic,tstval,rbigst,rsecnd)); end
            
            if rbigst>=small && rbigst>=tstval 
                if rbigst<2*rsecnd
                    nummed=nummed+1;
                end
                numbig=numbig+1;
                double=0;
                if rbigst<=bigfac*rsecnd || numbig>8
                    numbig=0;
                    nummed=nummed+1;
                    if nummed>=4 && iorder==4
                        nummed=0;
                        double=1;
                    elseif nummed>=8 && iorder>4
                        nummed=0;
                        double=1;
                    end
                end
                
                if pdebug, (fprintf('numbig %g, nummed %g',numbig,nummed)); end
                
                if double
                    if pdebug, (fprintf('smpmsh: nmsh %g, intref %g, numadd %g',nmsh,intref,numpt)); end
                    if (use_c && nodouble && ~ forcedouble)
                        double = 0;
                        [t,h,nmsh,maxmsh,double,told]=selcondmsh(problem,t,nmax,nmsh,omega,h,problem.a,problem.b,moncondmsh_par);
                        itcond=itcond+1;
                    else
                        [t,h,nmsh,maxmsh,told]=dblmsh(problem,t,nmax);
                        itcond=0;
                        % [t,h,nmsh,maxmsh]=smpmsh(problem,t,intref,numpt,nmax);
                    end
                    return
                else
                    if (use_c && nodouble && ~ forcedouble)
                        [t,h,nmsh,maxmsh,double,told]=smpselcondmsh(problem,t,intref,numpt,nmax,nmsh,omega,h,problem.a,problem.b,moncondmsh_par);
                        itcond=itcond+1;
                    elseif (forcedouble && use_c)
                        double = 1;
                        itcond=0;
                        [t,h,nmsh,maxmsh,told]=dblmsh(problem,t,nmax);
                    else
                        itcond=0;
                        numadd = numpt;
                        [t,h,nmsh,maxmsh]=smpmsh(problem,t,intref,numpt,nmax);
                    end
                    return
                end
                % we decided to double, instead of just adding
            end
        end
        if ~largefound
            if pdebug, disp('no significantly large value'); end
            numbig=0;
            nummed=0;
            double=1; %flag to next iteration that the mesh has been doubled
            
            if (use_c && nodouble &&  ~forcedouble )
                double = 0;
                [t,h,nmsh,maxmsh,double,told]=selcondmsh(problem,t,nmax,nmsh,omega,h,problem.a,problem.b,moncondmsh_par);
                itcond = itcond + 1;
            else
                if pdebug, disp('double mesh'); end
                [t,h,nmsh,maxmsh,told]=dblmsh(problem,t,nmax);
                itcond=0;
            end
        end
    end

end
function [y,t,err,succes,maxmsh,stabcond,condpar,problem,iseries]= twpbvp_l(problem,nmsh,nmax,ltol,tol,liseries)
%
%   Private function implementing the twpbvp_l and twpbvpc_l solvers 
%   based on the fortran solver twpbvplc
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
ntol = size(ltol,1);
ddouble = 0;
itcondmax = 5;           % different in twpbvpc
itcond=0;

condpar.kappa = realmax;
condpar.gamma1 = realmax;
condpar.kappa1 = realmax;
condpar.sigma= 0;
numbig=0;
nummed=0;
numpt=14;
bigfac=10;
nodouble=0;
small=1e-2;
bigfac=10;

onto6 = 0;

reaft6=0;
trst6=1;

er8old=realmax;

efact=1;
fxfct=10;
epsmch=eps/2;
maxmsh=0;
newton_par.grfct = 100;
newton_par.greps = 10;

moncondmsh_par.sfatt_alpha = 1e-5;
moncondmsh_par.sfatt_r1r2 = 0.65;
moncondmsh_par.sfatt_r2 = 1e-5;
compcond_par.fatt_nl_sigma = 10;
compcond_par.fatt_lin_sigma = 10;

[A_num1,C_num1] = stcon1();
[A_num2,B_num,C_num2] = stcon2();


etest6=10*ones(size(tol));


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

  
indnm = 0;   
nfail4 = 0; % number of newton failure starting from the initial guess

while ~done 
    if pdebug
        (fprintf('================ nmsh = %g', nmsh));
    end
    indnm = indnm+1;
    iseries(indnm) = nmsh;
    if indnm > liseries
       succes = 2;
       break
    end 
    ninter=nmsh-1;
    dc=zeros(problem.ncomp,ninter);
    if (linear)
        ludone = 0;
        [y,iflag,problem,rhs,~,J,~,chold] = lineq(problem,h,t,nmsh,y,dc,ludone);
        fty = problem.f(t,y);
        if ~problem.vectorized 
            problem.NFUN = problem.NFUN + length(t);
        else
            problem.NFUN = problem.NFUN + 1;
        end
    else
        [y,iflag,rhs,problem,fty,itnwt,J,~,chold,]=newteq(problem,h,t,y,dc,ltol,tol,newton_par);
    end

    if iflag==0
        if comp_c
            condpar_old=condpar;
            [omega,condpar,stabcond] = comp_condpar(J,problem.ncomp,nmsh,h,problem.a,problem.b,condpar_old,linear,compcond_par); % new modify add linear in input
            if pdebug
                (fprintf('sigma = %g gamma1 = %g kappa1 = %g kappa = %g',condpar.sigma,condpar.gamma1,condpar.kappa1,condpar.kappa));
                (fprintf('stab_kappa = %g stab_gamma1 = %g stab_kappa1 = %g stab_cond = %g ill_cond%g',stabcond.kappa,stabcond.gamma1,stabcond.kappa1,stabcond.ill));
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%% dfexcl substitue of conv4 %%%%%%%%%%%%%%%%%%%
        [def8,def6,Jflag,problem] = dfexcl_l(problem,h,nmsh,t,y,fty,ntol,ltol,tol,linear,A_num1,C_num1);
        if (reaft6)
            onto6= 1;
        end
    elseif iflag==-1          
        %%%%%%%%%%%%%%%%%%% fail4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if pdebug, disp('fail4'); end
        Jflag=2;
        onto6=0;
        if pdebug, (fprintf('singular matrix, double mesh %g',nmsh*2-1)); end
        nmold=nmsh;
        [t,h,nmsh,maxmsh,told]=dblmsh(problem,t,nmax);
        if maxmsh
            if pdebug, disp('too many meshpoints'); end
            onto6=0;
            onto8=0;
            done=1;
        end
        y=initu(problem,t,problem.ncomp);
        nfail4=0;
    else
        if comp_c
            condpar_old=condpar;
            [omega,condpar,stabcond] = comp_condpar(J,problem.ncomp,nmsh,h,problem.a,problem.b,condpar_old,linear,compcond_par);
            if pdebug
                (fprintf('sigma = %g gamma1 = %g kappa1 = %g kappa = %g',condpar.sigma,condpar.gamma1,condpar.kappa1,condpar.kappa));
                (fprintf('stab_kappa = %g stab_gamma1 = %g stab_kappa1 = %g stab_cond = %g ill_cond%g',stabcond.kappa,stabcond.gamma1,stabcond.kappa1,stabcond.ill));
            end
        end
        %%%%%%%%%%%%%%%%%%% fail4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if pdebug, disp('fail4'); end
        reaft6=0;
        onto6=0;
        Jflag=2;
        told=t;
        nfail4=nfail4+1;
        maxmsh=mshref(4);
       if ~maxmsh
            if linear || itnwt==0 ||nfail4 >= 3
                if pdebug, disp('initu'); end
                y=initu(problem,t,problem.ncomp);
                nfail4=0;
            else
                if pdebug, disp('interp'); end
                yold=y;
                y=interpu(told,yold,t);
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
%if onto6 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        if (Jflag == 1)
      
            if (pdebug), ( fprintf('JFLAG = %g', Jflag)); end
       
                [t,h,nmsh,maxmsh,told]=dblmsh(problem,t,nmax);
                ddouble = 1;
                if maxmsh 
                if pdebug, disp('too many meshpoints'); end
                done=1;
                end
                itcond = 0;
        
            yold = y;
          
            if(~maxmsh)
                y=interpu(told,y,t); 
            end                     
            onto6 = 0;
        elseif (Jflag == 0)
            if use_c
                nodouble =  (stabcond.stiff  && ~ stabcond.all) && (use_c);
            else
                nodouble=0;
            end
            forcedouble = 0;
            if (use_c) && (itcond == itcondmax)
                    itcond = 0;
                    forcedouble = 1;    
            end
            if pdebug, ( fprintf('JFLAG = %g', Jflag)); end
            %    find where biggest deferred correction is
            smaldef=1e40;
            bigdef=0;
            icmph = 1;
        
            ixx=1;
            for iv = 1: nmsh-1;
                for iu = 1 : ntol;
                    ipoint = ltol(iu);
                    holdef= abs(def8(ipoint,iv));
                    if(smaldef > holdef),  smaldef = holdef; end
                    if(holdef > bigdef)
                        bigdef = holdef;
                     
                        icmph = iu;
                        ixx = iv;
                        intol = iu;
                    end
                end       
            end          
            [dgtm,problem] = expl(problem,h,t,y,fty);
            ix=ixx;
            siz =abs(dgtm(icmph,ixx));
            rat=siz/bigdef;
            told=t;
            if(rat > 50 && bigdef > sqrt(tol(icmph)) &&  siz > 1)
                if use_c && stabcond.stiff
                   
                        drat = bigdef/(max(1, abs(y(icmph,ix)))*tol(intol));
                        if (pdebug),(fprintf('drat %g u %g tol%g',drat,y(icmph,ix),tol(intol))); end
                        numadd = 15;
                        [t,h,nmsh,maxmsh,ddouble,told]=smpselcondmsh(problem,t,ix,numadd,nmax,nmsh,omega,h,problem.a,problem.b,moncondmsh_par);
                       itcond=itcond+1;
                else
                        drat = bigdef/(max(1, abs(y(icmph,ix)))*tol(intol));
                        if (pdebug),(fprintf('drat %g u %g tol%g',drat,y(icmph,ix),tol(intol))); end
                        numadd = 15;
                        [t,h,nmsh,maxmsh]=smpmsh(problem,t,ix,numadd,nmax);
                end
                
             
                yold = y;
                if ~maxmsh
                  y=interpu(told,yold,t);
                end
                onto6 =0;
            else
                onto6 = 1;
                if (linear && ddouble)
                    reposs=1;
                end
            end
        end         %%%% end for if Jflage
%end
    if onto6

        if pdebug, disp('start 6th order'); end
        %%%%%%%%%%%%%%%%%%%%%%%logic for 6th order %%%%%%%%%%%%%%%%%
        %  if (iprint == 1), disp('start 6th order'), end
        def = def6;
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

        if (iflag==0)
            itcond=0;
            %%%%%%%%%%%%%%%%%%%%%% conv6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if pdebug, disp('conv6'); end
            succes=0;
            reaft6=0;
            [errok,errmax6]=errest(ltol,tol,etest6,y,yold);
            if pdebug, (fprintf('errmax6 %g (errok %g)',errmax6,errok)); end
            if errok && trst6 
                if pdebug, disp('succes!'); end
                err=errmax6;
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
            if pdebug, (fprintf('fail6: double %d',ddouble)); end
            if pdebug, (fprintf('fail6: reaft6 %d',reaft6)); end
            if ~reaft6 || ~ddouble
                reaft6=1;
                ysave=y;
                told=t;
                tsave = t;
              
                maxmsh=mshref(6);
                if ~maxmsh
                    y=interpu(told,y,t);
                end
            else
                adjrer=1;
                [rerr,remax,itlmx]=rerrvl(y,ysave,ltol,tol,adjrer);
                if pdebug, (fprintf('remax %g 8*tol %g',remax,8*tol(itlmx))); end
                if remax < 8*tol(itlmx)
                    err=remax;
                    succes=1;
                    done=1;
                else
                    t=told;
                    t=tsave;
                    nmold=size(t,2);
                    nmsh=nmold;
                    if use_c && stabcond.stiff
                        [t,h,nmsh,maxmsh,ddouble,told]=selconderrmsh(problem,ltol,tol,problem.fixpnt,4,nmax,t,ysave,rerr,omega(1:2:end),stabcond,moncondmsh_par);
                    else
                        [t,h,nmsh,maxmsh,ddouble,told]=selmsh(problem,ltol,tol,problem.fixpnt,4,nmax,t,ysave,rerr);
                    end
                    if ~maxmsh
                        if ~ddouble
                            y=interpu(told,ysave,t);
                        else
                            reaft6=1;
                            ysave=y;
                            tsave = t;
                            maxmsh=mshref(6);
                            if ~maxmsh
                                y=interpu(told,y,t);
                            end
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if linear, fty = problem.f(t,y); 
            if ~problem.vectorized 
                problem.NFUN = problem.NFUN + length(t);
            else
                problem.NFUN = problem.NFUN + 1;
            end
            end
            [def8,Jflag,problem] = df8cal_l(problem,h,nmsh,t,y,fty,ntol,ltol,tol,linear,A_num2,B_num,C_num2);
            def6 = def;   
            if (Jflag==1)
                nodouble = stabcond.stiff  && ~ stabcond.all  && use_c;
                forcedouble = 0;
                told=t;
                if (use_c)
                    if ( itcond == itcondmax)
                        itcond = 0;
                        forcedouble = 1;
                    else
                        itcond = itcond + 1;
                        forcedouble = 0;
                    end
                end
                if (nodouble && ~ forcedouble)
                    ddouble = 0;
                    [t,h,nmsh,maxmsh,ddouble,told]=selcondmsh(problem,t,nmax,nmsh,omega,h,problem.a,problem.b,moncondmsh_par);
                    
                else
                    [t,h,nmsh,maxmsh,told]=dblmsh(problem,t,nmax);
                    itcond = 0;
                end
                
                yold = y;
                if(~maxmsh)                
                    y=interpu(told,yold,t);
                end
                % restart order 4
                onto8 =0;
            end
        end
        if onto8
            if linear
                def=def8;
            else
                def =def + def8;
            end
          
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
                itcond = 0;
                errok = 0;
                %%%%%%%%%%%%%%%%%% conv8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if pdebug, disp('conv8'); end

                etest8=10*ones(size(tol));
                [errok,errmax8]=errest(ltol,tol,etest8,y,yold);
                if pdebug, (fprintf('errmax8 %g (errok %g)',errmax8,errok)); end
                if errok 
                    succes=1;
                else
                    if( nmsh < nmold &&  errmax8 > efact*er8old )
                        nmold=nmsh;
                        [t,h,nmsh,maxmsh,told]=dblmsh(problem,t,nmax);
                        
                        if ~maxmsh
                            er6old=errmax6;
                            er8old=errmax8;
                            if linear
                                y=initu(problem,t,problem.ncomp);
                            else
                                if pdebug, disp('interp'); end
                                y=interpu(told,yold,t);
                            end
                        end
                    else
                      
                        nmold=nmsh;
                       
                        er8old=errmax8;

                      
                            if use_c && stabcond.stiff
                                [t,h,nmsh,maxmsh,ddouble,told]=selconderrmsh(problem,ltol,tol,problem.fixpnt,6,nmax,t,y,def8,omega,stabcond,moncondmsh_par);
                            else
                                [t,h,nmsh,maxmsh,ddouble,told]=selmsh(problem,ltol,tol,problem.fixpnt,6,nmax,t,y,def8);
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
                % told=t;
                %yold = y;
                nmold=nmsh;
                if use_c && stabcond.stiff
                    [t,h,nmsh,maxmsh,ddouble,told]=selconderrmsh(problem,ltol,tol,problem.fixpnt,4,nmax,t,yold,def6,omega,stabcond,moncondmsh_par);
                else
                    [t,h,nmsh,maxmsh,ddouble,told]=selmsh(problem,ltol,tol,problem.fixpnt,4,nmax,t,yold,def6);
                end
                if ~maxmsh
                    if pdebug, disp('interp'); end
                    y=interpu(told,yold,t);
                end
            end
            if succes
                done=1;
                err=errmax8;
            elseif maxmsh
                if pdebug, disp('too many meshpoints'); end
                done=1;
            end
        end
    end
end

    function maxmsh=mshref(iorder)

        if pdebug, (fprintf('mshref: nummed %g numbig %g',nummed,numbig)); end
  
       nodouble = (iorder == 4) &&  (stabcond.stiff && ~stabcond.all && use_c);
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
                ddouble=0;
                if rbigst<=bigfac*rsecnd || numbig>8
                    numbig=0;
                    nummed=nummed+1;
                    if nummed>=4 && iorder==4
                        nummed=0;
                        ddouble=1;
                    elseif nummed>=8 && iorder>4
                        nummed=0;
                        ddouble=1;
                    end
                end
                largefound=1;
                if pdebug, (fprintf('numbig %g, nummed %g',numbig,nummed)); end
                if ~use_c && ~comp_c
                    nodouble=0;
                end
                if ddouble
                    if pdebug, (fprintf('meshref smpselcondmsh: nmsh %g, intref %g, numadd %g',nmsh,intref,numpt)); end
                    if (use_c && nodouble && ~ forcedouble)
                        ddouble = 0;
                        [t,h,nmsh,maxmsh,ddouble,told]=selcondmsh(problem,t,nmax,nmsh,omega,h,problem.a,problem.b,moncondmsh_par);
                        itcond=itcond+1;
                        
                    else
                        [t,h,nmsh,maxmsh,told]=dblmsh(problem,t,nmax);
                        itcond=0;
                        ddouble=1;
                    end
                    return
                else
                    if (use_c && nodouble && ~ forcedouble)
                        ddouble=0;
                        [t,h,nmsh,maxmsh,ddouble,told]=smpselcondmsh(problem,t,intref,numpt,nmax,nmsh,omega,h,problem.a,problem.b,moncondmsh_par);
                        itcond=itcond+1;
                    elseif (forcedouble && use_c)
                        ddouble = 1;
                        itcond=0;
                        [t,h,nmsh,maxmsh,told]=dblmsh(problem,t,nmax);
                    else
                        
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
          

            if (use_c && nodouble &&  ~forcedouble )
                ddouble = 0;
                [t,h,nmsh,maxmsh,ddouble,told]=selcondmsh(problem,t,nmax,nmsh,omega,h,problem.a,problem.b,moncondmsh_par);
                
                itcond = itcond + 1;
            else
                ddouble=1; %flag to next iteration that the mesh has been doubled
                if pdebug, disp('ddouble mesh'); end
                [t,h,nmsh,maxmsh,told]=dblmsh(problem,t,nmax);
                itcond=0;
            end
        end
    end
 
end


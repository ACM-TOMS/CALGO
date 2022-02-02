function [y,t,err,maxmsh,voldmsh,stabcond,condpar,nmaxmesh,problem]= bvpsol_a(problem,t,y,nmax,lambda,ltol,tol,voldmsh,nmaxmesh)
%
%   Private function implementing the deferred correction scheme for acdc and acdcc 
%   based on the fortran solver acdc
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

linear=problem.linear;
use_c=problem.conditioning;
comp_c=problem.conditioning;
inumb=0;
if ~comp_c
    use_c = 0;
end
stabcond.stiff =0;
stabcond.all=0;
stabcond.ill = 0;

ntol = size(ltol,1);

strctr=0;
condpar.kappa = realmax;
condpar.gamma1 = realmax;
condpar.kappa1 = realmax;
condpar.sigma= 0;
nvold=0;

phiold = 0;

quan6 = 0.1;
quan8=0.025;
fxfct=10;
epsmch=eps/2;
maxmsh=0;
newton_par.grfct = 100;
newton_par.greps = 10;

if (linear)
    moncondmsh_par.sfatt_alpha = 1e-5;
    moncondmsh_par.sfatt_r2 = 1d-4;
else
    moncondmsh_par.sfatt_alpha = 1e-5;
    moncondmsh_par.sfatt_r2 = 1d-5;
end
moncondmsh_par.sfatt_r1r2 = 0.5;

compcond_par.fatt_nl_sigma = 100;
compcond_par.fatt_lin_sigma = 100;

[A_num1,C_num1] = stcon1();
[A_num2,B_num,C_num2] = stcon2();

if linear
    etest6 = tol.^(1/3);
    etest6(etest6<quan6)=quan6;
    etest6=1./etest6;
else
    etest6=ones(size(tol));
end
etest6=10*ones(size(tol));
if linear
    etest8=tol.^(1/4);
    etest8(etest8<quan8)=quan8;
    etest8=1./etest8;
else
    etest8=ones(size(tol));
end

h=t(2:end)-t(1:end-1);
nmsh=size(t,2);
 
told=t;
done=0;
err=0;
ncomp=problem.ncomp;
%in case rhs is small enough in the first iteration
J=sparse(ncomp*nmsh,ncomp*nmsh);
bhold=zeros(ncomp,ncomp,nmsh-1);
chold=zeros(ncomp,ncomp,nmsh-1);
indnm = 0;


while ~done
    if problem.debug
        disp(sprintf('================ nmsh = %g\n', nmsh));
    end
    ninter=nmsh-1;
    nmaxmesh = max(nmaxmesh,nmsh);
    if problem.debug, disp('Start 4th Order'); end
    iorder = 4;
    dc=zeros(problem.ncomp,ninter);
    if (linear)
        ludone = 0;
        [y,iflag,problem,rhs,fty,J,bhold,chold] = lineq_a(problem,h,t,nmsh,y,lambda,dc,ludone);
        fty = problem.f(t,y,lambda);
        if ~problem.vectorized
            problem.NFUN = problem.NFUN + length(t);
        else
            problem.NFUN = problem.NFUN + 1;
        end
    else
        rhsgiv = 0;
        [y,iflag,rhs,problem,fty,itnwt,J,bhold,chold]=newteq_a(problem,h,t,y,lambda,dc,ltol,tol,newton_par,J,bhold,chold);
        if(problem.iatt==-1), problem.nits = max(1,itnwt); end
        if problem.debug, fprintf('itnwt = %g\n',itnwt); end
       
    end
     J4=J;
    if iflag==0
        % the computation of the conditioning parameter is made only if the
        % order 8 converge, since for the moment is used only in the monitor
        % function

        
        [def6,Jflag,problem] = dfexcl_a(problem,h,nmsh,t,y,lambda,fty,ntol,ltol,tol,linear,A_num1,C_num1);
        
    else
        problem.iflbvp = 1;
        return
    end
    
    if (Jflag~=0)
        problem.iflbvp = 1;
        return
    end
    if problem.debug, disp('Start 6th Order'); end
    iorder = 6;
    def = def6;
    yold=y;
    if (linear)
        ludone = 1;
        [y,iflag,problem] = lineq_a(problem,h,t,nmsh,y,lambda,def,ludone,J4);
    else
        [y,iflag,problem,rhs,rnsq,fty]=fixjac_a(problem,iorder,ltol,tol,t,y,lambda,fty,def,def,rhs,J4);
        if (iflag==-3 && rnsq<fxfct*epsmch), rhsgiv = 1; end
        if (iflag~=0)
            [y,iflag,rhs,problem,fty]=newteq_a(problem,h,t,y,lambda,def,ltol,tol,newton_par);
        end
    end
    if (iflag~=0)
        problem.iflbvp = 1;
        return
    elseif (problem.ifinal==1)
        [errok,errmax6]=errest(ltol,tol,etest6,y,yold);
        if (errok)
            if comp_c
                condpar_old=condpar;
                [omega,condpar,stabcond] = comp_condpar(J4,problem.ncomp,nmsh,h,problem.a,problem.b,condpar_old,linear,compcond_par);
                if problem.debug
                    fprintf('sigma = %g gamma1 = %g kappa1 = %g kappa = %g\n',condpar.sigma,condpar.gamma1,condpar.kappa1,condpar.kappa);
                    fprintf('stab_kappa = %g stab_gamma1 = %g stab_kappa1 = %g stab_cond = %g stiff_cond%g  \n',stabcond.kappa,stabcond.gamma1,stabcond.kappa1,stabcond.stiff);
                end
            end
     
            problem.iflbvp = 0;
            err=errmax6;
            return
        end
    end
    
    if problem.debug, disp('Start 8th Order'); end
    iorder = 8;
    if linear
        fty = problem.f(t,y,lambda);
        if ~problem.vectorized
            problem.NFUN = problem.NFUN + length(t);
        else
            problem.NFUN = problem.NFUN + 1;
        end
    end
    [def8,Jflag,problem] = df8cal_a(problem,h,nmsh,t,y,lambda,fty,ntol,ltol,tol,linear,A_num2,B_num,C_num2);
    
    if (Jflag~=0)
        problem.iflbvp = 1;
        return
    end
    if linear
        def=def8;
    else
        def = def + def8;
    end
    
    yold=y;
    if (linear)
        ludone = 1;
        [y,iflag,problem] = lineq_a(problem,h,t,nmsh,y,lambda,def,ludone,J4);
    else
        rhsgiv = 0;
        [y,iflag,problem,rhs,rnsq,fty]=fixjac_a(problem,iorder,ltol,tol,t,y,lambda,fty,def,def8,rhs,J4);
        
        if (iflag==-3 && rnsq<fxfct*epsmch), rhsgiv = 1; end
        if (iflag~=0)
            [y,iflag,rhs,problem,fty]=newteq_a(problem,h,t,y,lambda,def,ltol,tol,newton_par);
        end
    end
    if (iflag==0)
        if problem.debug, disp('conv8'); end
        maxmsh = 0;
        % the computation of the conditioning parameter is made only if the
        % order 8 converge, since for the moment is used only in the monitor
        % function
         if comp_c
             condpar_old=condpar;
             [omega,condpar,stabcond] = comp_condpar(J4,problem.ncomp,nmsh,h,problem.a,problem.b,condpar_old,linear,compcond_par);
             if problem.debug
                 fprintf('sigma = %g gamma1 = %g kappa1 = %g kappa = %g\n',condpar.sigma,condpar.gamma1,condpar.kappa1,condpar.kappa);
                 fprintf('stab_kappa = %g stab_gamma1 = %g stab_kappa1 = %g stab_cond = %g stiff_cond%g  \n',stabcond.kappa,stabcond.gamma1,stabcond.kappa1,stabcond.stiff);
             end
         end
         
        if (problem.ifinal ==1)
            
            
            if (linear && strctr && nmsh < nmold), etest8=ones(size(tol)); end
            
            etest8=10*ones(size(tol));
            
            [errok,errmax8]=errest(ltol,tol,etest8,y,yold);
            if problem.debug, fprintf('errmax8 %g errok %g\n',errmax8,errok); end
            if errok
                problem.iflbvp=0;
                err=errmax8;
                return
            end
        end
        
        if (problem.ifinal == 1 && problem.iatt >= 1)
            
            [t,h,nmsh,maxmsh,told]=dblmsh(problem,t,nmax);
            yold = y;
            if (~maxmsh && problem.iprec ~= 2)
                y=interpu(told,yold,t);
            end
        else
            
            if use_c && stabcond.stiff 
                [t,h,nmsh,maxmsh,ddouble,told,nvold,voldmsh,phiold,problem]=selconderrmsh_a(problem,nvold,voldmsh,phiold,ltol,tol,problem.fixpnt,6,nmax,t,told,y,def8,omega,moncondmsh_par);
            else
                [t,h,nmsh,maxmsh,ddouble,told,nvold,voldmsh,phiold,problem]=selmsh_a(problem,nvold,voldmsh,phiold,ltol,tol,problem.fixpnt,6,nmax,t,told,y,def8);
            end
            if   (~maxmsh && problem.iprec ~= 2)
                if (linear && (problem.ifinal == 1 || problem.iatt ~= 0))
                    y=initu(problem,t,problem.ncomp);
                else
                    yold = y;
                    y=interpu(told,yold,t);
                end
            end
            
        end
        if (problem.iprec == 2)
            if (problem.debug), fprintf('Mesh Cannot Be Defined Within The Bounds Imposed'); end
            problem.iflbvp = 1;
            return
        end
    else
        problem.iflbvp = 1;
        return
    end
    
    if ~(maxmsh)
        problem.iflbvp = 0;
        if (problem.ifinal == 1)
            if (problem.iatt >= 1)
                if (nmsh < nmax/2)
                    inumb = inumb+1;
                else
                    problem.iflbvp = 1;
                    return;
                end
            end
            if (problem.iback ~= 1 && inumb == 3)
                problem.iflbvp = 1;
                return
            end
        elseif (problem.iatt == 0)
            
            return
        end
        if (problem.debug), fprintf('iatt %g\n',problem.iatt),end
        problem.iatt = problem.iatt+1;
        
    else
        problem.iflbvp = 1;
        return
    end
end
end


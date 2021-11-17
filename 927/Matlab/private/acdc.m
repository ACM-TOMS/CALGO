function [y,t,err,maxmsh,stabcond,condpar,problem,nc,nmaxmesh] = acdc(problem,nmsh,nmax,lambda,lambdamin,ltol,tol,maxcon)
%
%   Private function implementing the acdc and acdcc solvers based on the
%   fortran code acdc
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



ntol = size(ltol,1);
linear = problem.linear;
if problem.debug
    if (linear)
        fprintf('ncomp = %g\n',problem.ncomp);
    else
        fprintf('ncomp = %g\n',problem.ncomp);
    end
    fprintf('a = %g b = %g\n',problem.a,problem.b);
    nfxpnt = length(problem.fixpnt);
    if (nfxpnt > 0)
        fprintf('nfxpnt = %g fxpnt = %g\n',nfxpnt,problem.fixpnt(1:nfxpnt));
        fprintf('tolerance = %g\n',ltol(1:ntol));
        fprintf('tolerance = %g\n',tol(1:ntol));
    end
end
voldmsh=0;
phi = zeros(1,3);
e = zeros(1,3);
itsaim = 7;
epsmch=eps/2;
emin = 1/lambdamin;
ep = 1/lambda;

problem.pmax(1) = 0;
problem.pmax(2) = 0;
problem.nits = 0;
lambdap = 0;
lambda_changed = 0;
istep = 0;
nc = 0;
nss = 0;
idc = 0;
iextrap = 0;
istuk = 0;
ilin = -1;
hrat = 1;
amax = realmax;
bmax = 0;
cmax = realmax;
estep = realmax;
iusep = 0;
tolmax = 1;
phimax = 1/epsmch;
% Calculate Phimax, Which Is The Largest Monitor Function Value
% That We Believe Is Permissible In Double Precision For The Given Tolerance.
tolmax = min([tol(1:ntol)',tolmax]);
phimax = phimax*tolmax;
phiaim = 0;
phialt = 0;
epold = 0;
hsmlold = 0;
htpv = 0;
% We Predict That For Epsilon = 1 The Maximum Value Of The
% Monitor Function Is Minimised By Phi(3). We Will Use This Value
% When Performing Extrapolation.
e(3) = 1;
phi(3) = 1/(problem.b-problem.a);

if isscalar(nmsh)
    t=unimsh(problem.a,problem.b,nmsh,problem.fixpnt);
else
    t=nmsh;
end
y=initu(problem,t,problem.ncomp);
nmaxmesh = length(t);

while (problem.ifinal==0)
    lambda = max(1/ep,lambdamin);
    problem.lambda=lambda;
    ep = 1/lambda;
    if (lambda < (1.00001e0*lambdamin)),  problem.ifinal = 1;  end
    % Increase The Continuation Step Counter Nc.
    nc = nc+1;
    if (problem.debug),  fprintf('nc = %g lambda = %g\n',nc,lambda), end
    
    % If We Have Reached The Maximum Number Of Continuation Steps
    % Then This Will Be The Last Problem We Attempt.
    if (nc == maxcon)
        if (problem.debug), fprintf('maxcon = %g\n',maxcon), end
        problem.iback = 1;
        problem.ifinal = 1;
        lambda_changed=1;
    end
    
    problem.iatt = -1;
    problem.iflbvp = 0;
    problem.iucond = 0;
    
    
    [y,t,err,maxmsh,voldmsh,stabcond,condpar,nmaxmesh,problem] = bvpsol_a(problem,t,y,nmax,lambda,ltol,tol,voldmsh,nmaxmesh);
  
    if (problem.iflbvp == 0)
        problem.yc=y;
        problem.xc=t;
        
        nss = nss + 1;
        % If The Problem Eps = Epsmin Has Been Solved Succesfully Then
        % We May Finish. Ifinal = 1 When Eps = Epsmin.
        if (problem.ifinal == 1)
            if (problem.debug)
                fprintf('nc = %g lambda%g\n',nc,lambda);
                %     fprintf('y = %g ltol = %g\n',y, ltol(1:ntol));
                jstep = max(nmsh/30,1);
                for  i = 1:jstep:nmsh-1
                    fprintf('i = %g t = %g y = %g\n',i,t(i), y(ltol(1:ntol),i));
                end
                fprintf('nmsh = %g t(nmsh)= %g\n', length(t),t(length(t)));
                fprintf('y = %g\n',y(ltol(1:ntol),length(t)));
            end
            problem.lambda=lambda;
            if (lambda_changed), problem.iflbvp = -1; end
            return
        end
        % When Bactracking The Program Should Find A Problem That It Can
        % Eventually Solve. The Relationship Between The Last Two
        % Successful Epsilon Values Is Used To Restrict Future Epsilon
        % Values. This Is The Purpose Of Estep. If Istep = 1, Then We
        % Must Restrict Future Values Of Ep.
        if (istuk == 1)
            istep = 1;
            estep = min(ep/e(3),estep);
        end
        istuk = 0;
        % Calculate The Best Approximation To The Maximum Value Of The
        % Monitor Function By Extrapolating The Maximum Value On The
        % First Two Meshes. The Best Approximation Is Stored In Phit.
        
        if (iextrap == 0)
            h1 = problem.hord(1)/problem.pmax(1);
            h2 = problem.hord(2)/problem.pmax(2);
            if (h1/h2 > 1.1)
                c1h = (problem.pmax(2) - problem.pmax(1))/(h2-h1);
                phit = problem.pmax(1) - c1h*h1;
            else
                phit = problem.pmax(2);
            end
        end
        % It Is Possible That We Are Attempting To Solve Problems
        % That Are Beyond The Bounds Imposed By The Machine Precision.
        % The User Is Warned Of This Possibility.
        if (problem.iprec == 0 && phit > phimax )
            if (problem.debug), fprintf('lambda = %g\n',lambda), end
            lambdap = lambda;
            problem.iprec = 1;
        end
        % Save Details Of Last Problem Solved Successfully In Case We
        % Need To Backtrack
        
        told = t;
        yold = y;
        % If Iextrap Equals 1 Then We Have Abandoned Monitor Function
        % Extrapolation As Our Parameter Selection Basis.
        if (iextrap == 1)
            e(1) = e(2);
            e(2) = e(3);
            e(3) = ep;
            if (istep == 1)
                ep = estep*ep;
            else
                ep = emin;
            end
            continue
        end
        % If We Have A Decreasing Monitor Function Twice Consecutively
        % Then Extrapolation Will Not Work. The Flag Iextrap Equals 1
        % When Extrapolation Does Not Work.
        if (phit <= phi(3))
            if (idc == 1 || phit < phi(2))
                iextrap = 1;
                e(1) = e(2);
                e(2) = e(3);
                e(3) = ep;
                if (istep == 1)
                    ep = estep*ep;
                else
                    ep = emin;
                end
                continue
            else
                idc = 1;
                e(3) = ep;
                phi(3) = phit;
            end
            % Otherwise Update Extrapolation Data
        else
            idc = 0;
            e(1:2) = e(2:3);
            phi(1:2) = phi(2:3);
            e(3) = ep;
            phi(3) = phit;
        end
        % The Variable Irest Is A Flag Which Informs Us Whether Our
        % Extrapolation Procedure Is Accurate Or Not. Irest = 1 Means
        % Extrapolation Is Not Reliable And So We Should Restrict The
        % Next Parameter Step.
        irest = 0;
        % Check Errors In Extrapolation Procedure And Calculate Hmult
        % Which Restricts The Size Of Parameter Steps If Extrapolation
        % Is Working Poorly.
        if (nss == 2 && phit > 2*phiaim)
            irest = 1;
            hh = ep - epold;
            errp = abs(phit - phiaim);
            fx = errp/hh^2;
            hmax = sqrt(0.2*phiaim/fx);
            hmult = (epold+hmax)/epold;
        end
        
        if (nss > 2 && phit > 2*phiaim)
            irest = 1;
            hh = ep - epold;
            errp = abs(phit - phiaim);
            errp1 = abs(phit - phialt);
            if (ilin == -1)
                fx = errp/hh^2;
                hmax = sqrt(0.2*phiaim/fx);
                hmult = (epold+hmax)/epold;
                fx = errp1/hh^3;
                hmax = (0.2*phiaim/fx)^(1/3);
                hmult1 = (epold+hmax)/epold;
            else
                fx = errp/hh^3;
                hmax = (0.2*phiaim/fx)^(1/3);
                hmult = (epold+hmax)/epold;
                fx = errp1/hh^2;
                hmax = sqrt(0.2*phiaim/fx);
                hmult1 = (epold+hmax)/epold;
            end
        end
        epold = ep;
        % If The Extrapolation Procedure Is Particularly Inaccurate Then
        % We Need To Restrict It Permanently By Setting Istep = 1.
        pinacc = 3*max([problem.pmax(1),phiaim,phialt]);
        if (nss > 2 && phit > pinacc)
            istep = 1;
            pestep = max(hmult,hmult1);
            estep = min(estep,pestep);
        end
        % Decide Whether Linear Or Quadratic Extrapolation Is Best.
        if (nss > 2)
            per1 = abs(phiaim-phit);
            per2 = abs(phialt-phit);
            if (per2 < per1)
                ilin = -ilin;
                if (phit > 2*phialt);
                    hmult = hmult1;
                else
                    irest = 0;
                end
            end
        end
        % If The Continuation Steps Are Becoming Very Small Then We Have
        % Reached Our Final Problem
        dele = e(3)-e(2);
        if (dele < 0.01*e(3) || length(t) > (3*nmax/4))
            ep = e(3);
            lambdamin = 1/ep;
            if (problem.debug)
                if (dele < 0.01*e(3))
                    fprintf('Continuation Steps Too Small, Change lambdamin  To = %g\n',lambdamin);
                else
                    fprintf('Storage Limit Being Approached, Change lambdamin To = %g\n',lambdamin);
                end
            end
            emin = ep;
            problem.iback = 1;
            lambda_changed=1;
            continue
        end
        % The Following Section Of Code Calculates The Desired Value
        %(Phiaim) That We Would Like The Maximum Value Of The Monitor
        % Function To Take.
        itru = 0;
        if (phit > problem.pmax(1) && problem.pmax(2) ~= phit), itru = 1; end
        if (itru == 1), amax = -phit*phit/(2*c1h); end
        bmax = phit*h1;
        if (problem.npr(1) > 3*problem.npr(2)/2)
            ccm = max(bmax-3,bmax/1.5);
            cmax = min(cmax,ccm);
        end
        if (itru == 0), amax = max(1.5*bmax,bmax+3); end
        if (nss >= 2),  hrat = max(h1/hsmlold,1); end
        hsmlold = problem.hsml;
        bbm = max(1.5*bmax,bmax+3);
        if (~linear)
            fmax = itsaim*bmax/problem.nits;
            if (fmax > bmax),  fmax = max(fmax,bbm); end
            if (fmax < bmax),  fmax = max([bmax-3,bmax/1.5,fmax]);end
        end
        htot = min([amax,bbm,cmax]);
        if (linear), fmax = htot; end
        if (nss > 1 && htot ~= cmax)
            if (itru == 0 || amax > htpv), htot = max(htot,htpv); end
            htot = min([htot,fmax,cmax]);
        end
        htpv = htot;
        phiaim = htot/(problem.hsml*hrat);
        phiaim = max(1.05*phit,phiaim);
        % Quadratic And Linear Extrapolation To Find The Value Of Ep That
        % Corresponds To Phiaim
        if (nss ~= 1)
            f1 = (phi(2)-phi(1))/(e(2)-e(1));
            f2 = (phi(3)-phi(2))/(e(3)-e(2));
            ff2 = (f2-f1)/(e(3)-e(1));
            if (ff2 > 0)
                Qa = ff2;
                Qb = f1-(e(1)+e(2))*ff2;
                Qc = phi(1)-e(1)*(f1-e(2)*ff2)-phiaim;
                Qd = Qb^2-4*Qa*Qc;
                if (ilin == 1)
                    ep = min(emin,(-Qb+sqrt(Qd))/(2*Qa));
                    if (ep == emin), phiaim = phi(1)+(ep-e(1))*(f1+(ep-e(2))*ff2); end
                    phialt = phi(2)+f2*(ep-e(2));
                else
                    ep = min(emin,e(2)+(phiaim-phi(2))/f2);
                    if (ep == emin), phiaim = phi(2)+f2*(ep-e(2)); end
                    phialt = phi(1)+(ep-e(1))*(f1+(ep-e(2))*ff2);
                end
            else
                g1 = (e(2)-e(1))/(phi(2)-phi(1));
                g2 = (e(3)-e(2))/(phi(3)-phi(2));
                gg2 = (g2-g1)/(phi(3)-phi(1));
                if (ilin == 1)
                    ep = min(emin,e(1)+(phiaim-phi(1))*(g1+(phiaim-phi(2))*gg2));
                    if (ep == emin)
                        Qa = gg2;
                        Qb = g1-(phi(1)+phi(2))*gg2;
                        Qc = e(1)-phi(1)*(g1-phi(2)*gg2)-ep;
                        Qd = Qb^2-4.0*Qa*Qc;
                        phiaim = (-Qb+sqrt(Qd))/(2*Qa);
                    end
                    phialt = phi(2)+f2*(ep-e(2));
                else
                    ep = min(emin,e(2)+g2*(phiaim-phi(2)));
                    if (ep == emin), phiaim = phi(2)+f2*(ep-e(2)); end
                    Qa = gg2;
                    Qb = g1-(phi(1)+phi(2))*gg2;
                    Qc = e(1)-phi(1)*(g1-phi(2)*gg2)-ep;
                    Qd = Qb^2-4.0*Qa*Qc;
                    phialt = (-Qb+sqrt(Qd))/(2*Qa);
                end
            end
            % Linear Extrapolation
        else
            f2 = (phi(3)-phi(2))/(e(3)-e(2));
            ep = min(emin,e(2)+(phiaim-phi(2))/f2);
        end
        % Extrapolation May Not Be Working Very Well - If Not, Adjust
        % Continuation Stepsize And Recalculate Phiaim.
        if ((irest == 1 || istep == 1) && nss ~= 1)
            if (irest == 1 && istep == 1)
                ealt = min(epold*hmult,epold*estep);
            elseif (irest == 1)
                ealt = epold*hmult;
            else
                ealt = epold*estep;
            end
            if (ep > ealt)
                ep = min(emin,ealt);
                philin = phi(2)+f2*(ep-e(2));
                if (ff2 > 0)
                    phiqua = phi(1)+(ep-e(1))*(f1+(ep-e(2))*ff2);
                else
                    Qa = gg2;
                    Qb = g1-(phi(1)+phi(2))*gg2;
                    Qc = e(1)-phi(1)*(g1-phi(2)*gg2)-ep;
                    Qd = Qb^2-4.0*Qa*Qc;
                    phiqua = (-Qb+sqrt(Qd))/(2*Qa);
                end
                if (ilin == 1)
                    phiaim = phiqua;
                    phialt = philin;
                else
                    phiaim = philin;
                    phialt = phiqua;
                end
            end
        end
        % *****************************************************************
        % End Of Logic For A Successful Continuation Step And
        % Beginning Of Logic For An Unsuccessful Continuation Step.
        % *****************************************************************
    else
        % If A Problem Is Not Solved Successfully Then Print Details And
        % Select An Alternative Value For The Continuation Parameter.
        % This Process Is Known As Backtracking. Istuk = 1 Indicates That
        % We Are Backtracking.
        istuk = 1;
        % If Iback = 1 Then The Final Problem Has Not Been Solved.
        % In This Case We Stop.
        if (problem.iback == 1)
            if (problem.debug)
                if (lambdap ~= 0)
                    fprintf('lambda = %g iflbvp = %g lambdap = %g\n',lambda,problem.iflbvp,lambdap);
                else
                    fprintf('lambda = %g iflbvp = %g\n',lambda,problem.iflbvp);
                end
            end
            return
        end
        % If Iprec = 2, Then We Know That We Cannot Define A Mesh On
        % Which The Current Eps Value Can Be Solved To The Requested
        % Tolerances. We Alter Epsmin Accordingly.
        if (problem.iprec == 2)
            lambdamin = 1/(max((ep+e(3))/2,0.9*ep));
            emin = 1/lambdamin;
            lambda_changed = 1;
        end
        % Insert Details For Backtracking
        if (problem.debug), disp('Failed Step - Bactracking For Larger Value '); end
        problem.ifinal = 0;
        ep = (ep+e(3))/2;
        if ((ep-e(3)) < 0.01*e(3))
            ep = e(3);
            lambdamin = 1/ep;
            lambda_changed = 1;
            if (problem.debug), fprintf('Continuation Steps Too Small Change lambdamin To %g\n',lambdamin); end
            emin = ep;
            problem.iback = 1;
            problem.iflbvp = 0;
        end
        
        if (nss == 1)
            phiaim = phi(2)+f2*(ep-e(2));
        elseif (nss > 1)
            if (iextrap ~= 1)
                philin = phi(2)+f2*(ep-e(2));
                if (ff2 > 0)
                    phiqua = phi(1)+(ep-e(1))*(f1+(ep-e(2))*ff2);
                else
                    Qa = gg2;
                    Qb = g1-(phi(1)+phi(2))*gg2;
                    Qc = e(1)-phi(1)*(g1-phi(2)*gg2)-ep;
                    Qd = Qb^2-4.0*Qa*Qc;
                    phiqua = (-Qb+sqrt(Qd))/(2*Qa);
                end
                
                if (ilin == 1)
                    phiaim = phiqua;
                    phialt = philin;
                else
                    phiaim = philin;
                    phialt = phiqua;
                end
            end
        end
        
        % Re-Insert Details From Last Problem Successfully Solved.
        if (nss > 0)
            %if length(t) < length(told)
            %   t = told;
            %end
            if linear
               y=initu(problem,t,problem.ncomp);
            else
               y = interpu(told,yold,t);
            end
            problem.yc=y;
            problem.xc=t;
        end
           
        if (nss > 0)
            givmsh = 1;
            if (linear), Giveu = 1; end
        end
        % *****************************************************************
        % End Of Logic For An Unsuccessful Continuation Step.
        % *****************************************************************
     
    end
end
end

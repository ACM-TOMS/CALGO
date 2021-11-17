function [X, Y1, Y2, errmax, info] = APGR(G, K1, projMap, options, Y1)
% [X, Y1, Y2, errmax, info] = APGR(G, K1, projMap, options, Y1)
%*************************************************************************
% APGR Algorithm
% Reference:
% N. Ito, S. Kim, M. Kojima, A. Takeda, and K.-C. Toh.
% A Sparse Doubly Nonnegative Relaxation of Polynomial Optimization Problems
% with Binary, Box and Complementarity Constraints,
% Research Rport B-48?, Department of Mathematical and Computing Sciences, 
% Oh-Okayama, Meguro-ku, Tokyo 152-8552, March 2018. 
%*************************************************************************
% min { 0.5*norm(projK2(G+Y),'fro')^2 | Y in K1star} 
% K1 = cartesian product of psd cones
% K2 = (cartesian product nonnegative cones)
%      \cap (Null space of overlapping blocks)
%*************************************************************************
% Originally 
% [X, Y1, Y2, errmax, info] =...
%    projK1K2_fista_dualK(G, K1, projMap, options, Y1)
%*************************************************************************

    if (nargin < 4); options = []; end
    tol      = 1e-8;
    maxiter  = 300; 
    printyes = 3;
    if isfield(options, 'tol'); tol = options.tol; end
    if isfield(options, 'maxiterAPGR'); maxiter = options.maxiterAPGR; end
    if isfield(options, 'printyes'); printyes = options.printyes; end
    if ~isfield(projMap, 'K2')
        projMap.K2 = @(X) X; 
        projMap.K2Star = @(X) 0; 
    end
    symflag = checksymmetric(G, K1);
    if ~all(symflag)
       error('G is not symmetric'); 
    end
    normG = max(1, mexFnormK(G));
    reltol = tol * normG;
 
    %% Initialization
    tstart = clock; 
    if (printyes >= 2)
        fprintf('\n----------------------------------------------'); 
        fprintf('\n APGR Allgorithm');
        fprintf('\n----------------------------------------------');
        if (printyes >= 3) 
            fprintf('\n  iter   errmax   errratio   normX    Xerrratio normXratio'); 
            fprintf('     obj        objdiff'); 
            fprintf('\n  [reltol = %3.2e]', reltol); 
        end
    end
    if (nargin < 5) || isempty(Y1)
        Y1 = sparse(size(G, 1), size(G, 2)); 
    end   
    Y1bar = Y1; t = 1; mineigXratio = 1; 
    Linit = options.Linit;%0.8; 
    L = Linit; 
    Lreduce_iter = 0; 
    breakyes = 0; err1=0; err4=0; err5=0; err6=0; rankX = inf; obj = inf;
    rankY = [];
   
    %% FISTA
    restart_iter = 1;
    restart_interval = 100;
    recount = 0;
    for iter = 1:maxiter
        told = t; objold = obj; Lold = L; Y1old = Y1;  
       
        %% Step 1
        tau = 1 / L; 
        GpY1bar = G + Y1bar;
        DfYbar= projMap.K2(GpY1bar); % gradient
        Gbar  = Y1bar - tau * DfYbar; 

        % Use the Lanczos method if the rank is very small
        [Y1, ~, rankY] = blkprojSDPBP(Gbar, K1, rankY);  %%TKC

        %% Step 2
        t = 0.5*(1+sqrt(1+4*told^2));
        beta = (told-1)/t; 
       
        %% Step 3
        Y1bar = (1+beta)*Y1 - beta*Y1old;
       
        %% Y2 and Multiplier X
        GpY1 = G+Y1;
        Xtmp = DfYbar; 
        normY = mexFnormK(Y1); 
        err2tmp = abs(Xtmp(:)' * Y1(:)) / (1 + mexFnormK(Xtmp) + normY); 
        if false && (err2tmp > 20*reltol) % if KKT residual seems to be large
            Y2 = DfYbar - GpY1bar; % avoid computing projection to save computation time
        else
            Y2 = projMap.K2Star(-GpY1);
        end
        if (rem(iter, 10)==1)
            [~, Y2] = checksymmetric(Y2, K1);
        end
        X = GpY1 + Y2; %% important
        normX  = mexFnormK(X);
        normY2 = mexFnormK(Y2);
        obj = 0.5 * normX^2;
 
        %% Restart 
        objdiff = obj-objold;
        if objdiff > (10 * 1e-13) && (iter - restart_iter) > restart_interval
            restart_iter = iter;
            restart_interval = restart_interval * 2;
            t = 1; Y1 = Y1old; Y1bar = Y1old; obj = objold;
            recount = recount + 1;
            if (L < 1); L = min(1, L / 0.9); end
            if (printyes >= 2) 
                fprintf('\n      [restart : iter=%4.0f, objdiff=%- 3.2e, Lold=%4.3f, Lnew=%4.3f]',iter,objdiff,Lold,L);
            end
        elseif (objdiff < -1e-8) && (L > Linit && Linit < 1) && (iter - Lreduce_iter > 200)
            L = 0.97 * L;
            Lreduce_iter = iter;
            if (printyes >= 2) 
                fprintf('\n      [reduce L: iter=%4.0f, objdiff=%- 3.2e, Lold=%4.3f, Lnew=%4.3f]',iter,objdiff,Lold,L);
            end
        end

        %% Compute KKT residuals
        % important to update err4, err5 more frequently
        err2 = abs(X(:)' * Y1(:))  / (1+normX+normY); 
        err3 = abs(X(:)' * Y2(:)) / (1+normX+normY2); 
        err236 = max([err2,err3,err6]); 
        % If kkt residual (err236) is large, 
        % skip computing err4,5 which are expensive.
        if (err236 < 0.001*reltol) 
            compeigX = 1;
        elseif (err236 < 0.01*reltol) 
            compeigX = 2;  
        elseif (err236 < 0.1*reltol)
            compeigX = 5;
        elseif (err236 < 1.0*reltol)
            compeigX = 10;
        elseif (err236 < 10*reltol)
            compeigX = 20; 
        else
            compeigX = 50; 
        end
        recomperr45 = 0;
        if (rem(iter,compeigX)==0) 
            err5 = mexFnormK(projMap.K2Star(-X)) / (1 + normX);
            eigX = blkEigBP(X, K1); 
            err4 = norm(eigX(eigX < 0)) / (1 + normX); 
            mineigXratio = abs(min(eigX)) / normX; %%max(eigX);
            rankX = length(find(abs(eigX) > 0.1*sqrt(tol)));
            recomperr45 = 1; 
        end
        errmax = full(max([err1, err2, err3, err4, err5, err6])); 
        runhist.normX(iter) = normX;
        runhist.err(iter) = errmax;
        if (iter > 50)
            indx = (iter - 9):iter; 
            errratio = exp(mean(log(runhist.err(indx) ./ runhist.err(indx - 30))));
            normXratio = exp(mean(log(runhist.normX(indx) ./ runhist.normX(indx - 30)))); 
        else
            errratio = 1; normXratio = 1; 
        end
        Xerrratio = normX / max(1e-20, errmax);
        
        %% Stopping criteria
        if (iter > 50) && (errmax < reltol) && (Xerrratio > 300) ...
          && (normX > 100*reltol) && (normXratio > 0.9) && options.heuristicFISTA
            breakyes = 1; % infeasible
        elseif (normX < 0.5*reltol) && (errmax < reltol) 
            breakyes = 2; % feasible
        elseif (normX < 1.0*reltol) && (errmax < 0.01*reltol)
            breakyes = 3; % feasible
        elseif (normX < 2.0*reltol) && (errmax < 0.1*reltol) ...
          && (rankX==0) && (obj < 0.01*tol) 
            breakyes = 4; % feasible
        elseif (errmax < 1e-4*reltol)
            breakyes = 5; % infeasible
        end
        
        %% Heuristic stopping criteria (infeasible)
        if options.heuristicFISTA
            if toc(options.startingTime) > options.maxtime
                breakyes = 79;
            end
            if (iter > 300) && (errmax < sqrt(reltol)) && (Xerrratio > 1e4) ...
              && (normX > 100*sqrt(reltol)) && (normXratio > 0.995) 
                breakyes = 88; %%new
            end
%             if (iter > 300) && (errmax < 1000*reltol) && (Xerrratio > 2e3) ...
%               && (normX > 5*10000*reltol) && (normXratio > 0.997) ...
%               && (errratio > 0.95)
%                 %%breakyes = 89; %%new
%             end
            if (iter > 1*100) && (errmax < 100*reltol) && (Xerrratio > 5e3)...
              && (normX > 10000*reltol) && (normXratio > 0.995) 
                breakyes = 90; %%new
            end
            if (iter > 1*100) && (errmax < 10*reltol) && (Xerrratio > 500)...
              && (normX > 10*reltol) && (normXratio > 0.995) ...
              && (errratio > 0.5)
                breakyes = 91;
            end
            if (iter > 1*100) && (errmax < 100*reltol) && (Xerrratio > 500)...
              && (normX > 100*reltol) && (normXratio > 0.995) ...
              && (errratio > 0.5)
                breakyes = 92;
            end
            if (iter > 300) && (errmax < 0.1*reltol) && (Xerrratio > 500)...
              && (normX > 10*reltol) && (normXratio > 0.99) ...
              && (errratio > 0.5)
                breakyes = 93;
            end
            if (iter > 300) && (errmax < 0.01*reltol) && (Xerrratio > 500)...
              && (normX > 1*reltol) && (normXratio > 0.995) ...
              && (errratio > 0.5)
                breakyes = 94;
            end
            if (iter > 500) && (errmax < 0.1*reltol) && (Xerrratio > 300)...
              && (normX > 10*reltol) && (normXratio > 0.98) ...
              && (errratio > 0.8) && (mineigXratio < 1e-3)
                breakyes = 95; 
            end
            if (iter > 500) && (errmax < 10*reltol) && (Xerrratio > 500)...
              && (normX > 100*reltol) && (normXratio > 0.97) ...
              && (errratio > 0.8) && (mineigXratio < 5e-4)
                breakyes = 96; 
            end
            if (iter > 300) && (errmax < 2e4*reltol) && (Xerrratio > 400)...
              && (normX > 1e6*reltol) && (normXratio > 0.999) ...
              && (errratio > 0.999)
                breakyes = 97; %new
            end 
            if (iter > 300) && (errmax < 1e3*reltol) && (Xerrratio > 500)...
              && (normX > 5*reltol) && (normXratio > 0.999) ...
              && (errratio > 0.99) && (mineigXratio < 1e-4)
                breakyes = 99;
            end
        end
        
        %% Print info
        if (breakyes > 80) && (printyes >= 2)
            fprintf('\n      [stagnate: %2.0f, errratio=%5.4f, normXratio=%8.7f]',breakyes,errratio,normXratio);           
        end
        if (iter <= 200) 
            print_iter = 50; 
        elseif (iter <= 2000)
            print_iter = 100; 
        else
            print_iter = 200; 
        end   
        if (rem(iter,print_iter)==1 || breakyes || iter==maxiter) && (printyes >= 3)
            fprintf('\n  %4.1d| %3.2e',iter,errmax); 
            fprintf(' %3.2e %3.2e %3.2e %8.7f|',errratio,normX,Xerrratio,normXratio);
            if (iter > 1)
                fprintf(' %7.6e %- 3.2e| L=%3.2f|',obj,objdiff,L);
            end
            if (recomperr45) && (err4>0)
                fprintf(' [mineigXratio=%3.2e]',mineigXratio);
            end
        end
        
        %%
        if (breakyes); break; end
    end 
    [eigX] = blkEigBP(X,K1); 
    [Y1,eigY] = blkprojSDPBP(Y1,K1); normY = norm(eigY); 
    Y2 = projMap.K2Star(-G-Y1);  %%new
    err = zeros(1,7); 
    err(1) = mexFnormK(X-G-Y1-Y2)/(1+normX);
    err(2) = abs(sum(sum(X.*Y1)))/(1+normX+normY);
    err(3) = abs(sum(sum(X.*Y2)))/(1+normX+normY2); 
    err(4) = norm(eigX(eigX < 0))/(1+normX);
    err(5) = mexFnormK(projMap.K2Star(-X))/(1+normX);
    err(6) = norm(eigY(eigY < 0))/(1+normY); 
    err(7) = mexFnormK(projMap.K2(-Y2))/(1+normY2);
    ttime = etime(clock,tstart);
    if (printyes >= 2)
        rankX = length(find(eigX > 0.1*sqrt(tol)));
        mineigX = min(eigX); maxeigX = max(eigX); 
        Xmin  = full(min(min(X))); Xmax = full(max(max(X)));
        fprintf('\n----------------------------------------------');
        fprintf('\n  number iter = %2.0d,  time = %3.2f',iter,ttime);
        fprintf('\n  0.5*norm(X-B)^2 = %9.8e',0.5*norm(X-G,'fro')^2);
        fprintf('\n  normB = %2.1e, B(1,1) = %2.1e',mexFnormK(G),full(G(1,1))); 
        fprintf('\n  min(eigX) = %3.2e, max(eigX) = %3.2e,',mineigX,maxeigX); 
        fprintf(' rankX = %1.0f, rankY = %1.0f',rankX,sum(rankY));
        fprintf('\n  min(X) = %3.2e, max(X) = %3.2e',Xmin,Xmax);
        fprintf('\n  relnormX = %3.2e',normX/normG);
        fprintf('\n  reltol = %3.2e',reltol);
        fprintf('\n  relerr = %3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e,%3.2e',err/normG);
        fprintf('\n  relative residual in KKT conditions = %3.2e',max(err)/normG);
        fprintf('\n  breakyes = %2.0f',breakyes);
        fprintf('\n  ||X||/||B-X|| = %3.2e',normX/mexFnormK(G-X));
        fprintf('\n  ||G+Y1+Y2||/||G|| = %3.2e',mexFnormK(G+Y1+Y2)/normG);
        fprintf('\n  ||G||,||X||,||Y||,||Y2|| = %3.2e,%3.2e,%3.2e,%3.2e',normG,normX,normY,normY2); 
        fprintf('\n----------------------------------------------\n');
    end
    info.iter = iter;
    info.breakyes = breakyes;
    info.recount = recount;
% if (false)
%    subplot(121); eigBpY = blkEigBP(BpY,K1); plot(eigBpY(2:end)); 
%    rankBpY = length(find(eigBpY > -1e-16)); 
%    fprintf('\n rank(BpY) = %3.0f, n=%3.0f',rankBpY,sum(K1.s)); 
%    axis('square'); grid;
%    subplot(122); sortY2 = sort(Y2(:)); plot(sortY2(1:end-1));
%    len = length(find(Y2(:) > -1e-16)); 
%    fprintf('\n percentage of nz(Y) = %3.2f',len/sum(K1.s)^2); 
%    axis('square'); grid;
%    pause
% end
%     if options.plot
%         semilogy(1:iter, runhist.normX(1:iter), '-', 'LineWidth', 2, 'MarkerSize', 10);
%         hold on
%         semilogy(1:iter, runhist.err(1:iter), '--', 'LineWidth', 2, 'MarkerSize', 10);
%         hold off
%         legend('normX','errmax','Location', 'best');
%         set(gca, 'FontSize', 20);
%         keyboard
%     end
end%function 

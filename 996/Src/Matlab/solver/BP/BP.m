function [sol, info] ...
    = BP(Q0Vect, H0Vect, K1, projMap, LB, UB, params)
%BP Algorithm to solve 
%
%   maximize    y0
%   subject to  Q0Vect - y0 * H0Vect = Y1 + Y2
%               Y1 \in K1star,  Y2 \in K2star
%
%   maximize    y0 + mu * UbdIX
%   subject to  Q0Vect - y0 * H0Vect - mu * I = Y1 + Y2
%               Y1 \in K1star,  Y2 \in K2star,  mu \leq 0
% 
% Usage:
%   [sol, info] = BP(Q0Vect, H0Vect, K1, projMap, LB, UB, params);
% 
% Reference:
% N. Ito, S. Kim, M. Kojima, A. Takeda, and K.-C. Toh.
% A Sparse Doubly Nonnegative Relaxation of Polynomial Optimization Problems
% with Binary, Box and Complementarity Constraints,
% Research Rport B-48?, Department of Mathematical and Computing Sciences, 
% Oh-Okayama, Meguro-ku, Tokyo 152-8552, March 2018. 
    
    %% check input data
    errorSW = 0;
    if nargin < 7
        errorSW = errorSW + 1;
        fprintf('Some of input QlamVect, H0Vect, K1, K2, LB, UB, params is missing?\n');
    end
    
    function out = toScalar(in)
        if ~isscalar(in)
            errorSW = errorSW + 1;   
            out = in;
        else
            out = full(in);
        end
    end
    LB = toScalar(LB);
    UB = toScalar(UB);
    
    matSize2 = sum(K1.s .* K1.s);
    function Xout = checkVect(Xin)
        Xin = Xin(:)';
        if length(Xin) ~= matSize2
            errorSW = errorSW + 1;
            fprintf('length(Xin) = %d ~= matSize*matSize = %d\n', size(Xin, 2), matSize2);
            Xout = Xin;
        else
            [symflag, Xout] = checksymmetric(Xin, K1);
            if ~all(symflag)
                errorSW = errorSW + 1;
                fprintf('Xin is not symmetric\n');
            end
            Xout = full(Xout);
        end
    end
    Q0Vect = checkVect(Q0Vect);
    H0Vect = checkVect(H0Vect);

    if errorSW > 1
        fprintf('\n!!! Input data QlamVect, H0Vect, matSize, params are inconsistent !!!\n');
        fprintf('!!! At least %d errors are found !!!', errorSW);
        error('Error at BP_');
    end
    
    %% fill unspecified params
    params = fillUnspecifiedParams(params, defaultParamBP());

    %% initialization of LB, UB, and y0
    LBv = -1e21;
    y0init = UB;
    
    %% scaling the given COP
    % This scaling is fixed throughout the BP Algorithm
    if params.scale_data == 1
        Q0Vectscale = max(1e-3, norm(Q0Vect));
    elseif params.scale_data == 2
        Q0Vectscale = max(1e-3, sqrt(norm(Q0Vect)));
    else
        Q0Vectscale = 1;
    end
    Q0Vectscale = full(Q0Vectscale);

    Q0Vect = Q0Vect / Q0Vectscale;
    sol.y0init = y0init / Q0Vectscale;
    sol.UB = UB / Q0Vectscale;
    sol.LB = LB / Q0Vectscale;
    sol.LBv = LBv / Q0Vectscale;
    
    if params.printyes >= 2
        fprintf(' Original: LB=%+7.6e, y0init=%+7.6e, UB=%+7.6e\n', LB, y0init, UB);
        fprintf(' Scaled  : LB=%+7.6e, y0init=%+7.6e, UB=%+7.6e\n', sol.LB, sol.y0init, sol.UB);
        fprintf(' Q0Vectscale=%7.6e\n', Q0Vectscale);
    end
    
    %% Allocate
    runhist = zeros(10, 1 + params.maxiterBP);  
    runhist(:, 1) = [sol.y0init; sol.LBv; sol.LB; sol.UB; sol.UB-sol.LB; inf; inf; -1; 0; 0];
    
    %% Initialize
    y0 = sol.y0init;
    G = y0 * H0Vect - Q0Vect;
    Y1 = projMap.K1Star(-G);
%     Yu = Y1;
%     Yl = zeros(size(Y1));
    
    %% Gscale --- scaling regression problem
    % This scaling is updated at each iteration of the BP Algorithm
    Gscale = 1;
    if params.Gscale_yes
        Gscale = max(1e-3, sqrt(norm(G, inf))); % Shrink large G
        G = G / Gscale;
        Y1 = Y1 / Gscale;
%         Yu = Yu / Gscale;
%         Yl = Yl / Gscale;
    end  


%%
    if params.printyes >= 1
        fprintf('\n----------------------------------------------\n'); 
        fprintf(' BP Algorithm');
        fprintf('\n----------------------------------------------'); 
    end
    startingTime = tic;
    params.startingTime = startingTime;
    break_yes = 0;
%     l_updated = false;
%     updated_enough = false;
%     Rcount = 0;
    recount = 0;
    APGiter = 0;
    for iteration = 1:params.maxiterBP
        %% Solve regression
        %if paramBP.newBPswitch && paramBP.betterInit; Y1 = Yu; end
        [X, Y1, Y2, ~, infoProj] = APGR(G, K1, projMap, params, Y1);
%         Rcount = Rcount + infoProj.Rcount;
        recount = recount + infoProj.recount;
        APGiter = APGiter + infoProj.iter;
        if params.printyes >= 2; fprintf('\n'); end
        
        %% Computing sol.LBv (valid lower bound) 
        validy0 = 0;
        % validy0 = 0 ---> dual infeasible (no update sol.LB), no update sol.LBv
        % validy0 = 1 ---> dual feasible (update sol.LB), no update sol.LBv
        % validy0 = 2 ---> dual infeasible (no update sol.LB), update sol.LBv
        % validy0 = 3 ---> dual feasible (update sol.LB), update sol.LBv
%         rhoYeigMin = inf;
%         if params.validSW || (sol.UB-sol.LB > 10*params.delta1*max([1e-10, abs(sol.UB), abs(sol.LB)]))
            YeigMin = min(min(blkEigBP(-G-Y2, K1)), 0);
            rhoYeigMin = Gscale*params.UbdIX*YeigMin;
            newLBv = y0 + rhoYeigMin;
            updated_enough = (newLBv > (sol.LB + 0.1 * (sol.UB - sol.LB)));
            if newLBv > sol.LBv
                validy0 = 2;
                sol.LBv = newLBv;
                if params.printyes >= 2
                    fprintf(' updated LBv= %+7.6e\n', sol.LB);
                end
            end
%         end

        %% Check feasibility
        relnormX = norm(X, 'fro') / max(1e-3, norm(G, 'fro')); % Does the scaling work?
        if (relnormX < params.validTOL) % Does it work when breakyes = 4?
%             l_updated = true;
%             Yl = Y1;
%            sol.y0 = y0; 
            sol.LB = max(y0, sol.LBv);
            if params.printyes >= 2
                fprintf(' updated LB = %+7.6e\n', sol.LB);
                fprintf(' y0=%+7.6e is valid, relnormX= %3.2e, rhoYeigMin= %3.2e\n',y0,relnormX,rhoYeigMin);
            end
            validy0 = validy0 + 1;
        else
%             Yu = Y1;
            sol.LB = max(sol.LB, sol.LBv);
%             if params.validSW && params.stabilize && updated_enough
            if params.stabilize && updated_enough
                if params.printyes >= 2
                    fprintf(' y0=%+7.6e may be invalid, relnormX= %3.2e, rhoYeigMin= %3.2e\n',y0,relnormX,rhoYeigMin);
                    fprintf(' But skipped updating UB.');
                end
            else
                sol.UB = y0;
                if params.printyes >= 2
                    fprintf(' updated UB= %+7.6e\n', sol.UB);
                    fprintf(' y0=%+7.6e is invalid, relnormX= %3.2e, rhoYeigMin= %3.2e\n',y0,relnormX,rhoYeigMin);
                end
            end
        end
        
        %% Store and print informations
        runhist(:, 1+iteration) = ...
            [y0; sol.LBv; sol.LB; sol.UB; sol.UB-sol.LB; relnormX; rhoYeigMin; validy0; infoProj.iter; infoProj.breakyes];
        if params.printyes == 1
            fprintf('\niter=%2d ', iteration);
        end    
        if params.printyes >= 2
            fprintf('\nOriginal\niter=%2d:y0=%+7.6e,LBv=%+7.6e,[LB,UB]=[%+7.6e,%+7.6e]\n',...
                iteration, y0 * Q0Vectscale, sol.LBv * Q0Vectscale, ...
                sol.LB * Q0Vectscale, sol.UB * Q0Vectscale); 
            fprintf('Scaled\niter=%2d:y0=%+7.6e,LBv=%+7.6e,[LB,UB]=[%+7.6e,%+7.6e]\n',...
                iteration, y0,sol.LBv, sol.LB, sol.UB);
        end

        %% Stopping criteria
        if (sol.UB - sol.LB) * Q0Vectscale < params.delta
            break_yes = 1; % absolute error is small
            break; 
        elseif (sol.UB - sol.LB) < params.delta1 * max([1e-10, abs(sol.UB), abs(sol.LB)]) 
            break_yes = 2; % relative error is small
            break; 
        % elseif (params.validSW == 1) && (runhist(1, iteration) < runhist(4, iteration)) ...
        elseif (runhist(1, iteration) < runhist(4, iteration)) ...
                && ((sol.LBv - runhist(2, iteration)) < params.delta2 * max([1e-10, abs(sol.UB), abs(sol.LBv)])) ...
                && ((runhist(4, iteration) - sol.UB) < params.delta2 * max([1e-10, abs(sol.UB), abs(sol.LBv)]))
            break_yes = 3; % convergence slows down
            break; 
        end
        if toc(params.startingTime) > params.maxtimeBP
            break_yes = -1;
            break;
        end

        %% update solution and scale
        y0 = (1 - params.weightLBvsUB) * sol.LB + params.weightLBvsUB * max(y0, sol.UB);
        G = y0 * H0Vect - Q0Vect;
        GscaleOld = Gscale;
        if params.Gscale_yes
            Gscale = max(1e-3, sqrt(norm(G, inf))); %max(1, sqrt(norm(G, inf)));
            G = G / Gscale;
            Y1 = (GscaleOld / Gscale) * Y1;
%             Yl = (GscaleOld / Gscale) * Yl;
%             Yu = (GscaleOld / Gscale) * Yu;
        end
%         if params.newBPswitch && l_updated
%             params.gammaR = params.gamma * norm(Yl - Yu)^2;
%         end
    end

    timeBisection = toc(startingTime);
    
    %% Reallocate
    runhist = runhist(:, 1:1+iteration);

    %% Print information
    runhist(1:5, 1:1+iteration) = runhist(1:5, 1:1+iteration) * Q0Vectscale;
    fprintf('\n     LBv           y0            UB           UB-LBv   relnormX iter b_yes');
        fprintf('\n%2d  %+7.6e %+7.6e %+7.6e %3.2e ', 0, runhist([2, 1, 4, 5], 1));
    for j = 2:iteration+1 % size(runhist,1)
        % fprintf('\n%2d  %+7.6e %+7.6e %+7.6e %3.2e %3.2e  %d     %4d %d',j-1,runhist([2,1,4,5,6,8,9,10],j));
        fprintf('\n%2d  %+7.6e %+7.6e %+7.6e %3.2e %3.2e %4d %d',j-1,runhist([2,1,4,5,6,9,10],j));
    end
    fprintf('\n timeBP (Excecution time) = %4.2f,  termcodeBP = %d\n',...
        timeBisection, break_yes);

    %% Rescale and output
    sol.LB = Q0Vectscale * sol.LB;
    sol.LBv = Q0Vectscale * sol.LBv;
    sol.UB = Q0Vectscale * sol.UB;
    sol.Y1 = Y1(:)*Gscale*Q0Vectscale;  
    sol.Y2 = Y2(:)*Gscale*Q0Vectscale;  
    info.recount = recount;
    info.iter = iteration;
    info.Q0Vectscale = Q0Vectscale;
    info.Gscale = Gscale;  
    info.APGiter = APGiter;
    info.params = params;
    info.runhist = runhist;
    info.break_yes = break_yes;
    info.timeBisection = timeBisection;
end


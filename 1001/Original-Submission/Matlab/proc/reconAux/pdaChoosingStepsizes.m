%% pdaChoosingStepsizes
%
% Choose the primal step size tau and the dual step size sigma for
% primal-dual algorithm.
%
%% Syntax
%
%   [tau,sigma] = pdaChoosingStepsizes(stepsize,L,ipda,xnpRVDstore,xnRVDstore,KcompNorm,tau,sigma,vartheta,seti,JA,JB)
%
%% Description
% |[tau,sigma] = pdaChoosingStepsizes(stepsize,L,ipda,xnpRVDstore,xnRVDstore,KcompNorm,tau,sigma,vartheta,seti,JA,JB)|
% computes the  primal step size tau and the dual step size sigma for
% primal-dual algorithm in dependence of the method in |stepsize|.
% This is a internal function and used in pda.m
%
%
%% Input Arguments
%
% * stepsize    :   Method to compute the step sizes
%                   ('fix' (default) or 'adaptive')
%                   (stored in seti.pdaStepsize)
% * L           :   Operator norm L = ||K||, see <opNormNum.html>
% * ipda        :   iteration number of inner iteration (pda iteration)
% * xnpRVDstore :   x_{n+1} as real vector downscaled
%                   (real i.e. via transformation from complex to real).
% * KcompNorm   :   components of K to compute the operator norm L = ||K|| in pda
%                   or to compute the adaptive step sizes.
%                   See subfunction KcomponentsStruct in <pda.html> for
%                   details.
% * tau         :   primal step size tau from iteration before
% * sigma       :   dual step size sigma from iteration before
% * vartheta    :   vartheta $\in (0,1)$ (vartheta = seti.vartheta in pda.m)
% * seti        :   struct with seti.dVinv for norminv2, see <norminv2.html>
% * JA and JB   :   auxilary matrices for Jacobian matrix, see <mimo.html>
%
%
%% Output Arguments
%
% * tau         :   primal step size
% * sigma       :   dual step size
%
%% More About
%
% See also Section 4.8 in [2].
%
% *Fix step size*
%
% $\tau = 1/L$
%
% $\sigma = \tau$
%
% *Adaptive step size*
%
% This method is from p. 6 in [3].
% They discussed a sufficient estimation in the proof of convergence in 
% Theorem 1 in [1] that leads to the following choose.
%
% * $\chi := \|x_d\|_{\mathrm{roi},2}/K_\mathrm{comp}(x_d)$ 
%   with difference $x_d := x^{(n)} - x^{(n-1)}$
%
% Code:
%
%   chi = norminv2(xdRVD,seti)/KcompNorm(xdRVD,JA,JB);
%
% Note the relation to the operator norm in opNormNum
%   
%   opnorm1(i) = KcompNorm(h,JA,JB)/norminv2(h,seti);
%
%
% * $\rho := \sigma \tau$
%
% * $\texttt{scriptS} = \bf{P}(\chi)$ with cases
%
% $\bf{P}(\chi) := \chi \quad \mathrm{ if\ } \chi \leq \vartheta \rho$,
%
% $\bf{P}(\chi) := \sqrt{\vartheta \rho} 
%  \quad \mathrm{ if\ } \vartheta \rho < \chi \leq \rho$,
%
% $\bf{P}(\chi) := \sqrt{\rho} \quad \mathrm{ if\ } \rho < \chi.$
%
% * Then we get as new step sizes 
%  with balancing parameter $\eta$ (currently $\eta = 1$)
%
% $\sigma = \texttt{scriptS} \cdot \eta$
%
% $\tau = \texttt{scriptS}/ \eta$
%
%
%
%% References
%
% * [1] Antonin Chambolle and Thomas Pock. A first-order primal-dual algorithm for convex problems with applications to imaging. _Journal of Mathematical Imaging and Vision_, 40(1):120-145, 2011.
% * [2] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
% * [3] Kristian Bredies and Martin Holler.
%  A TGV-based framework for variational image decompression, zooming, and reconstruction. Part II: Numerics. 
%  _SIAM Journal on Imaging Sciences_, 8(4):2851-2886, 2015.
%
%% See Also
%
% * <opNormNum.html>
% * <pda.html>
%
%% Code
%
function [tau,sigma] = pdaChoosingStepsizes(stepsize,L,ipda,xnpRVDstore,xnRVDstore,KcompNorm,tau,sigma,vartheta,seti,JA,JB)

switch stepsize
%%
% *method: fix step size*
%
case 'fix'
    tau = 1/L;
    sigma = tau;

    % test 1:
    % tau = 1;
    % sigma = 1/L^2;
    
    % test 2:
    % factor = 1/1E-5; % strong sparsity influence...
    % tau = 1/(L*sqrt(factor)); sigma = factor*tau;
    
    % fprintf('   stepsizes: tau = %g, sigma = %g\n',tau,sigma);

%%
% *method: adaptive step size*
%
case 'adaptive'
    if ipda < 3 % start value
        tau = 1/L;
        sigma = tau;
    else
        % parameter set:
        eta = 1; % balancing parameter
        % p. 6 (BrediesHoller_2015B.pdf): 
        % eta balances primal stepsize tau and dual step size sigma;
        % in application: choose eta = 1
        
        xdRVD = xnpRVDstore-xnRVDstore; % x difference: x^n - x^(n-1) (here: x^n = xnp = x^(n+1) and x^(n-1) = xn = x^n)
        % Kristian Bredies and Martin Holler: A TGV-based framework for variational image decompression, zooming and reconstruction. Part II: Numerics
        % eq. (9)

        chi = norminv2(xdRVD,seti)/KcompNorm(xdRVD,JA,JB);
        % compare with opNormNum: opnorm1(i) = KcompNorm(h,JA,JB)/normroi2(h,seti);
        fprintf('chi = %g | ',chi);
        
        % sigmatau = sigma * tau
        % sigma_n = sigman, sigma_{n+1} = sigmanp
        % notation does not fit to xn and xnp (it's shifted...)

        sigmataun = sigma*tau;
        fprintf('sigmataun = %g | ',sigmataun);
      
        %scriptS(sigma*tau,chi)
        if chi <= vartheta*sigmataun % vartheta = 1 was set
            scriptS = chi;
            fprintf('NOT as fix stepsize | ');
            % error('stop'); % ------------------------------------
        elseif chi > sigmataun
            scriptS = sqrt(sigmataun); % as fix stepsize...
        else % vartheta*sigmataun < chi =< sigmataun % because vartheta = 1 this case will never be used
            scriptS = sqrt(vartheta*sigmataun);
            fprintf('NOT as fix stepsize | ');
            % error('stop'); % -----------------------------------
        end

        sigmanp = scriptS*eta;
        taunp = scriptS/eta;

        % now use sigmanp and taunp
        sigma = sigmanp;
        tau = taunp;
        fprintf('tau = %g, sigma = %g\n',tau,sigma);
    end
end

end

%% ADFFqFast
% Fast computation of the adjoint of the derivative if Jacobi
% matrices JA and JB were already computed.
%
%% Syntax
%
%   ADFFq = ADFFqFast(FFqmF,JA,JB,seti)
%
%% Description
%
% |ADFFq = ADFFqFast(FFqmF,JA,JB,seti)| is a fast computation of the
% adjoint of the derivative, ADFFq, if the Jacobi matrices JA and JB were 
% already computed.
%
% |ADFFq| gives the same as |ADFFq| in |[FFqmF, ADFFq] = mimo(qROI,FmeasDelta)|.
%
%% Input Arguments
%
% The following input arguments are described in <mimo.html>:
%
% * FFqmF
% * JA, JB
% * seti.model
% * seti.dSMeas
% * seti.dV
%
%% Output Arguments
%
% * ADFFq : adjoint of derivative 
%           (complex matrix written as vector of size seti.nROI^seti.dim x 1)
%           (see also <mimo.html>).
%
%% More About
%
% $\texttt{ADFFq} = [\mathcal{F}'(q)]^\ast [\mathcal{F}(q) -
% F_\mathrm{meas}^\delta]$
%   is the adjoint of derivative.
%
% For fast computation use, see [1], Sec. 3.6, eq. (35):
%
% $[\underline{\mathcal{F}}'(\underline{q})]^\ast \underline{H}$
% $= \sum_{j = 1}^{N_{\mathrm{Sca}}} \frac{\omega_j^{\mathrm{Sca}}}{\texttt{seti.dV}}$
% $\sum_{\ell = 1}^{N_{\mathrm{Inc}}}
%  \underline{H}_{j,\ell}\, \overline{A_{N_D,N_{\mathrm{Sca}}}(j,\cdot)} \, 
%  \overline{B_{N_D,N_{\mathrm{Inc}}}(\cdot,\ell)}$
% $\quad \mathrm{for\ } \underline{H} \in \bf{C}^{N_{\mathrm{Sca}} \times N_{\mathrm{Inc}}}$
%
% with 
%
% * $N_{\mathrm{Inc}} = \texttt{seti.incNb}$,
% * $N_{\mathrm{Sca}} = \texttt{seti.measNb}$,
% * $\omega^{\mathrm{Sca}} = \texttt{seti.dSMeas}$,
% * $N_D = \texttt{seti.nROI}^\texttt{seti.dim}$.
%
%% References
%
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
%
% * <mimo.html>
%
%% Code
%
function ADFFq = ADFFqFast(FFqmF,JA,JB,seti)

DG2 = zeros(size(JA,2),1);
% ADFFq has the same size as qROI as vector

if strcmp(seti.model,'helmholtz2D') || strcmp(seti.model,'helmholtz3D')
    for iRX = 1:size(JA,1)
        for iTX = 1:size(JB,2)
            DG2 = DG2 + FFqmF(iRX,iTX)*conj(JA(iRX,:).'.*JB(:,iTX));
        end
    end
%     DGAux = JA'*FFqmF*JB';
%     disp('Test in DFFfast: ')
%     norm(DG2-DGAux)/norm(DG2)
elseif strcmp(seti.model,'helmholtzHMode2D')
    for iRX = 1:size(JA(:,:,1),1)
        for iTX = 1:size(JB(:,:,1),2)
            DG2 = DG2 + FFqmF(iRX,iTX)*conj(JA(iRX,:,1).'.*JB(:,iTX,1))...
                      + FFqmF(iRX,iTX)*conj(JA(iRX,:,2).'.*JB(:,iTX,2));
        end
    end
else
    %fprintf(strcat('Error in ADFFqFast - pda not implemented for model ',  seti.model))
    fprintf('Error in DFFfast - pda not implemented for model %s.\n',  seti.model)
end

DG2 = DG2*seti.dSMeas/seti.dV;

ADFFq = DG2;

end

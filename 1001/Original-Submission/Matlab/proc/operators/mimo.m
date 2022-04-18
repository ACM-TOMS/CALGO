%% mimo
% multiple-input-multiple-output
%
%
%% Syntax
%
%   [FFqmF,ADFFq] = mimo(seti,qROI,FmeasDelta)
%   [FFqmF,ADFFq] = mimo(seti,qROI,FmeasDelta,'adjOfDer')
%   [FFqMeas,~,FFqROI] = mimo(seti,qROI,'simo')
%   [FFqMeas,ADFFq,FFqROI] = mimo(seti,qROI,'simo')
%   [JA,JB] = mimo(seti,qROI,'jacobian')
%
%
%% Description
%
% *Adjoint of derivative ('adjOfDer')*
%
%   [FFqmF, ADFFq] = mimo(seti,qROI,FmeasDelta)
%
% and 
%
%   [FFqmF, ADFFq] = mimo(seti,qROI,FmeasDelta,'adjOfDer')
%
% * computes the difference of the evalutated forward operator and the 
%   measured data, |FFqmF|, and
% * computes the adjoint of derivative, |ADFFq|.
% * You can leave |'adjOfDer'| because it is used automatically if the last 
%   argument |op| is empty.
%
% *Single-input-multiple-output ('simo')*
%
%   [FFqMeas,~,FFqROI] = mimo(seti,qROI,'simo')
%
% * computes the scattered field evaluated on receivers' positions for each transmitter, |FFqMeas|, and
% * the scattered field on ROI for each transmitter, |FFqROI|.
% * ~ is the adjoint of derivative, |ADFFq|, but this is uninteresting because it is 0.
%   (This is expected, because |simo| is used with input of FmeasDelta = 0 in |mimocalc|.)
%
% *Auxilary matrices JA and JB to compute the Jacobian matrix ('jacobian')*
%
%   [JA,JB] = mimo(seti,qROI,'jacobian')
%
% computes auxiliary matrices JA and JB, such that the Frechet derivative 
% of $\mathcal{F}$ at $q$ evaluated on $h$, i.e. $\mathcal{F}'(q)[h]$, that
% is in finite-dimensional spaces is represented by the Jacobian matrix of 
% $\mathcal{F}$ at $q$ evaluted on $h$ can easily computed by
%
% $$\texttt{DFF} = \mathcal{F}'(q)[h] = \texttt{JA*diag(h)*JB}$$
%
%
%% Example
%
% We give one example to compute the difference of the evalutated forward 
% operator and the measured data, |FFqmF|, and the adjoint of derivative, |ADFFq|.
%
%   init;
%   seti.contrast = 'corner2D';
%   seti = setGeomSim(seti);
%   qROI = seti.qROIexact;
%   FmeasDelta = zeros(seti.measNb,seti.incNb);
%   [FFqmF,ADFFq] = mimo(seti,qROI,FmeasDelta,'adjOfDer');
%   figure(1); imagesc(real(FFqmF)); axis xy; colorbar;
%   figure(2); imagesc(real(seti.G(ADFFq))); axis xy; colorbar;
%
%
%% Input Arguments
%
% * |seti|       :    struct, details see below.
% * |qROI|       :    contrast in ROI (region of interest)
%                     (matrix as vector of size seti.nROI^seti.dim x 1)
% * |FmeasDelta| :    $\mathcal{F}_\mathrm{meas}^\delta$ measurement data (on receivers) with noise level $\delta$
%                     (matrix of size seti.measNb x seti.incNb)
%
% *struct |seti|*
%
% As in convenience functions 
% <adjOfDer.html>, <forward.html> and <derivative.html>
% all necessary fields (and more) are set by 
%
%   seti = setGeomSim(seti);
%
% A list of all necessary (and not more) fields is below.
%
% seti in <setGeomSim.html>:
%
% * |seti.mCD|      : internal parameter, set seti.mCD = 0 
%                     (no coarse grid is used)
% * |seti.tol|      : internal parameter, set 1E-6, is used in simo and solveLippmannSchwinger
%
% seti in <setGrid.html>:
%
% * |seti.dim|     : dimension of the problem (2 or 3)
% * |seti.nCD|     : number of discretization points for each dimension
%                    of computational domain (CD)
% * |seti.nROI|    : discretization points for each dimension
%                    of region of interest (ROI)
% * |seti.gridROI| : grid of region of interest (ROI) 
%                    (matrix of size seti.dim x seti.nROI^seti.dim)
% * |seti.ROImask| : biggest square in 2D (cube in 3D) in ball with radius rCD/2
%                    (logical array of size seti.nCD x seti.nCD)
% * |seti.dV|      : pixel/voxel area/volume (number)
% 
% seti in <setKernel.html>:
%
% * |seti.model|   : model (helmholtz2D or helmholtz3D)
% * |seti.k|       : wave number
% * |seti.kHat|    : Fourier coefficients of discretized integral operator
%                    (seti.nCD x seti.nCD)
%
% seti in <expSetup.html> (subfiles <pntsGeometry.html>,
% <setIncField.html>, <setMeasKer.html>):
%
% * |seti.measNb|   : number of receivers (number of measurements)
% * |seti.incNb|    : number of transmitters (number of incident fields)
% * |seti.dSInc|    : approximation of the infinitesimal element of closed contour with control points in case of transmitters
% * |seti.dSMeas|   : approximation of the infinitesimal element of closed contour with control points in case of receivers
% * |seti.incField| : incident field, <setIncField.html>
%                     (matrix of size seti.nROI^seti.dim x seti.incNb)
% * |seti.measKer|  : matrix to compute simulated measurement data from solution, <setMeasKer.html>
%                     (size seti.measNb x seti.nROI^seti.dim)
%
%
%% Output Arguments
%
% More details of output arguments are in section _More About_.
% type: vector, matrix, etc...
%
% * |FFqMeas|  :    scattered field evaluated on the receivers' positions
%                   for each transmitter
%                   (complex matrix of size seti.measNb x seti.incNb)
% * |FFqmF|    :    |FFqMeas| minus noisy data
%                   (complex matrix of size seti.measNb x seti.incNb)
% * |ADFFq|    :    adjoint of derivative
%                   (complex matrix written as vector of size seti.nROI^seti.dim x 1)
% * |FFqROI|   :    scattered field evaluated on ROI (region of interest)
%                   stored as vector for each transmitter
%                   (complex matrix of size seti.nROI^seti.dim x seti.incNb)
% * |JA|, |JB| :    Auxilliary matrices
%                   (JA is a complex matrix of size seti.measNb x seti.nROI^seti.dim,
%                    JB is a complex matrix of size seti.nROI^seti.dim x seti.incNb)
%
%% Best Practice
%
% If you are only interested in the scattered field evaluated on receivers' 
% positions (i.e. measurements), |FFqMeas|, use
%
%   [FFqMeas,~,~] = mimo(seti,qROI,'simo')
%
% instead of
%
%   [FFqMeas,~] = mimo(seti,qROI,0,'adjOfDer')
%
% or
%
%   [FFqMeas,~] = mimo(seti,qROI,zeros(seti.measNb,seti.incNb),'adjOfDer')
%
% The latter with |'adjOfDer'| computes the adjoint of derivative additionally, 
% which is more expensive than |'simo'|.
%
%
%% More About
%
% * $\texttt{FF} = \mathcal{F}$ : contrast-to-measurement operator
% * $\mathcal{F}(q)$ is the matrix of TX-RX scattered fields
%   (TX: transmitters, RX: receivers), i.e. multiple-input-multiple-output
%   setting, i.e. the k-the column is the RX-field for TX(k).
% * $\texttt{FmeasDelta} = F_\mathrm{meas}^\delta$, i.e. 
%   with noise level $\delta$ noisy simulated data (measurement) $F_\mathrm{meas}$.
% * $f_\mathrm{dis} = (1/2)\ \|\mathcal{F}(q) -
% F_\mathrm{meas}^\delta\|_\mathrm{dis}^2$ : discrepancy
% * $\texttt{FFqMeas}$ is the scattered field evaluated on receivers'
% positions (measurements) (also called $\texttt{uScattRX}$).
%   Note that $\texttt{FFqMeas} = \mathcal{F}(q)$.
% * $\texttt{FFqmF} = \texttt{FFqMeas} - \texttt{FmeasDelta} =
% \mathcal{F}(q) - F_\mathrm{meas}^\delta$.
%
% * $\texttt{FFqROI}$ is the scattered field evaluated on ROI 
%   stored for each transmitter for each transmitter 
%   (for one transmitter it is also called $\texttt{uScattROI}$).
% * Scattered field $\texttt{FFqROI} = u_\mathrm{ROI}^s = V(q \cdot
% u_\mathrm{ROI}^t)$ with volume potential $V$.
% * Remind: FFqROI is on ROI, but FFqmF is on RX (receivers positions).
%
% * $\texttt{DFFqh} = \texttt{JA*diag(h)*JB} = \mathcal{F}'(q)[h]$ 
%   is the Frechet derivative 
%   with auxiliary matrices $\texttt{JA}$ and $\texttt{JB}$.
% * $\texttt{ADFFq} = [\mathcal{F}'(q)]^\ast [\mathcal{F}(q) -
% F_\mathrm{meas}^\delta]$
%   is the adjoint of derivative.
%
% * |[JA,JB] = mimo(seti, qROI,'jacobian')| 
%   is used in case of primal-dual algorithm.
% * |[FFqmF,ADFFq] = mimo(seti,qROI,seti.FmeasDelta)| 
%   is used in case of soft-shrinkage (not available in public version).
%
%% See Also
% For easy usage we offer some convenience functions for mimo:
%
% * adjOfDer: <adjOfDer.html>
% * simo:     <forward.html>
% * jacobian: <derivative.html>
%

%% Code: mimo
function varargout = mimo(seti, qROI, varargin)

if ~ischar(varargin{end})
    op = 'adjOfDer';
else
    op = varargin{end};
end

if strcmpi(op, 'adjOfDer') % old name: gradient
    [varargout{1}, varargout{2}] = mimocalc(seti, qROI, varargin{1}, op);
    % Note: varargin{1} = FmeasDelta.
    % Output documentation:
    % varargout{1} = FFqmF
    % varargout{2} = ADFFq
    
elseif strcmpi(op, 'simo')
    [varargout{1}, varargout{2}, varargout{3}] = mimocalc(seti, qROI, 0, op);
    % Output documentation:
    % varargout{1} = uScattRX = FFqMeas
    %     (and not FFqMeas - FmeasDelta, because FmeasDelta was set to 0).
    % varargout{2} = ADFFq = 0 (in this case)
    % varargout{3} = uScattROI = FFqROI
    
elseif strcmpi(op, 'jacobian')
    if strcmp(seti.model, 'helmholtz3D') ...
            || strcmp(seti.model, 'helmholtz2D') ...
            || strcmp(seti.model, 'helmholtzHMode2D')
        [~, ~, ~,varargout{1},varargout{2}] = mimocalc(seti, qROI, 0, op); 
        % Output documentation:
        % varargout{1} = JA
        % varargout{2} = JB
    else
        disp('mimo.m: Error - Jacobian is not (yet ... ?) implemented');
        varargout{1} = [];
        varargout{2} = [];
    end
end

end

%% Code: mimocalc
function [FFqmF, ADFFq, FFqROI, JA, JB] = mimocalc(seti, qROI, FmeasDelta, op)
%Note: HMode is not available in public version.

nnROI = length(seti.gridROI(1,:));
ADFFq = 0;
FFqROI = zeros(nnROI, seti.incNb);
FFqmF = zeros(seti.measNb, seti.incNb);

if strcmpi(op, 'adjOfDer')
    ADFFq = zeros(size(qROI));
end
if strcmpi(op, 'jacobian') 
    if strcmp(seti.model, 'helmholtz3D') || strcmp(seti.model, 'helmholtz2D')
        JA = zeros(seti.measNb, nnROI); JB = zeros(nnROI, seti.incNb);  
    elseif strcmp(seti.model, 'helmholtzHMode2D')
        JA = zeros(seti.measNb, nnROI,2); JB = zeros(nnROI, seti.incNb,2);
    else
        disp('mimo.m: Error - Jacobian is not (yet ... ?) implemented');
    end
end

dSInc = seti.dSInc.*ones(seti.incNb,1);

[~, S, VGStar, VStar, TStar, ~, QStarU, ~, ~] = intOpsFuncs(qROI,seti);

for jj = 1:seti.incNb
    SL = dSInc(jj)*seti.incField(:,jj);
    % SL = uIncROI (i.e. single layer potential/Herglotz wave function)
    % SL \phi = uIncROI
    [uMeas, uScattROI] = S(SL); % uMeas is uScattRX (i.e. u^s on measPnts)
        
    h = uMeas;
    
    if strcmpi(op, 'adjOfDer')
        h = h - FmeasDelta(:,jj); % result: h = FFqmF = FF(q) - FmeasDelta
    end;
    FFqmF(:,jj) = h; % In case of adjOfDer (default): FFqmF = FF(q) - FmeasDelta
                     % otherwise so just FFq:         FFqmF = FF(q)
    h( isnan(h) ) = 0; % workaround
    uTotROI = uScattROI + SL; % uTotROI = uScattROI + uIncROI
    if strcmpi(op, 'adjOfDer')
        if strcmp(seti.model, 'helmholtz3D') || strcmp(seti.model, 'helmholtz2D')
            f = VGStar(h);
            v = conj(uTotROI).*(f + VStar(TStar(QStarU(f)))); %VStar NOT VqStar
            ADFFq = ADFFq + v;
            %figure(101); x = ADFFq; imagesc(real(extendROItoCD(x,seti.ROImask))); axis xy; colorbar; caxis([-3,3]); error('stop');
        elseif strcmp(seti.model, 'helmholtzHMode2D')
            ADFFq = ADFFhelmHMode2D(h,seti,uTotROI);
        else
            disp('mimo.m: Error - model is not (yet ... ?) implemented');
        end
    end
    if strcmpi(op, 'simo')
        FFqROI(:,jj) = reshape(uScattROI,[nnROI 1]);
        ADFFq = 0;
    end
    if strcmpi(op, 'jacobian')
        if strcmp(seti.model, 'helmholtz3D') || strcmp(seti.model, 'helmholtz2D')
            JB(:,jj) = uTotROI;
        elseif strcmp(seti.model, 'helmholtzHMode2D')
            [uTotROI1, uTotROI2] = gradientHMode(uTotROI,seti);
            JB(:,jj,1) = uTotROI1;
            JB(:,jj,2) = uTotROI2;
        else
            disp('mimo.m: Error - Jacobian is not (yet ... ?) implemented');
        end 
    end
end
if strcmpi(op, 'jacobian')
    dSMeas = seti.dSMeas.*ones(seti.measNb,1);
    if strcmp(seti.model, 'helmholtz3D') || strcmp(seti.model, 'helmholtz2D')
        for jj = 1:seti.measNb
            f = VGStar( ((1:seti.measNb)==jj)' * seti.dV/dSMeas(jj) );
            JA(jj,:) = conj(f + VStar(TStar(QStarU(f)))); %VStar NOT VqStar
        end;
    elseif strcmp(seti.model, 'helmholtzHMode2D')
        for jj = 1:seti.measNb 
            [f1,f2] = helmholtzHMode2DZStar( ((1:seti.measNb)==jj)' * seti.dV/dSMeas(jj), seti); % Z*(g_j) 
            qf1 = conj(seti.qROIexact).*f1; qf2 = conj(seti.qROIexact).*f2;
            
%             qf1Ext = 1.0i*pi/seti.rCD * fft2(extendROItoCD(qf1,seti.ROImask)) * seti.gradMat;
%             qf2Ext = 1.0i*pi/seti.rCD * seti.gradMat * fft2(extendROItoCD(qf2,seti.ROImask));
%             qf1 = restrictCDtoROI(ifft2(qf1Ext+qf2Ext),seti.ROImask);

            qf1 = divHMode(seti,qf1,qf2);

            qf1 = TStar(qf1);
        
            [qf1, qf2] = helmholtzHMode2DNablaVStar(qf1,seti);
            
            JA(jj,:,1) = conj(f1 + qf1); JA(jj,:,2) = conj(f2 + qf2); 

        end;    
    else
        disp('mimo.m: Error - Jacobian is not (yet ... ?) implemented');
    end 
end

end

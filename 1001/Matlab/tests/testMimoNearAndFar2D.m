%% testMimoNearAndFar2D
% Test near field data against far field data.
%
%% Syntax
%
%   testMimoNearAndFar2D
%
%% Description
%
% |testMimoNearAndFar| compares by <mimo.html> computed scattered fields
% for near field and far field (at receivers positions (RX) and on region of
% interest (ROI)) in 2D case.
%
% * testset 1: compares <mimo.html> against reference solution 
%              (as in <testMimo.html>) and further compares 
%              near field (with huge radius radSrc) against far field results.
%              For comparison the relative error is computed.
% * testset 2: compares by <mimo.html> computed near with far field solution 
%              in case of Fresnel op. 1 setting (only the setting, not the
%              data is required).
%              (More about Fresnel settings, see <readRAWData.html> and
%              <loadData.html>).
%              (There is no comparison with reference solution).
%              For comparison the relative error is computed.
%
% This test should run in environment runtests, see <runtests.html> via
%
%   runtests('mnf');
%
%% Input Arguments
%
% See <testMimo.html>
%
%% Output Arguments
%
% Terminal output.
%
%% More About
%
% To compare the far field data with near field, we 
% compute the field $u(x)$ from far field $u_\infty(x)$.
%
% The relation is, see [1, Th. 2.6]: 
%
% $u(x) = \exp(\mathrm{i} k |x|)/(|x|^{(d-1)/2}) 
%  (u_\infty(\hat{x} + \mathcal{O}(|x|^{-1})
%  \quad \mathrm{if} \quad |x| \rightarrow \infty$.
%
%% References
%
% * [1] David Colton and Rainer Kress. _Inverse Acoustic and Electromagnetic Scattering Theory_. Springer, New York, 2013.
%
%% See Also
%
% * <runtests.html>
% * <testMimo.html>
% * <referenceData2D.html>
% * <referenceData3D.html>
% * <readRAWData.html>
% * <loadData.html>
%
%% Code
%
% *Only in case of standalone (not recommended)*

% init;

%%
% *testsets*

for testset = 1:2
    
    disp(' ');
    fprintf('-- testset %d --\n',testset);
    disp(' ');

    seti.dim = 2;
    seti.model = 'helmholtz2D';

    if testset == 1 % planeWave and incNb = 1 to compare with reference solution
        seti.nCD = 2^9;
        seti.contrast = 'referenceBall2D'; % contrast function (shape)
        seti.incType = 'planeWave'; %use planeWave(!) (referenceData2D.m has only planeWave as reference data)
        seti.incNb = 1; %1! otherwise test will not run of course (referenceData is plane wave incidence with direction (1,0))
        seti.k = 3; % wave number
        seti.rCD = 2.5;
        seti.qBall = 0.5; % contrast of ball
        seti.rBall = 0.5; % radius of ball
        seti.measNb = 2;
        seti.radMeas = 200;
        % radMeas must be high (otherwise you see the difference between
        % usFarNear bzw. usFarRefNear againt usNearRef)
        seti.radSrc = 4;
    elseif testset == 2 % Fresnel
        seti.nCD = 2^8;
        seti.contrast = 'fresnel_op1_twodielTM';
        seti.incType = 'pointSource'; % can not be compared with reference solution, because there ist not reference solution for point sources.
        seti.incNb = 36;
        f = 4E9; % frecuency 4 GHz
        % if f is smaller, the difference between near and far field on RX
        % descreases...
        % e.g. 4 GHz:
        %usNearROI against usFarROI   :   0.000000
        %usNearRX  against usFarRXNear:   0.040871
        % e.g. 1 GHz:
        %usNearROI against usFarROI   :   0.000000
        %usNearRX  against usFarRXNear:   0.014894
        % Result:
        % to divide Fresnel data with factor (see below) and use them as
        % far field data does not seem to be a good approximation
        seti.k = 2*pi/(3E8/f);
        seti.rCD = 0.2;
        seti.measNb = 72;
        seti.radSrc = 0.72;
        seti.radMeas = 0.76;
        %seti.radMeas = 200; % to test... and yes, does fit better to farField...
    end

    for i = 1:2

        if i == 1
            seti.measType = 'nearField';
        elseif i == 2
            seti.measType = 'farField';
        end

        seti = setGeomSim(seti);

        N = seti.nROI;

        if strcmp(seti.model, 'helmholtz2D') && testset == 1
            % && strcmp(seti.incType,'planeWave')
            % && strcmp(seti.incNb,1) (to be exactly; it is testset = 1)
            % but planeWave indicates that we want to compare with reference
            % solution
            disp('Compute reference solution')
            tic
            [uScattRXRef, uScattROIRef] = referenceData2D(seti);
            toc
        elseif strcmp(seti.model, 'helmholtz2D') && testset == 2
            disp('No reference solution available (so no comparison with ref. sol.).')
        else
            disp('testMimo: Unknown Model!');
            return;
        end

        disp('Compute solution')
        tic
        [uScattRX,~,uScattROI] = mimo(seti, seti.qROIexact, 'simo');
        toc

        if i == 1
            usNearROI = uScattROI;
            usNearRX = uScattRX;
            if testset == 1
                usNearROIRef = uScattROIRef;
                usNearRXRef = uScattRXRef;
            end
            measPntsNear = seti.measPnts;
        elseif i == 2
            usFarROI = uScattROI;
            usFarRX = uScattRX;
            if testset == 1
                usFarROIRef = uScattROIRef;
                usFarRXRef = uScattRXRef;
            end

            % usFarRXNear: Compute field u(x) from farField u_\infty(x)
            % to compare the far field data with near field
            % u(x) = exp(i k |x|)/(|x|^{(d-1)/2}) (u_\infty(\hat{x} + \mathcal{O}(|x|^{-1})
            % if |x| \rightarrow \infty
            usFarRXNear = zeros(size(usFarRX));
            usFarRXRefNear = zeros(size(usFarRX));
            RX = measPntsNear; % use measPnts from nearField, because farField is on circle with radius 1...
            k = seti.k;
            for j = 1:seti.measNb
                absx = sqrt(RX(1,j)^2+RX(2,j)^2);
                factor = exp(1i*k*absx)/(absx)^((seti.dim-1)/2);
                usFarRXNear(j,:) = factor*usFarRX(j,:);
                if testset == 1
                    usFarRXRefNear(j,:) = factor*usFarRXRef(j,:);
                elseif testset == 2
                    usFarRXNear(j,:) = factor*usFarRX(j,:);
                end
            end

        end

    %    errorROI = norm( seti.ballMask(seti.ROImask).*(uScattROI-uScattROIRef) )...
    %        /norm( seti.ballMask(seti.ROImask).*uScattROIRef );
    %    errorCD = norm( uScattRX-uScattRXRef )/norm( uScattRXRef );
    %    fprintf('%4d^%i  |  %10f |  %10f |  %10f | %f\n',NMat, seti.dim, errorROI,errorCD,tgrid,t);

    end

    if testset == 1
        disp('testset = 1: planeWave, incNb = 1')
        disp('compare scatterd near field with reference solution')
        disp('compare scatterd far field with reference solution')
        disp('compare scatterd near and far field (near field with huge radius)')

        errorROI = @(uScattROI,uScattROIRef)...
                   norm( seti.ballMask(seti.ROImask).*(uScattROI-uScattROIRef) )...
                   /norm( seti.ballMask(seti.ROImask).*uScattROIRef );
        errorCD = @(uScattRX,uScattRXRef) norm( uScattRX-uScattRXRef )/norm( uScattRXRef );

        fprintf('ROI: usNear       against usNearRef: %10f\n',errorROI(usNearROI,usNearROIRef)')
        fprintf('RX : usNear       against usNearRef: %10f\n',errorCD(usNearRX,usNearRXRef))

        fprintf('ROI: usFar        against usFarRef : %10f\n',errorROI(usFarROI,usFarROIRef))
        fprintf('RX : usFar        against usFarRef : %10f\n',errorCD(usFarRX,usFarRXRef))

        fprintf('ROI: usFar        against usNearRef: %10f\n',errorROI(usFarROI,usNearROIRef))
        fprintf('RX : usFarNear    against usNearRef: %10f\n',errorCD(usFarRXNear,usNearRXRef))

        fprintf('ROI: usFarRef     against usNearRef: %10f\n',errorROI(usFarROIRef,usNearROIRef))
        fprintf('RX : usFarRefNear against usNearRef: %10f\n',errorCD(usFarRXRefNear,usNearRXRef))

        %figure(101);plot(real(uScattRXRef))
        %figure(102);plot(real(uScattRX))
        % test factor error: uScattRX./uScattRef
    elseif testset == 2
        disp('testset = 2: pointSource; Fresnel, op. 1')
        disp('compare scatterd near and far field')
        errorROI = @(uScattNearROI,uScattFarROI)...
                   norm( seti.ballMask(seti.ROImask).*(uScattNearROI-uScattFarROI) )...
                   /norm( seti.ballMask(seti.ROImask).*uScattFarROI );
        errorCD = @(uScattNearRX,uScattFarRX) norm( uScattNearRX-uScattFarRX )/norm( uScattFarRX );

        fprintf('usNearROI against usFarROI   : %10f\n',errorROI(usNearROI,usFarROI)')
        fprintf('usNearRX  against usFarRXNear: %10f\n',errorCD(usNearRX,usFarRXNear))

    end
    
end

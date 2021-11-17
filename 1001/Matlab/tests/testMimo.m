%% testMimo
% Test of function mimo (multiple input, multiple output).
%
%% Syntax
%
%   testMimo
%
%% Description
%
% |testMimo| compares by <mimo.html> computed scattered fields
% (at receivers positions and on region of interest )
% with reference solutions
% (<referenceData2D.html> and <referenceData3D.html>).
% This comparison is done for several grid sizes.
%
% This test should run in environment runtests, see <runtests.html> via
%
%   runtests('m')
%
%% Input Arguments
%
% * seti.dim    : dimension of the problem (2 or 3)
% * seti.model  : model 'helmholtz2D' or 'helmholtz3D'
%                 (other options are not supported in the public version.)
%
% Further parameters are desribed in the code in the 
% section "General settings", "Discretization parameters", and
% "Experimental set-up".
%
%% Output Arguments
%
% *Figures*, see <runtests.html>.
%
%% More About
%
% * If the error is high, set a higher "nMax = 20" in referenceData*D.m.
%
%% See Also
%
% * <runtests.html>
% * <referenceData2D.html>
% * <referenceData3D.html>
%
%% Code
%
%%
% *Run this file standalone... (not recommended)*
%
% If you call testMimo 'standalone' you can use the following inside if

if 0
    init;
    seti = struct; % old: setInput;
    seti.dim = 2;
    seti.model = 'helmholtz2D';
end

%%
% *General settings*
%
% * seti.measType   :   type of measurement: 'farField' or 'nearField',
%   see <expSetup.html>.
% * 'farField' is available in 2D and 3D.
% * 'nearField' is only available in 2D (because 3D reference data are not implemented).
%   See <referenceData2D.html>, <referenceData3D.html>.
%
% * useReferenceData3D  :   Use reference data in case of 3D (0 or 1) 
%                           (It is expensive to us it).
%                           Note: If reference data is not used rel. errors are 0.
%

seti.measType = 'farField'; % (can be choosen in 2D and 3D)
%seti.measType = 'nearField'; % (can be choosen only in 2D!)

useReferenceData3D = 1;

HMode = 0;      % if 1: using helmholtzHMode instead of helmholtz (not available in public version).
useCoarse = 0;  % if 1: using coarse grid yes or no (logical) (not available in public version and not recommended).

%%
% *Setting: check consistency*

if ~exist('seti', 'var') || ~isfield(seti,'dim')
    seti.dim = 2;
    if HMode
        seti.model = ['helmholtzHMode' int2str(seti.dim) 'D'];
    else
        seti.model = ['helmholtz' int2str(seti.dim) 'D'];
    end
end
seti.fileSuffix = sprintf('_testMimo%dD',seti.dim);

%% 
% *Discretization parameters*
%
% * NMat    :   vector with entries for different seti.nCD 
%               (number of discretization points for each dimension of CD).

% Note: Mexp = -1 -> M = floor(1/seti.nCD) = 0 -> twogrid is not used
if seti.dim == 3
    longtest = 1;
    if longtest
        NMat = sort([2.^(6:8) 2.^(6:8)-1]);
    else
        NMat = sort([2.^(6:6) 2.^(6:6)-1]);
    end
elseif seti.dim == 2
    if HMode % HMode not supported in public version
        NMat = sort([2.^(7:10) 2.^(7:10)-1]);
    else
        NMat = sort([2.^(5:9) 2.^(5:9)-1]);
    end
end

if useCoarse % useCoarse is not supported in public version.
    Mexps = [-1, 5/6, 2/3, 1/2];
else
    Mexps = 0; % Something with length 1.
    % mCD = 0 will be set later, so the value Mexps will not be used
    % but the length is used in a for-loop.
end

%%
% *Experimental set-up*
%
% _Changeable parameters_
%
% * seti.k          :   wave number
%
% * seti.qBall      :   contrast value of the ball
% * seti.rBall      :   radius of the ball
% * seti.rCD        :   size of computational domain [-rCD,rCD)^dim
%
% * seti.measNb         :   number of receivers, default: 35
% * seti.radSrc         :   radius of circle for transmitters, default: 5
%                           (not necessary in case of incType='planeWave').
% * seti.radMeas        :   radius of circle for receivers, default: 5
%                           (not necessary in case of measType='farField')
%
% _Fixed parameters (for this test)_
%
% * seti.contrast   :   'referenceBall2D' or 'referenceBall3D'
%                       (currently no reference data is available for other
%                       contrasts.)
% * seti.incType    :   type of incident wave; in this case 'planeWave' is
%                       required, because reference data referenceData2D.m 
%                       and referenceData3D.m requires a plane wave.
% * seti.incNb      :   number of transmitters; in this case 1 is required
%                       because referenceData requires an incident plane 
%                       wave with direction (1,0) in 2D or (1,0,0) in 3D.
 
seti.k = 2;

% helmholtz*D
seti.qBall = 0.8;
seti.rBall = 0.015;
seti.rCD = 0.2;

% helmholtzHMode*D (not available in public version)
%seti.qBall = 3.0; % contrast of ball
%seti.rBall = 0.5; % radius of ball
%seti.rCD = 2.5;

seti.measNb = 2;
seti.radSrc = 4; 
seti.radMeas = 4; 

% Fixed parameters:
seti.contrast = ['referenceBall' int2str(seti.dim) 'D'];
seti.incType = 'planeWave';
seti.incNb = 1;

%%
% *setGeomSim*

disp('# setGeomSim')
seti = setGeomSim(seti); % set geometry and simulation
%seti = setGeomSim(seti,1); % to debug with output...

%%
% *testing mimo*

fprintf('testing mimo for dim = %i, model = %s\n', seti.dim, seti.model);
fprintf('\n        |         | rel. error  | rel. error  | grid. gen.  |        ');
fprintf('\nPoints  | Points  | of u_s      | of u_s      | +ref.u gen. | comp.  ');
fprintf('\nin CD(N)| in CD(M)| on ROI      | at RX       | time        | time.  ');
fprintf('\n------- + ------- + ----------- + ----------- + ----------- + --------\n');

for j = 1:length(NMat)
    if j > 1 && length(Mexps) > 1 fprintf('------- + ------- + ----------- + ----------- + ----------- + --------\n');end%#ok<SEPEX>
    
    tic;
    seti.nCD = NMat(j); % number of collocation points in 1 dimension of CD (computational domain)

    seti = setGeomSim(seti);

    if strcmp(seti.model, 'helmholtz2D')
        [uScattRXRef, uScattROIRef] = referenceData2D(seti);
    elseif strcmp(seti.model, 'helmholtzHMode2D')
        [uScattRXRef, uScattROIRef] = referenceDataHMode2D(seti);
    elseif strcmp(seti.model, 'helmholtz3D')
        % If you want to skip the reference solution
        % and compare with the first computation for this N, 
        % then change useReferenceData3D to false (0).
        % Note that the reference solution is expensive.
        if useReferenceData3D
            [uScattRXRef, uScattROIRef] = referenceData3D(seti);
        else
            uScattRXRef = 0;
            uScattROIRef = 0;
        end
    else
        disp('testMimo: Unknown Model!');
        return;
    end
    tref = toc;

    for k = 1:length(Mexps)
        tic;
        
        if useCoarse
            seti.mCD = floor(seti.nCD^Mexps(k));
            [~,seti] = evalc('setCoarse(seti);');
        else
            seti.mCD = 0;
        end
        
        E = @(x) x;
        N = seti.nROI;
        %E = @(x) extendROItoCD(x,seti.ROImask);
        %N = seti.nCD;
        
        tgrid = toc;
        tic;
        
        [FFqMeas,~,FFqROI] = mimo(seti, seti.qROIexact, 'simo');
        uScattRX = FFqMeas;
        if seti.incNb == 1
            uScattROI = FFqROI(:,1);
        else
            error('testMimo does not support more than one transmitter yet.')
        end
        
        t = toc;
        
        % If reference solution is not set use this one
        if uScattROIRef == 0
            uScattRXRef = uScattRX;
            uScattROIRef = uScattROI;
        end
        
        errorROI = norm( seti.ballMask(seti.ROImask).*(uScattROI-uScattROIRef) )...
            /norm( seti.ballMask(seti.ROImask).*uScattROIRef );
        errorRX = norm( uScattRX-uScattRXRef)/norm(uScattRXRef);
        
        fprintf(' %4d^%i |  %4d^%i |  %10f |  %10f |  %10f | %f\n',seti.nCD, seti.dim, seti.mCD, seti.dim, errorROI, errorRX, tgrid+tref, t);
               
        part = @real;
        s = N*ones(1,seti.dim);
        
        EuScattROI = reshape(part(E(uScattROI)), s);
        EuScattROIRef = reshape(part(E(uScattROIRef)), s);
        Ediff = reshape(part(E(uScattROI-uScattROIRef)), s);
        
        switch seti.dim
            case 2
                sliceWarn = '';
            case 3
                % Multiple slices
                v = round(linspace(1, N, min(10,N)));
                
                figure(28); slice(EuScattROI,v,N,1); title('Computed scattered field'); colormap(litman); colorbar; shading flat;
                figure(29); slice(EuScattROIRef,v,N,1); title('Reference scattered field'); colormap(litman); colorbar; shading flat;
                figure(30); slice(Ediff,v,N,1); title('Difference'); colormap(litman); colorbar; shading flat;
                
                % Display a slice of the data.
                sliceWarn = ' (slice z = 0)';
                N2 = ceil((N-1)/2);
                EuScattROI = EuScattROI(:,:,N2);
                EuScattROIRef = EuScattROIRef(:,:,N2);
                Ediff = Ediff(:,:,N2);
        end
        
        if seti.dim == 2
            figNo1 = 21; figNo2 = 22; figNo3 = 23;
        elseif seti.dim == 3
            figNo1 = 25; figNo2 = 26; figNo3 = 27;
        end
        figure(figNo1); imagesc(EuScattROI); axis xy; title(['Computed scattered field' sliceWarn]);colormap(litman);colorbar;
        figure(figNo2); imagesc(EuScattROIRef); axis xy; title(['Reference scattered field' sliceWarn]);colormap(litman);colorbar;
        figure(figNo3); imagesc(Ediff); axis xy; title(['Difference' sliceWarn]);colormap(litman);colorbar;
        drawnow;
        savePngFig(figNo1,0,seti);
        savePngFig(figNo2,0,seti);
        savePngFig(figNo3,0,seti);


    end
end
fprintf('\n');

%% setContrast
% Set contrast function and evaluates it on CD and ROI.
%
% Warning: This function uses feval. Be careful when using.
%
%% Syntax
%
%   seti = setContrast(seti)
%   seti = setContrast(seti,dispDepth)
%   seti = setContrast(seti,dispDepth,out)
%
%% Description
% |seti = setContrast(seti)| sets the contrast defined in |seti.contrast|
% on the computational domain (CD) and the region of interest (ROI) and
% stores them as vectors in |seti.qCDexact| and |seti.qROIexact|.
%
% |seti = setContrast(seti,dispDepth)| does the same but allows to control
% displayed messages by |dispDepth|.

% |seti = setiContrast(seti,dispDepth,out)| does the same but deals with the output 
% depth |out|: no figure (0), plot figure (1), plot and save figure (2).
%
% * If |seti.rotation| contains a value, the contrast is rotated around
%   this number in degrees in mathematical positive sense.
%
% *Structure of the program*
%
% * Checks consistency of the field contrasts, i.e. essentially sets a 
%   valid seti.contrast (name of file in folder incontrasts).
%
% - If no real-world (experimentally measured) data is used (|expData| is empty) and no
% contrast is defined, set default contrasts in dependence of dimension
% (|seti.dim|).
%
% - If no real-world data is used and the contrast is
% |referenceBall|, choose automatically referenceBall2D or referenceBall3D 
% in dependence of the choosen dimension (seti.dim).
%
% * Deals with real-world data (if seti.expData = 'fresnel')
%
% - If |seti.fresnelFile| contains "dielTM" 
%   choose |seti.contrast = 'fresnel_op1_dielTM'|.
%
% - If |seti.fresnelFile| contains "twodielTM"
%   choose |seti.contrast = 'fresnel_op1_twodielTM'|.
%
% - If |seti.fresnelFile| contains "rectTM_cent"
%   choose |seti.contrast = 'fresnel_op1_rectTM_cent'|.
%
% * Rotate grid to rotate the contrast.
%
% * Evalutes the predefined contrast in CD and ROI.
%
% * Restricts the contrast:
%   Contrast has to be restricted to B(0,rCD/2) which is seti.ballMask
%   (seti.ROImask is rectangular with length rCD/2).
%
% * Saves the contrast as a vector (in ROI and CD).
%
% * Plots the predefined contrast, see <plotPredefinedContrast.html>, if
% |out| was set to 1 or 2. In case of |out| = 2 the plot is additionally
% saved.
%
%% Example
%
%   init;
%   seti.dim = 2;                           % dimension of the problem
%   seti = setGrid(seti);                   % set grids CD and ROI
%   seti.contrast = 'cornerBallSparse2D';   % predefined contrast
%   seti.rotation = 20;                     % rotate the contrast 20 degrees
%
%   seti.plotVisible = 'on';                % required in plotPredefinedContrast
%   seti.usecbarlim = 0;                    % required in plotPredefinedContrast
%   seti.plotPublish = 0;                   % required in plotPredefinedContrast
%
%   seti = setReshapeVecMat(seti);          % definition of seti.G is required
%   seti = setContrast(seti,1,1);
%
%% Input Arguments
%
% * seti    :   structure array
% 
% * |seti.contrast|   :   name of predefined contrast function
% * |seti.rotation|   :   mathematical positive rotation of the contrast in degrees (only in 2D).
%                         (If the field rotation does not exist, no rotation
%                         is done.)
%
% *Optional Input Arguments*
%
% * |dispDepth|   : Depth of displayed messages (0: no, 1 or greater: yes).
% * |out|         : output depth: no figure (0), plot figure (1), 
%                   plot and save figure (2).
%
%% Output Arguments
%
% * |seti.qCDexact|     :   predefined contrast evalutated on CD
%                           (computational domain) 
%                           (as vector of size seti.nCD^seti.dim x 1).
% * |seti.qROIexact|    :   predefined contrast evalutated on ROI
%                           (region of interest) 
%                           (as vector of size seti.nROI^seti.dim x 1).
%
%% More About
%
% * Contrasts are defined in the folder "incontrasts". To use the 
%   corresponding contrast it is sufficient to write the name of the 
%   function in |seti.contrast| (a file path is not necessary).
% * A list of all predefined contrasts is available in
% <incontrastsRef.html>.
% * IPscatt only uses the contrast evaluated on ROI. 
% If you use |seti.qCDexact| and your own predefined contrasts, note the following:
% The predefined contrasts are defined on ROI via the scaling factor |R = seti.rCD/2|.
% (To be exactly, the mathematical sensible region is the open circle with radius
% |seti.rCD/2| and for sake of simplicity the region of interest is the 
% biggest square inside this mathematical sensible region.) This function
% provides the contrasts' evaluation on CD too. If you want to use it and
% you use your own predefined contrast, make sure that your contrast
% vanishes outside the open ball of radius rCD/2.
%
%% See Also
% * <setGrid.html>
% * <plotPredefinedContrast.html>
% * <incontrastsRef.html>
%
%% Code
%
function seti = setContrast(seti,varargin)

if nargin == 2
    dispDepth = varargin{1};
    out = 0;
elseif nargin == 3
    dispDepth = varargin{1};
    out = varargin{2};
else
    dispDepth = 0;
    out = 0;
end

% Check consistency
seti = contrastConsis(seti,dispDepth);

% Deal with real-world data: choose the corresponding contrast
if isfield(seti,'expData')
    seti = dealWithExpData(seti,dispDepth);
end

% Rotate grid to rotate the contrast
[gridRot, gridRotROI] = rotateGrid(seti);
seti = evalContrast(gridRot,gridRotROI,seti);

% Contrast has to be restricted to B(0,rCD/2) which is seti.ballMask
% (seti.ROImask is rectangular with length rCD/2)
seti = restrictContrast(seti);

% Save contrast as a vector
seti.qCDexact = seti.qCDexact(:);
seti.qROIexact = seti.qROIexact(:);

% Plot predefined contrast
if out >= 1
    plotPredefinedContrast(seti,out);
end

end

%% Code: Subfunctions

%%
% *Code: subfunction: contrastConsis*
%
% Check consistency of seti.contrast.

function seti = contrastConsis(seti,dispDepth)
if ~isfield(seti,'expData') && ~isfield(seti,'contrast')
    if seti.dim == 2
        seti.contrast = 'cornerBallSparse2D';
    elseif seti.dim == 3
        seti.contrast = 'cross3D';
    end
    setmessage(seti,'contrast',dispDepth);
end

if ~isfield(seti,'expData') && (strcmp(seti.contrast,'referenceBall'))
    if seti.dim == 2
        seti.contrast = 'referenceBall2D';
    else
        seti.contrast = 'referenceBall3D';
    end
end

end

%%
% *Code: subfunction: dealWithExpData*
%
% Deals with real-world data: sets the corresponding predefined
% contrast in case of Insitute's Fresnel data if predefined contrast is available.
%
function seti = dealWithExpData(seti,dispDepth)
if strcmp(seti.expData,'fresnel')
    %%
    % *---- expData: Fresnel*
    
    filename = seti.fresnelFile;
    if regexpi(filename,'/dielTM') > 0
        seti.contrast = 'fresnel_op1_dielTM';
        if dispDepth >= 1
            fprintf('   Using %s as exact contrast.\n',seti.contrast);
        end
    elseif regexpi(filename,'/twodielTM') > 0
        seti.contrast = 'fresnel_op1_twodielTM';
        if dispDepth >= 1
            fprintf('   Using %s as exact contrast.\n',seti.contrast);
        end
    elseif regexpi(filename,'/rectTM_cent') > 0
        seti.contrast = 'fresnel_op1_rectTM_cent';
        if dispDepth >= 1
            fprintf('   Using %s as exact contrast.\n',seti.contrast);
        end
    else
        disp('No exact contrast implemented for this data yet.')
        %error('op. 1: rectTM_cent, rectTM_dece, recttTR and uTM not implemented yet; op. 2 too...');
        disp('setting contrast to zero because seti.contrast == "fresnel"');
        seti.qCDexact = zeros(seti.nCD^seti.dim,1); % no exact contrast available on CD
        seti.qROIexact = zeros(seti.nROI^seti.dim,1); % no exact contrast available on ROI
        return;
    end
    %% 
    % *---- expData: Simonetti (out of order...)*
    
%     elseif strcmp(seti.expData,'simonetti')
%         if isfield(seti,'contrast')
%             seti.qCDexact = feval(seti.contrast, seti.grid(1,:), seti.grid(2,:), seti);
%             seti.qROIexact = feval(seti.contrast, seti.gridROI(1,:), seti.gridROI(2,:), seti);
%             seti.qCDexact = seti.qCDexact(:);
%             seti.qROIexact = seti.qROIexact(:);
%         end

%    else
%        error('seti.expData is not implemented in setContrast');
end

end

%%
% *Code: subfunction: rotateGrid*
%
% Rotation of grid to rotate the contrast (in 2D).
%
% The contrast is usually evalutated on seti.grid (CD) respectively
% seti.gridROI (ROI). 
% But we will evaluate it on rotated grid to rotate the contrast.
%
% Feature in future: rotation in 3D case. Make sure to adapt the slice
% through the object then too.

function [gridRot, gridRotROI] = rotateGrid(seti)

if seti.dim == 2 && isfield(seti,'rotation') && seti.rotation ~= 0
  
    % Rotation matrix in R^2
    R = @(alpha) [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)];

    % Compute alpha in rad and set minus
    % such that seti.rotation is mathematical positive rotation in degrees.
    alpha = -seti.rotation/360*2*pi;

    gridRot = R(alpha)*seti.grid;
    gridRotROI = R(alpha)*seti.gridROI;
else
    gridRot = seti.grid;
    gridRotROI = seti.gridROI;
end

end

%%
% *Code: subfunction: evalContrast*
%
% Evaluate the contrast in grid respectively the rotated grid.

function seti = evalContrast(gridRot,gridRotROI,seti)
% Make sure that a file exists and that the path to the file contains the string /incontrasts/.
% For Windows we look for the string \incontrasts\.
if exist(seti.contrast,'file') == 2 && ...
        ( ~isempty(regexp(which(seti.contrast),'/incontrasts/','once')) || ~isempty(regexp(which(seti.contrast),'\incontrasts\','once')) )
    switch seti.dim
        case 2
            seti.qCDexact = feval(seti.contrast, gridRot(1,:), gridRot(2,:), seti);
            seti.qROIexact = feval(seti.contrast, gridRotROI(1,:), gridRotROI(2,:), seti);
        case 3
            seti.qCDexact = feval(seti.contrast, gridRot(1,:), gridRot(2,:), gridRot(3,:), seti);
            seti.qROIexact = feval(seti.contrast, gridRotROI(1,:), gridRotROI(2,:), gridRotROI(3,:), seti);
    end
else
    fprintf('Error: %s is not a file in the folder /incontrasts/ or its subfolders.\n',seti.contrast);
end
end

%%
% *Code: subfunction: restrictContrast*
%
% Restrict contrast
%
% * on ballMask in case of CD and 
% * on ballMaskROI in case of ROI.
%
% Contrast has to be restricted to B(0,rCD/2) which is seti.ballMask
% (seti.ROImask is rectangular with length rCD/2).
%
% Note that the anisotropic case is not supported in the public version.
%
function seti = restrictContrast(seti)

% if... else...: distinguish anisotropic and not anisotropic contrast
if numel(seti.qCDexact) == seti.nCD^seti.dim
    %%
    % No anisotropic contrast
    if norm(seti.qCDexact(:)-seti.qCDexact(:).*seti.ballMask(:))>1e-10
        seti.qCDexact = seti.qCDexact(:).*seti.ballMask(:);
        seti.qROIexact = seti.qROIexact(:).*seti.ballMaskROI(:);
        disp('Warning - contrast has been restricted to B(0,rCD/2)');
    end
else
    %%
    % Anisotropic contrast (not supported in public version)
    disp('Note: Detected anisotropic contrast.');
    
    if norm(seti.qCDexact(:,1)-seti.qCDexact(:,1).*seti.ballMask(:))>1e-10
        seti.qCDexact(:,1) = seti.qCDexact(:,1).*seti.ballMask(:);
        seti.qROIexact(:,1) = seti.qROIexact(:,1).*seti.ballMaskROI(:);
        disp('Warning - contrast has been restricted to B(0,rCD/2)');
    end
    
    if norm(seti.qCDexact(:,2)-seti.qCDexact(:,2).*seti.ballMask(:))>1e-10
        seti.qCDexact(:,2) = seti.qCDexact(:,2).*seti.ballMask(:);
        seti.qROIexact(:,2) = seti.qROIexact(:,2).*seti.ballMaskROI(:);
        disp('Warning - contrast has been restricted to B(0,rCD/2)');
    end
    
    if norm(seti.qCDexact(:,3)-seti.qCDexact(:,3).*seti.ballMask(:))>1e-10
        seti.qCDexact(:,3) = seti.qCDexact(:,3).*seti.ballMask(:);
        seti.qROIexact(:,3) = seti.qROIexact(:,3).*seti.ballMaskROI(:);
        disp('Warning - contrast has been restricted to B(0,rCD/2)');
    end
end

end

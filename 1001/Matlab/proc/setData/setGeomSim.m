%% setGeomSim
% Set geometry and simulation.
%
%% Syntax
%
%   seti = setGeomSim(seti)
%   seti = setGeomSim(seti,dispDepth)
%   seti = setGeomSim(seti,dispDepth,out)
%
%% Description
% |seti = setGeomSim(seti)|
% sets geometry and simulation, i.e. essentially set
% general settings for figures, 
% the grids (computational domain (CD) and region of interest (ROI)), 
% the kernel, 
% the experimental set-up, 
% and the predefined contrast.
%
% |seti = setGeomSim(seti,dispDepth)| does the same, but allows to control
% the depth of displayed messages by |dispDepth|.
%
% |seti = setGeomSim(seti,dispDepth,out)| additionally plots the figures in case
% of |out = 1| and also saves them in case of |out = 2|.
%
%% Example
%
%   init;
%   seti = struct;
%   seti = setGeomSim(seti);
%
%% Input Arguments
%
% * seti    :   structure array, see below in "Output Arguments" and also <setData.html>.
%
% *Optional Input Arguments*
%
% * dispDepth   :   Depth of displayed messages (between 0 and 5).
% * out         :   Depth of output for figures and files: generate no plots (0),
%               generate plots (1), generate plots and save them (2).
%
%% Output Arguments
%
% * seti        :   structure array
%
% Following parameters can be set by user (otherwise default is set):
%
% * seti.tol    :   tolerance for GMRES in <solveLippmannSchwinger.html>
%                   (default: 1E-6)
% * seti.mCD    :   coarse grid size (default: 0). 
%                   If mCD = 0, then two-grid is disabled.
%                   Note that coarse gris is not supported in public version.
%
% * Fields of seti, that are directly processed in setGeomSim, are
% described in this file.
% * Other fields are explained in the corresponding functions
% (<setGrid.html>, <setReshapeVecMat.html>,
% <setGridScale.html>, <setIdImagReal.html>, <setKernel.html>,
% <expSetup.html>, <setContrast.html>).
% An assignment of the fields to the functions is given in <setData.html>.
% 
% *Subfunction: setFigureSettings*
%
% Following parameters can be set by user (otherwise default is set):
%
% * seti.plotFreq       :   (default: 1) 
%                           0: do not plot; 
%                           n > 0: plot each n-th outer iteration.
% * seti.usecbarlim     :   (default: 1) 
%                           0: colorbar limits are set automatically; 
%                           1: set the limits manually (for real and imag colorbars in 2D)
% * seti.cbarlim        :   (default: [-0.2 1.4]) vector [cbarmin,cbarmax].
% * seti.plotPublish    :   (default: 0) 
%                           0 or 1, save figures designed for publications,
%                           i.e. no title, no axis, bigger font size for colorbar.,
%                           see <plot2DstylePublish.html> and <plot3DstylePublish.html>.
% * seti.pubFontSize    :   (default: 20) font size for publishing figures.
% * seti.plotVisible    :   (default: 'off') 
%                           figure visibility by MATLAB figure property
%                           "Visible", see [1].
% * seti.savepng        :   (default: 1) 0 or 1, save figures as *.png.
% * seti.saveepsc       :   (default: 0) 0 or 1, save figures as colored *.eps.
% * seti.savefig        :   (default: 0) 0 or 1, save figures as *.fig.
% * seti.plotFreqiPda   :   (default: 0) Frequency of plots inside pda 
%                           A folder is made for every outer iteration.
%                           (This option is not available in public version.)
%
%% References
% * [1] MathWorks. MATLAB documentation.
%   https://de.mathworks.com/help/matlab/ref/figure-properties.html.
%   Accessed: 2016-10-07.
%
%% See Also
%
% * <setData.html>
%
% * <solveLippmannSchwinger.html>
%
% * <setGrid.html>
% * <setReshapeVecMat.html>
% * <setGridScale.html>
% * <setIdImagReal.html>
% * <setKernel.html>
% * <expSetup.html>
% * <setContrast.html>
%
% * <plot2DstylePublish.html>
% * <plot3DstylePublish.html>
%
%% Code: setGeomSim
%
function seti = setGeomSim(seti,varargin)

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

seti = checkfield(seti,'tol',1E-6,dispDepth); % used in solveLippmannSchwinger (as GMRES tolerance)

seti = setFigureSettings(seti,dispDepth);

dispMessage(' - setGrid - ',dispDepth);
seti = setGrid(seti,dispDepth);

dispMessage(' - setReshapeVecMat - ',dispDepth);
seti = setReshapeVecMat(seti);

dispMessage(' - setGridScale - ',dispDepth);
seti = setGridScale(seti,dispDepth);

dispMessage(' - setIdImagReal - ',dispDepth);
seti = setIdImagReal(seti);

dispMessage(' - setKernel - ',dispDepth);
seti = setKernel(seti,dispDepth);

dispMessage(' - expSetup - ',dispDepth);
seti = expSetup(seti,dispDepth,out); 

dispMessage(' - setContrast - ',dispDepth);
seti = setContrast(seti,dispDepth,out);
if dispDepth >= 1
    fprintf('   qROIexact: real [%+3.2f,%+3.2f], imag [%+3.2f,%+3.2f].\n',...
        min(real(seti.qROIexact)),max(real(seti.qROIexact)),min(imag(seti.qROIexact)),max(imag(seti.qROIexact)));
end

seti = checkfield(seti,'mCD',0,dispDepth);
if seti.mCD ~= 0 % no coarse grid is used (default)
    dispMessage(' ',dispDepth);
    dispMessage(' - set coarse grid - ',dispDepth);
    dispMessage('Warning: changes in seti have to be copied to seti.setiM from now on!) --',dispDepth);
    dispMessage('---- setCoarse start ------------------------------',dispDepth);
    seti = setCoarse(seti);
    dispMessage('---- setCoarse end --------------------------------',dispDepth);
    dispMessage(' ',dispDepth);
end

end

%% Code: subfunction: dispMessage
%
function dispMessage(string,dispDepth)
if dispDepth >= 1
    disp(string)
end
end

%% Code: subfunction: setFigureSettings
% Set general settings for figures.
%
function seti  = setFigureSettings(seti,varargin)

if nargin == 2
    dispDepth = varargin{1};
else
    dispDepth = 0;
end

seti = checkfield(seti,'plotFreq',1,dispDepth); % of outer iteration...
seti = checkfield(seti,'usecbarlim',1,dispDepth);

if seti.usecbarlim == 1
    seti = checkfield(seti,'cbarlim',[-0.2,1.4],dispDepth); % for contrasts
end

seti = checkfield(seti,'plotPublish',0,dispDepth);
seti = checkfield(seti,'pubFontSize',20,dispDepth);

seti = checkfield(seti,'plotVisible','off',dispDepth); % on or off (default)

% Image formats
seti = checkfield(seti,'savepng',1,dispDepth);
seti = checkfield(seti,'saveepsc',0,dispDepth);
seti = checkfield(seti,'savefig',0,dispDepth);

% Frequency of plots inside pda (not in public version... is experimentally...)
% (is used in pda.m to call plotInsidePda.m)
seti = checkfield(seti,'plotFreqiPda',0,dispDepth);

end

%% setGrid
% Essentially set grids on ROI and CD.
%
%% Syntax
%
%   seti = setGrid(seti)
%   seti = setGrid(seti,dispDepth)
%
%% Description
%
% |seti = setGrid(seti)| essentially sets grids 
% on ROI and CD in |seti.grid| and |seti.gridROI|
% for |seti.dim| dimensions, 
% |seti.nCD| discretization points for each dimension of CD, which
% has a size of |[-seti.rCD,seti.rCD)^seti.dim|.
%
% |seti = setGrid(seti,dispDepth)| does the same, but allows to control the
% depth of displayed messages by |dispDepth|.
%
% ROI is the biggest square in 2D (and cube in 3D)
% in the mathematical sensible region (open ball with radius rCD/2).
%
%% Example
%
% *Example 1*
%
%   init;
%   seti.dim = 2;           % dimension 2
%   seti.rCD = 0.2;         % size of computational domain [-rCD,rCD)^dim
%   seti.nCD = 256;         % number of discretization points for each dimension of CD
%   seti = setGrid(seti);
%
% *Example 2*
%
%   init;
%   seti.dim = 2;           % dimension 2
%   seti.rCD = 0.2;         % size of computational domain [-rCD,rCD)^dim
%   seti.nCD = 256;         % number of discretization points for each dimension of CD
%   seti = setGrid(seti,1);
%
%% Input Arguments
%
% * seti    :   structure array
%
% If input arguments are not set, default values will be used.
%
% * seti.dim     :  dimension of the problem (2 or 3) (default: 2)
% * seti.nCD     :  number of discretization points for each dimension
%                   of computational domain (CD) (in samples)
%                   (default: 256).
% * seti.rCD     :  Size of computational domain [-rCD,rCD)^dim
%                   (default: 0.2) (in meters)
%
% *Optional Input Argument*
%
% * dispDepth   : Depth of displayed messages (0: no, 1 or greater: yes).
%
%% Output Arguments
%
% * seti    :   structure array
%
% * seti.h          :   length of the infinitesimal element of CD.
% * seti.dV         :   area/volume of the infinitesimal element (pixel/voxel) of CD
% * seti.grid       :   grid of computational domain (CD)
%                       (seti.dim x seti.nCD^seti.dim)
% * seti.ballMask   :   mask (logical matrix of size seti.nCD x seti.nCD)
%                       to restrict (later) the contrast in CD 
%                       to the mathematical sensible region.
% * seti.nROI       :   discretization points for each dimension
%                       of region of interest (ROI) (in samples)
% * seti.gridROI    :   grid of region of interest (ROI) 
%                       (matrix of size seti.dim x seti.nROI^seti.dim)
% * seti.ROImask    :   mask (logical matrix of size seti.nCD x seti.nCD)
%                       to restrict (later) the contrast in CD
%                       to the region of interest (ROI).
% * seti.ballMaskROI :  mask (logical matrix of size seti.nROI x seti.nROI)
%                       to restrict (later) the contrast in ROI
%                       to the mathematical sensible region.
%                       (This mask is currently useless, see "More About".)
%
%% More About
%
% * CD: computational domain: square [-rCD,rCD)^dim
% * mathematical sensible region: open ball with radius rCD/2
% * ROI: region of interest: biggest square in 2D (and cube in 3D)
%   in the mathematical sensible region, i.e. for $x = (x_1, x_2)$ in ROI:
%
% $-\frac{r_\mathrm{CD}}{2} \cdot \frac{1}{\sqrt{2}} < x_i < \frac{r_\mathrm{CD}}{2} \cdot \frac{1}{\sqrt{2}}$
%
% The areas are also defined and illustrated in Section 3 in [1].
% (Note that in [1] is a small mistake: The interval must be open to be
% inside the mathematical sensible region $B_R$.) 
% Note that the radius $R$ in [1] is related to the mentioned
% $\texttt{rCD}$ by $\texttt{rCD} = 2 R$.
%
% *The masks...*
%
% * The mask "ballMask" is defined to restrict
%   the contrast in CD to the mathematical sensible region.
%
% * The mask "ballMaskROI" is defined to restrict
%   the contrast in ROI to the mathematical sensible region. 
%   Note that this does only make sense, if the alternative
%   definition of ROI is used, because otherwise ROI is already inside the
%   mathematical sensible region. The alternative definition of ROI is:
%
%   % Alternative:
%   indx = (-seti.rCD/2<x1)&(x1<seti.rCD/2); % ROI with length rCD
%
% * The mask "ROImask" restricts the contrast in CD to the biggest square
%   inside the mathematical sensible region.
%
%
%% References
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
% * <setGeomSim.html>
%
%% Code
function seti = setGrid(seti,varargin)

if nargin == 2
    dispDepth = varargin{1};
else 
    dispDepth = 0;
end

%%
% *Check consistency*

if ~(isfield(seti,'dim') && seti.dim >= 2 && seti.dim <= 3)
    seti.dim = 2;
    setmessage(seti,'dim',dispDepth);
end

seti = checkfield(seti,'rCD',0.2,dispDepth);

seti = checkfield(seti,'nCD',256,dispDepth);

%%
% *Element size and volume*

% Length of the infinitesimal element of CD
seti.h = 2*seti.rCD/seti.nCD;

% Volume of the infinitesimal element of CD
seti.dV = seti.h^seti.dim;

%%
% *Set grid on CD in |seti.grid|*

x1 = (-seti.rCD:seti.h:seti.rCD-seti.h);
% x1: discr. intervall [-rCD,rCD)
% in setKernel [-N/2, N/2) is used analog

if seti.dim == 2
    [X1CD, X2CD] = meshgrid(x1,x1);
    seti.grid = [X1CD(:).'; X2CD(:).'];
elseif seti.dim == 3
    [X1CD,X2CD,X3CD] = meshgrid(x1,x1,x1);
    seti.grid = [X1CD(:).'; X2CD(:).'; X3CD(:).'];
else
    error('Error - please choose seti.dim = 2 or 3.');
end

%%
% *Set ballMask on CD in |seti.ballMask|...*
%
% ... to restrict (later) the contrast to the mathematical sensible region
% (ball with radius rCD/2).

switch seti.dim
    case 2
        seti.ballMask = (hypot(X1CD,X2CD) < seti.rCD/2);
    case 3
        h3 = @(a,b,c) sqrt(abs(a).^2+abs(b).^2+abs(c).^2);
        seti.ballMask = (h3(X1CD,X2CD,X3CD) < seti.rCD/2);
end

seti.ballMask = logical(seti.ballMask); % logical is important(!)

%%
% *Set grid |seti.gridROI| and for restriction |seti.ballMaskROI| on ROI*

% ROI will be the biggest square in 2D (and cube in 3D) in ballMask:
indx = (-seti.rCD/2*1/sqrt(2)<x1)&(x1<seti.rCD/2*1/sqrt(2));

x1ROI = x1(indx);
seti.nROI = length(x1ROI);

if seti.dim == 2
    [X1ROI,X2ROI] = meshgrid(x1ROI,x1ROI);
    seti.gridROI = [X1ROI(:).'; X2ROI(:).'];
    
    seti.ROImask = zeros(seti.nCD,seti.nCD);
    seti.ROImask(indx,indx) = 1;

    seti.ballMaskROI = (hypot(X1ROI,X2ROI) < seti.rCD/2);

elseif seti.dim == 3
    [X1ROI,X2ROI,X3ROI] = meshgrid(x1ROI,x1ROI,x1ROI);
    seti.gridROI = [X1ROI(:).'; X2ROI(:).'; X3ROI(:).'];

    % To reconstruct X1ROI, X2ROI and X3ROI later you could use:
    % X1ROI = reshape(seti.gridROI(1,:),seti.nROI,seti.nROI,seti.nROI);
    % X2ROI = reshape(seti.gridROI(2,:),seti.nROI,seti.nROI,seti.nROI);
    % X3ROI = reshape(seti.gridROI(3,:),seti.nROI,seti.nROI,seti.nROI);

    seti.ROImask = zeros(seti.nCD,seti.nCD,seti.nCD);
    seti.ROImask(indx,indx,indx) = 1;

    seti.ballMaskROI = (h3(X1ROI,X2ROI,X3ROI) < seti.rCD/2);
    
else
    error('Error - please choose seti.dim = 2 or 3.'); % checked in checkConsistency
end

seti.ROImask = logical(seti.ROImask); % logical is important(!)
end

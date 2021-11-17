%% expSetup
% Set the experimental set-up as well as incident fields and measurements.
%
%% Syntax
%
%   seti = expSetup(seti)
%   seti = expSetup(seti,dispDepth)
%   seti = expSetup(seti,dispDepth,out)
%
%% Description
%
% |seti = expSetup(seti)| does the following
%
% # Check consistency of input in struct |seti|
%   (otherwise default values are set) (subfunction |expSetupCons|)
% # Define incident and measurement points
%   (i.e. transmitters and receivers positions)
%   (functions <setIncPnts.html> and <setMeasPnts.html>)
%   (If real-world data is loaded this is skipped).
% # Evalute incident fields on ROI (function <setIncField.html>)
% # Set up measurement kernels (such that measurement = k^2 * kernel * solution * voxelVolume)
%   (function <setMeasKer.html>)
% # Plot experimental set-up (function |plotExpSetup|) (in case of |out| is
% 1 or greater.
%
% |seti = expSetup(seti,dispDepth)| does the same but allows to control the
% depth of displayed messages by |dispDepth|.
%
% |seti = expSetup(seti,dispDepth,out)| does the same but allows to plot
% and additionally save figures as well as files controlled by |out|.
%
%% Input Arguments
%
% * |seti|  :   structure array
%
% Several fields in |seti| are required. Because this is an internal
% function we do not list them.
%
% *Optional Input Arguments to differ from default values*
%
% For details of the the following fields look inside subfunction 
% |expSetupCons| below and the functions  <setIncPnts.html>, 
% <setMeasPnts.html>, <pntsGeometry.html> and <pntsGeometry3D.html>:
%
% * seti.incType        :   type of incident field, default: 'pointSource'
%                           ('planeWave' or 'pointSource')
% * seti.measType       :   type of measurement, default: 'nearField'
%                           ('nearField' or 'farField')
% * seti.incNb          :   number of transmitters, default: 35
% * seti.measNb         :   number of receivers, default: 35
% * seti.radSrc         :   radius of circle for transmitters, default: 5
%                           (not necessary in case of incType='planeWave').
% * seti.radMeas        :   radius of circle for receivers, default: 5
%                           (not necessary in case of measType='farField')
% * seti.incPntsType    :   string with type of geometry for transmitters, default: 'circle'
% * seti.measPntsType   :   string with type of geometry for receivers, default: 'circle'
%
% *Optional Input Arguments*
%
% * dispDepth   : Depth of displayed messages (0: no, 1 or greater: yes).
% * out         : Output depth: generate no plots (0),
%                 generate plots (1), generate plots and save them (2).
%
%% Output Arguments
%
% * |seti|  :   structure array
%
% The following fields in |seti| are defined. Look in documentation of 
% corresponding functions and subfunctions (subfunctions are in the Section
% Code).
%
% For details of the the following fields see the functions <setIncPnts.html>, 
% <setMeasPnts.html>, <pntsGeometry.html>:
%
% * seti.incPnts        :   coordinates of transmitters, size seti.dim x seti.incNb
% * seti.dSInc          :   Approximation of the infinitesimal element of closed contour with control points.
% * seti.measPnts       :   coordinates of receivers, size seti.dim x seti.measNb
% * seti.dSMeas         :   Approximation of the infinitesimal element of closed contour with control points.
% 
% Further details of the following parameter are in <setIncField.html>:
%
% * seti.incField       :   incident fields, complex matrix of size 
%                           seti.nROI^seti.dim x seti.incNb
%
% Further details of the following parameter are in <seMeasKer.html>:
%
% * seti.measKer        :   measurement kernel on region of interest (ROI) 
%                           for each receiver 
%                           (complex matrix of size seti.measNb x seti.nROI^seti.dim)
%
%% More About
%
% For a convenience function to plot the transmitters and receivers
% positions see <plotExpSetup.html>.
%
%% See Also
% * <setIncPnts.html>
% * <setMeasPnts.html>
% * <plotExpSetup.html>
% * <setIncField.html>
% * <setMeasKer.html>
%
%% Code: expSetup
%
function seti = expSetup(seti,varargin)

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

seti = expSetupCons(seti,dispDepth);

if ~isfield(seti,'expData') || ( isfield(seti,'expData') && strcmp(seti.expData,'fresnel') )
    % Future work: implementation of incPnts and measPnts type 'manually' in 3D
    seti = setIncPnts(seti,dispDepth);
    seti = setMeasPnts(seti,dispDepth);
end

% Evaluate incident fields on ROI
seti = setIncField(seti);

% and set up measurement kernels seti.measKer
% (such that measurement = k^2 * kernel * solution * voxelVolume, i.e.
%  measurement = uScattRX = FmeasDelta = seti.k^2*seti.measKer*qROI.*(uIncROI+uScattROI)*seti.dV)
seti = setMeasKer(seti);

% Plot experimental set-up (function |plotExpSetup|)
if out >= 1
    plotExpSetup(seti,out);
end

end

%% Code: subfunction: expSetupCons
% Check consistency of input (or set default values)
%
% *Syntax*
%
%   seti = expSetupCons(seti,dispDepth)
%
% *Output Arguments*
%
% * seti.incType    :   type of incident wave
%                       ('planeWave' or 'pointSource')
% * seti.measType   :   type of measurement data ('nearField' or 'farField')
% * seti.incNb      :   number of transmitters
% * seti.measNb     :   number of receivers
% * seti.radSrc     :   radius of sphere transmitters are arranged on
%                       (useless in case of 'planeWave')
%                       Must be *> rCD/2* (is checked and adapted).
% * seti.radMeas    :   radius of sphere receivers are arranged on
%                       (useless in case of farField)
%                       Must be *> rCD/2* (is checked and adapted.)
%
%
% *More About*
%
% * radSrc  :   useless in case of planeWave
% * radMeas :   useless in case of farField
%
% _A stricter restriction than necessary for transmitters' and receivers' positions_
% 
% In the source code we use a stricter restriction for the transmitters' 
% and receivers' positions than necessary.
%
% Transmitters:
% 
% In fact, the singularity of the fundamental solution requires the absence
% of point sources inside the region of interest ROI. This restriction in 
% the continuous formulation of the single-layer potential drops in the 
% discretized version until the transmitter is not nearby a grid point.
%
% In the code we simply require |seti.radSrc > rCD/2| to omit the
% implementation of a nearby condition. Remember that we have chosen ROI as
% biggest square inside the mathematical sensible region, that is the open
% ball with radius rCD/2. Hence, for the current choice of ROI it would 
% be sufficient to check if the
% receivers' points are outside this square, but in case of an adapted ROI
% it is safer to omit the hole mathematical sensible region.
%
% Receivers:
%
% Actually, there is no restriction for receivers' positions until they are 
% not inside the support of the contrast |q|.
%
% In the code we simply require |seti.radMeas > rCD/2|, because actually the true
% contrast is unknown.
%
% *Code*
%
function seti = expSetupCons(seti,dispDepth)

seti = checkfield(seti,'incType','pointSource',dispDepth); % incType: planeWave or pointSource
seti = checkfield(seti,'measType','nearField',dispDepth); % measType: nearField or farField

seti = checkfield(seti,'incNb',35,dispDepth);
seti = checkfield(seti,'measNb',35,dispDepth);

seti = checkfield(seti,'radSrc',5,dispDepth);  % (useless in case of planeWave)
seti = checkfield(seti,'radMeas',5,dispDepth); % (useless in case of farField)

% Check: Is incNb an integer?
if round(seti.incNb)~=seti.incNb || seti.incNb<1
   seti.incNb = round(seti.incNb);
   if seti.incNb<1; seti.incNb = 1; end
   disp('Parameter "incNb" was either non-integer or negative - corrected.')
end

% Check: Is measNb an integer?
if round(seti.measNb)~=seti.measNb || seti.measNb<1
   seti.measNb = round(seti.measNb);
   if seti.measNb<1; seti.measNb = 1; end
   disp('Parameter "measNb" was either non-integer or negative - corrected.')
end

% Check: Is radMeas > rCD/2? (This is required before setting fields!)
if strcmp(seti.measType,'nearField') && isfield(seti,'radMeas')
    if seti.radMeas <= seti.rCD/2
        seti.radMeas = seti.rCD/2;
        disp('Parameter "radMeas" was <= rCD/2. Set radMeas = rCD');
    end
elseif strcmp(seti.measType,'nearField') && ~isfield(seti,'radMeas')
    seti.radMeas = 2*seti.rCD;
    disp('Parameter "radMeas" was not set. Set radMeas = 2*rCD');
end

% Check: Is radSrc > rCD/2? (This is required before setting fields!)
if strcmp(seti.incType,'pointSource') && isfield(seti,'radSrc')
    if seti.radSrc <= seti.rCD/2
        seti.radSrc = seti.rCD/2;
        disp('Parameter "radSrc" was <= rCD/2. Set radSrc = rCD');
    end
elseif strcmp(seti.incType,'pointSource') && (~isfield('radSrc',seti))
    seti.radSrc = 2*seti.rCD;
    disp('Parameter "radSrc" was not set. Set radSrc = 2*rCD');
end

end

%% setMeasPnts
% Set the positions of the receivers.
%
%% Syntax
%
%   seti = setMeasPnts(seti)
%   seti = setMeasPnts(seti,dispDepth)
%
%% Description
% |seti = setMeasPnts(seti)| set the positions |seti.measPnts| of the 
% receivers in dependence of input in struct seti. The number of
% receivers is stored in |seti.measNb|. The approximation of the 
% infinitesimal element of closed contour with control points is stored in
% |seti.dSMeas|.
%
% |seti = setMeasPnts(seti,dispDepth)| does the same but allows to control
% the displayed messages by |dispDepth|.
%
% * If a closed contour does not make sense, |seti.dSMeas| is set to 1.
%
% *In case of 2D* (seti.dim == 2)
%
% * If *near field* is used (seti.measType = 'nearField'), all
% geometries in seti.measPntsType can be used 
% (i.e. 'manually', 'circle', 'square', 'line', 'borehole').
% * If *far field* is used (seti.measType = 'farField'), only geometry
% circle is supported (circle with radius 1 is automatically used in case of far field).
%
% *In case of 3D* (seti.dim == 3)
%
% * Only a ball is supported in case of near field or far field (in
% latter case, the radius is set to 1).
%
%
%% Example
%
% Arrange 5 point sources on a circle with radius 3. Coordinates are stored
% in seti.incPnts.
%
%   init;
%   seti.dim = 2;
%   seti.measType = 'nearField';
%   seti.measPntsType = 'circle';
%   seti.measNb = 5;
%   seti.radMeas = 3;
%   seti = setMeasPnts(seti);
%
%
%% Input Arguments
%
% * seti.dim            :   dimension of the problem (2 or 3)
% * seti.measType       :   Type of measurement field: 'nearField' (default) or 'farField'
%
% *Input Arguments in case of 2D (seti.dim = 2) and near field (seti.incType = 'nearField')*
%
% * seti.measPntsType   :   string with type of geometry:
%                           'manually', 'circle' (default), 'square', 'line', 'borehole'
%
% It follows a *list of type-depending input parameters* to describe the
% details of the geometry. The names are analog to the parameters in 
% <pntsGeometry.html> (mostly with the prefix meas and a big letter,
% exceptions are marked with emphasized).
% Examples are in <pntsGeometry.html>.
% If these parameters are not set, default values are set automatically
% (see inside the code).
%
% _varargin in case of seti.measPntsType = 'manually'_
%
% * seti.measPnts        :   coordinates of points (real matrix of size 2 x Nb).
%
% _varargin in case of seti.measPntsType = 'circle'_
%
% * seti.measNb             :   Number of points (transmitters)
% * _seti.radMeas_          :   radius of circle
%
% _varargin in case of seti.measPntsType = 'square'_
%
% * seti.measNbEdge      :  measNbEdge+1 points are on one edge
%                           (+1 because a point in a corner belongs to two corners).
%                           (Total number of points in square is |4*measNbEdge|.)
% * seti.measEdgeLength  :  Length of one edge.
%
% _varargin in case of seti.measPntsType = 'line'_
%
% * seti.measNbLine      :   Number of points on line.
% * seti.measLineLength  :   Length of line.
% * seti.measLinePos     :   Distance of line right to origin.
%
% _varargin in case of seti.measPntsType = 'borehole'_
%
% Note: Option 'borehole' is only useful if seti.incPntsType is also 'borhole'.
%
% * _seti.boreNbLine_      :   Number of points on line.
% * _seti.boreLineLength_  :   Length of line.
% * _seti.boreLinePos_     :   Distance of line right to origin.
%
% *Input Arguments in case of 3D (seti.dim = 3)*
%
% * seti.incPntsType    :   string with type of geometry:
%                           'manually', 'sphereLatLon', 'sphereFibo' (default).
%
% It follows a *list of type-depending input parameters* to describe the
% details of the geometry. The names are analog to the parameters in 
% <pntsGeometry3D.html> (mostly with the prefix meas and a big letter,
% exceptions are marked with emphasized).
% Examples are in <pntsGeometry3D.html>.
% If these parameters are not set, default values are set automatically
% (see inside the code).
%
% _varargin in case of seti.measPntsType = 'manually'_
%
% * seti.measPnts        :   coordinates of points (real matrix of size 3 x Nb).
%
% _varargin in case of seti.measPntsType = 'sphereLatLon' or 'sphereFibo'_
%
% * seti.measNb             :   Number of points (receivers)
% * _seti.radMeas_          :   radius of sphere
%
% *Optional Input Argument*
%
% * dispDepth   : Depth of displayed messages (0: no, 1 or greater: yes).
%               
%% Output Arguments
%
% * seti.measPnts    :   Coordinates of points (real matrix of size seti.dim x Nb).
% * seti.measNb      :   Number of points.
% * seti.dSMeas      :   Approximation of the infinitesimal element of closed 
%                        contour with control points.
%                        A closed contour does not make sense in the cases
%                        'manually', 'line', and 'borehole', then it is set to 1.
%
% For further information see the corresponding parameters 
% Pnts, Nb, and dS in <pntsGeometry.html> and <pntsGeometry3D.html>.
%
%
%% See Also
%
% * <expSetup.html>
% * <pntsGeometry.html>
% * <pntsGeometry3D.html>
% * <dS2D.html>
% * <dS3D.html>
% * <setIncPnts.html>
%
%
%% Code: function: setIncPnts
%
%%
function seti = setMeasPnts(seti,varargin)

if nargin == 2
    dispDepth = varargin{1};
else
    dispDepth = 0;
end

if seti.dim == 2
    
    seti = checkfield(seti,'measPntsType','circle',dispDepth); % default

    % incPntsType and measPntsType must be the same in case of borehole
    if isfield(seti,'incPntsType') && strcmp(seti.incPntsType,'borehole')
        seti.measPntsType = 'borehole';
        setmessage(seti,'measPntsType');
    end

    if strcmp(seti.measType,'farField')
        % currently only circle case... makes sense...
        [seti.measPnts,seti.measNb,seti.dSMeas] = pntsGeometry(seti,'circle',seti.measNb,1); % radMeas: 1
        
    elseif strcmp(seti.measType,'nearField')
        switch seti.measPntsType
            case 'manually'
                [seti.measPnts,seti.measNb,seti.dSMeas] = pntsGeometry(seti,'manually',seti.measPnts);
            case 'circle'
                [seti.measPnts,seti.measNb,seti.dSMeas] = pntsGeometry(seti,'circle',seti.measNb,seti.radMeas);
            case 'square'
                seti = checkfield(seti,'measNbEdge',9);
                seti = checkfield(seti,'measEdgeLength',7);
                [seti.measPnts,seti.measNb,seti.dSMeas] = pntsGeometry(seti,'square',seti.measNbEdge,seti.measEdgeLength);
            case 'line'
                seti = checkfield(seti,'measNbLine',35);
                seti = checkfield(seti,'measLineLength',10);
                seti = checkfield(seti,'measLinePos',-5);
                [seti.measPnts,seti.measNb,seti.dSMeas] = pntsGeometry(seti,'line',seti.measNbLine,seti.measLineLength,seti.measLinePos);
            case 'borehole'
                % checkfield is in setIncPnts: use the same for inc and meas Pnts
                [seti.measPnts,seti.measNb,seti.dSMeas] = pntsGeometry(seti,'borehole',seti.boreNbLine,seti.boreLineLength,seti.boreLinePos,'meas');
            otherwise
                error('No valid measPntsType.')
        end

    else
        error('measType must be farField or nearField.')
    end

elseif seti.dim == 3
    
    seti = checkfield(seti,'measPntsType','sphereFibo',dispDepth); % default
    % seti.measPntsMode = 'scattFullSphere'; % old for sphereLatLon
    % sphereLatLon was default until August 2017.
    
    if strcmp(seti.measType,'farField')
        % currently only type sphere (longiute-latitude or Fibonacci)
        switch seti.measPntsType
            case 'sphereLatLon'
                [seti.measPnts,seti.measNb,seti.dSMeas] = pntsGeometry3D(seti,'sphereLatLon',seti.measNb,1,dispDepth);
            case 'sphereFibo'
                [seti.measPnts,seti.measNb,seti.dSMeas] = pntsGeometry3D(seti,'sphereFibo',seti.measNb,1,dispDepth);
            otherwise
                error('No valid measPntsType.')
        end
    
    elseif strcmp(seti.measType,'nearField')

        switch seti.measPntsType
            case 'manually'
                [seti.measPnts,seti.measNb,seti.dSMeas] = pntsGeometry3D(seti,'manually',dispDepth);
            case 'sphereLatLon'
                [seti.measPnts,seti.measNb,seti.dSMeas] = pntsGeometry3D(seti,'sphereLatLon',seti.measNb,seti.radMeas,dispDepth);
            case 'sphereFibo'
                [seti.measPnts,seti.measNb,seti.dSMeas] = pntsGeometry3D(seti,'sphereFibo',seti.measNb,seti.radMeas,dispDepth);
            otherwise
                error('No valid incPntsType.')
        end
    
    end

else
    disp('   Error: seti.dim is not equal to two or three.')
end
end
